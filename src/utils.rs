// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Misc utility functions.

use std::{fs::File, io::Write, path::Path, time::Duration};

pub const INV_ROOT2: f64 = 0.7071067811865475;


/// Format a duration to a fixed width.
pub fn format_elapsed_time(duration: Duration) -> String {
    const DECIMAL_DIGITS: usize = 4;
    let round_mul: f64 = 10.0_f64.powi(DECIMAL_DIGITS as i32);
    let secs = duration.as_secs_f64(); 

    if secs >= 1.0 {
        format!("{:<1$} s ", (secs * round_mul).floor() / round_mul, DECIMAL_DIGITS + 4)
    } else if secs >= 1e-3 {
        format!("{:<1$} ms", (secs * 1e3 * round_mul).floor() / round_mul, DECIMAL_DIGITS + 4)
    } else if secs >= 1e-6{
        format!("{:<1$} µs", (secs * 1e6 * round_mul).floor() / round_mul, DECIMAL_DIGITS + 4)
    } else {
        format!("{:<1$} ns", (secs * 1e9 * round_mul).floor() / round_mul, DECIMAL_DIGITS + 4)
    }
}
/// XOR two u64 slices in place.
pub fn xor_in_place(a: &mut [u64], b: &[u64]) {
    for (b1, b2) in a.iter_mut().zip(b.iter()) {
        *b1 ^= *b2;
    }
}

// Calculates fast ceil(log2) of integer.
pub fn fast_log2(in_int: u64) -> u32 {
    if in_int == 0 {
        return 0;
    }
    64 - (in_int - 1).leading_zeros()
}

/// Create 24-bit color .ppm image from byte vec.
/// pixels must contain height * width * 3 bytes.
/// Useful for visually checking for patterns in data.
pub fn create_ppm(
    file_path: &str,
    width: usize,
    height: usize,
    image_data: &[u8],
) -> std::io::Result<()> {
    assert_eq!(image_data.len(), height * width * 3);
    let path = Path::new(file_path);
    let mut file = File::create(path)?;
    let header = format!("P6 {} {} 255\n", width, height);
    file.write_all(header.as_bytes())?;
    file.write_all(image_data)?;
    Ok(())
}

/// Format a number of bytes into a pretty String.
/// e.g. 1048576 is 1 MiB
pub fn format_byte_count(num_bytes: usize) -> String {
    // 2**30 = 1073741824
    if num_bytes >= 1073741824 {
        format!("{:.2} GiB", (num_bytes as f64 / 1073741824.0))
    // 2**20 = 1048576
    } else if num_bytes >= 1048576 {
        format!("{:.2} MiB", (num_bytes as f64 / 1048576.0))
    // 2**10 = 1024
    } else if num_bytes >= 1024 {
        format!("{:.2} KiB", (num_bytes as f64 / 1024.0))
    } else {
        format!("{:.2} B", num_bytes as f64)
    }
}

/// Print a binary matrix represented as list of u32 ints.
pub fn print_matrix(matrix: &[u32]) {
    for &row in matrix {
        println!("{:032b}", row);
    }
}

/// Calculate the rank of a 32x32 binary matrix.
/// Assuming all calculations take place over GF(2).
/// Alternative procedure compared to the one specified
/// in Appendix F of NIST Special Publication 800-22.
/// Speedup of around 2x observed.
pub fn rank_binary_matrix(matrix: [u32; 32]) -> usize {
    // Matrix must be square MAXTRIX_SIZE x MAXTRIX_SIZE
    const MAXTRIX_SIZE: usize = 32;
    let mut mat = matrix;
    let mut rank = 0;

    for col_index in 0..MAXTRIX_SIZE {
        let mask: u32 = 1 << (MAXTRIX_SIZE - 1 - col_index);
        // Find the pivot row in the current rank or below
        if let Some(pivot_row) = (rank..MAXTRIX_SIZE).find(|&r| (mat[r] & mask) != 0) {
            // Swap the pivot row with the current rank row
            mat.swap(rank, pivot_row);
            let pivot_val = mat[rank];

            // Eliminate this column only in rows below the pivot row
            for row in mat.iter_mut().take(MAXTRIX_SIZE).skip(rank + 1) {
                if (*row & mask) != 0 {
                    *row ^= pivot_val;
                }
            }

            rank += 1;
        }
    }

    rank
}

/// Calculate the rank of a 32x32 binary matrix.
/// Procedure from Appendix F of NIST Special Publication 800-22
pub fn rank_binary_matrix_nist(matrix_input: [u32; 32]) -> usize {
    const MAXTRIX_SIZE: usize = 32;
    // Matrix must be square MAXTRIX_SIZE x MAXTRIX_SIZE
    let mut matrix = matrix_input;
    for col in 0..MAXTRIX_SIZE {
        let col_mask: u32 = 1 << (MAXTRIX_SIZE - col - 1);
        // Check if entry at col,col is zero
        if col_mask & matrix[col] == 0 {
            // Search following rows for one at row,col
            for row in col + 1..MAXTRIX_SIZE {
                if col_mask & matrix[row] != 0 {
                    // Swap rows
                    matrix.swap(row, col);
                    break;
                }
            }
        }
        // Check if entry at col,col is now one
        if col_mask & matrix[col] != 0 {
            // Checking for ones in col in following rows
            for row in col + 1..MAXTRIX_SIZE {
                if col_mask & matrix[row] != 0 {
                    matrix[row] ^= matrix[col];
                }
            }
        }
    }
    // Reverse step
    for col in (0..MAXTRIX_SIZE).rev() {
        let col_mask: u32 = 1 << (MAXTRIX_SIZE - col - 1);
        if col_mask & matrix[col] == 0 {
            for row in (0..col).rev() {
                if col_mask & matrix[row] != 0 {
                    // Swap rows
                    matrix.swap(row, col);
                    break;
                }
            }
        }
        if col_mask & matrix[col] != 0 {
            // Checking for ones in col in following rows
            for row in (0..col).rev() {
                if col_mask & matrix[row] != 0 {
                    matrix[row] ^= matrix[col];
                }
            }
        }
    }
    // Count zero rows
    let mut rank: usize = MAXTRIX_SIZE;
    for row in matrix {
        if row == 0 {
            rank -= 1;
        }
    }
    rank
}
#[cfg(test)]
mod tests {
    use crate::testdata;

    use super::*;

    #[test]
    fn binary_matrix_rank_test_nist() {
        for (i, test_matrix) in testdata::matrix_test::TEST_MATRICES.iter().enumerate() {
            println!("Matrix: {}", i);
            assert_eq!(
                rank_binary_matrix_nist(test_matrix.matrix),
                test_matrix.rank
            );
        }
    }
    #[test]
    fn binary_matrix_rank_test() {
        for (i, test_matrix) in testdata::matrix_test::TEST_MATRICES.iter().enumerate() {
            println!("Matrix: {}", i);
            assert_eq!(rank_binary_matrix(test_matrix.matrix), test_matrix.rank);
        }
    }
}
