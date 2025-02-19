// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Misc utility functions.

use std::{fs::File, io::Write, path::Path};

pub const ROOT2: f64 = 1.4142135623730951;
pub const INV_ROOT2: f64 = 0.7071067811865475;

pub fn xor_in_place(a: &mut [u64], b: &[u64]) {
    for (b1, b2) in a.iter_mut().zip(b.iter()) {
        *b1 ^= *b2;
    }
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
    if num_bytes > 1073741824 {
        format!("{:.2} GiB", (num_bytes as f64 / 1073741824.0))
    // 2**20 = 1048576
    } else if num_bytes > 1048576 {
        format!("{:.2} MiB", (num_bytes as f64 / 1048576.0))
    // 2**10 = 1024
    } else if num_bytes > 1024 {
        format!("{:.2} KiB", (num_bytes as f64 / 1024.0))
    } else {
        format!("{:.2} B", num_bytes as f64)
    }
}
