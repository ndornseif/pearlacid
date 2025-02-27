// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of methods for statistical analysis.

//TODO:
// Interesing tests:
// - Seed output difference
// - Birthday spacings test
// - Blocks average hamming distance.

use core::f64;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::{rngs::RNG, utils};

/// Generate 'sample size' u64s using the supplied rng.
///     -> generates 'sample_size' * 8 bytes.
/// Write them out to the supplied file path.
pub fn fill_test_file(
    file_path: &str,
    test_rng: &mut impl RNG,
    sample_size: usize,
) -> std::io::Result<()> {
    let file = File::create(file_path)?;
    let mut writer = BufWriter::new(file);
    for _ in 0..sample_size {
        let sample = test_rng.next().to_le_bytes();
        writer.write_all(&sample)?;
    }
    Ok(())
}

/// Generate a vector of lenght 'sample_size'
/// filled with u64 generated using the supplied RNG.
/// Measures the time taken to generate the specified amount of samples.
/// Returns RNG speed in bytes per second.
pub fn generate_test_data(test_rng: &mut impl RNG, sample_size: usize) -> (Vec<u64>, f64) {
    let mut testdata: Vec<u64> = vec![];
    let start = std::time::Instant::now();
    for _ in 0..sample_size {
        testdata.push(test_rng.next());
    }
    let timer = start.elapsed();
    (
        testdata,
        ((sample_size as f64) * 8.0) / ((timer.as_nanos() as f64) / 1e9),
    )
}

/// Generate a ppm image and fill it with random data from supplied RNG.
pub fn fill_test_image(
    file_path: &str,
    test_rng: &mut impl RNG,
    width: usize,
    height: usize,
) -> std::io::Result<()> {
    let path = Path::new(file_path);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let header = format!("P6 {} {} 255\n", width, height);
    writer.write_all(header.as_bytes())?;
    for _ in 0..height * width / 2 {
        let sample = test_rng.next().to_le_bytes();
        writer.write_all(&sample[0..6])?;
    }
    Ok(())
}

/// Measures the distribution among the bytes.
/// Returns p value based on the chi2 statistic.
pub fn byte_distribution_test(test_data: &[u64]) -> f64 {
    if test_data.is_empty() {
        return 0.0;
    }
    let mut counts: [usize; 256] = [0; 256];
    for block in test_data.iter() {
        let sample = block.to_le_bytes();
        for by in sample {
            counts[by as usize] += 1;
        }
    }
    let expected: f64 = (test_data.len() as f64 * 8.0) / 256.0;
    let mut chi_squared: f64 = 0.0;
    for value in counts {
        chi_squared += (value as f64 - expected).powi(2) / expected;
    }
    if chi_squared == 0.0 {
        return 0.0;
    }
    statrs::function::gamma::gamma_lr(chi_squared / 2.0, 255.0 / 2.0).clamp(0.0, 1.0)
}

/// Examines the average distance between u64 values with 'zero_count' leading zeroes.
/// Returns p value based on the chi2 statistic.
pub fn leading_zeros_frequency_test(test_data: &[u64]) -> f64 {
    const BIN_COUNT: usize = 256;
    const EXPECTED_SAMPLE_COUNT: u64 = 16384;

    if test_data.is_empty() {
        return 0.0;
    }
    // Adjust leading zero threshold so the correct amount of distance are expected.
    let zero_count: u32 = utils::fast_log2(test_data.len() as u64 / EXPECTED_SAMPLE_COUNT).max(1);
    let expected_spacing: usize = 1 << zero_count;
    let max_bin: usize = 4 * expected_spacing;
    let base_p: f64 = 1.0 / expected_spacing as f64;
    let bin_spacing: f64 = max_bin as f64 / BIN_COUNT as f64;

    let geometric_cdf = |x: f64| 1.0 - (1.0 - base_p).powf(x);
    let mut bins: [f64; BIN_COUNT] = [0.0; BIN_COUNT];
    let mut expected: [f64; BIN_COUNT] = [0.0; BIN_COUNT];
    let mask: u64 = u64::MAX >> (64 - zero_count);
    let mut current_distance: usize = 0;

    for &sample in test_data {
        if (sample & mask) == 0 {
            let bin_index = (current_distance as f64 / bin_spacing).floor() as usize;
            bins[bin_index.min(BIN_COUNT - 1)] += 1.0;
            current_distance = 0;
        } else {
            current_distance += 1;
        }
    }

    let total_samples: f64 = bins.iter().sum();
    if total_samples == 0.0 {
        return 0.0;
    }
    for (i, entry) in expected.iter_mut().enumerate() {
        *entry = if i == BIN_COUNT - 1 {
            (1.0 - geometric_cdf(bin_spacing * i as f64)) * total_samples
        } else {
            (geometric_cdf(bin_spacing * (i + 1) as f64) - geometric_cdf(bin_spacing * i as f64))
                * total_samples
        };
    }
    let chi_squared: f64 = bins
        .iter()
        .zip(expected.iter())
        .map(|(bin, exp)| (*bin - exp).powi(2) / exp)
        .sum();
    if chi_squared == 0.0 {
        return 0.0;
    }
    statrs::function::gamma::gamma_lr((BIN_COUNT as f64 - 1.0) / 2.0, chi_squared / 2.0)
        .clamp(0.0, 1.0)
}
/// Measures the difference between the number of ones and zeros generated.
/// NIST Special Publication 800-22 Test 2.1
/// Returns p value
pub fn monobit_test(test_data: &[u64]) -> f64 {
    if test_data.is_empty() {
        return 0.0;
    }
    let mut difference: i64 = 0;
    for sample in test_data.iter() {
        difference += (sample.count_ones() as i64) - 32;
    }
    statrs::function::erf::erfc(
        (difference.abs() as f64 / f64::sqrt(test_data.len() as f64 * 64.0)) * utils::INV_ROOT2,
    )
    .clamp(0.0, 1.0)
}

/// Measures the difference between the number of ones and zeroes in the bitstream.
/// An excess of ones is indicated by a positive value.
pub fn count_excess_ones(test_data: &[u64]) -> f64 {
    if test_data.is_empty() {
        return 0.0;
    }
    let mut difference: i64 = 0;
    for sample in test_data.iter() {
        difference += (sample.count_ones() as i64) - 32;
    }
    difference as f64
}

/// Measures the ratio of ones and zeroes in each u64
/// NIST Special Publication 800-22 Test 2.2
/// Returns p value
pub fn u64_block_bit_frequency_test(test_data: &[u64]) -> f64 {
    if test_data.is_empty() {
        return 0.0;
    }
    let mut chi_squared: f64 = 0.0;
    let expected: f64 = 0.5;
    for sample in test_data.iter() {
        chi_squared += ((sample.count_ones() as f64) / 64.0 - expected).powi(2);
    }
    if chi_squared == 0.0 {
        return 0.0;
    }
    chi_squared *= 4.0 * 64.0;
    statrs::function::gamma::gamma_lr((test_data.len() as f64) / 2.0, chi_squared / 2.0)
        .clamp(0.0, 1.0)
}

/// Meansures the number of unintterupted sequence of ones/zeroes.
/// NIST Special Publication 800-22 Test 2.3
/// Returns p value
pub fn runs_test(test_data: &[u64]) -> f64 {
    if test_data.is_empty() {
        return 0.0;
    }
    let mut runs: f64 = 0.0;
    let mut last_bit = (test_data[0] >> 63) & 1; // Extract the MSB of the first word
    let excess_ones = count_excess_ones(test_data);

    for &sample in test_data.iter() {
        let transitions = sample ^ (sample >> 1);
        runs += transitions.count_ones() as f64;

        let first_bit = sample & 1;
        if first_bit != last_bit {
            runs += 1.0; // Count transition across words
        }

        last_bit = (sample >> 63) & 1; // Store last bit for next iteration
        if last_bit != 0 {
            runs -= 1.0;
        }
    }
    let num_bits: f64 = test_data.len() as f64 * 64.0;
    let ones_ratio: f64 = ((num_bits / 2.0) + excess_ones) / num_bits;
    statrs::function::erf::erfc(
        (runs - (2.0 * ones_ratio * num_bits * (1.0 - ones_ratio))).abs()
            / (2.0 * f64::sqrt(2.0 * num_bits) * ones_ratio * (1.0 - ones_ratio)),
    )
    .clamp(0.0, 1.0)
}

/// Divide stream into 8192-bit (1 kiB, 128*u64)blocks.
/// Discarding excess bits.
/// Save the longest run of ones in the block
/// Produces bad results with test data shorter than 100 kiB.
/// NIST Special Publication 800-22 Test 2.4
/// Returns p value
pub fn longest_ones_run(test_data: &[u64]) -> f64 {
    const BIN_COUNT: usize = 5;
    const PI_TABLE: [f64; BIN_COUNT + 1] = [
        0.1344793662428856,
        0.23272062093019485,
        0.2389770820736885,
        0.17245227843523026,
        0.10381045937538147,
        0.11756019294261932,
    ];
    if test_data.is_empty() {
        return 0.0;
    }
    let mut last_bit = 0;
    let mut current_run = 0;
    // The max_runs values are binned as follows:
    // =<10, 11, 12, 13, 14, >=15.
    let mut bins: [f64; BIN_COUNT + 1] = [0.0; BIN_COUNT + 1];

    for chunk in test_data.chunks_exact(128) {
        let mut longest_run = 0;

        for &sample in chunk {
            let mut value = sample;
            if sample == 0 {
                current_run = 0;
                last_bit = 0;
            }

            while value != 0 {
                let ones = value.trailing_ones();

                if last_bit == 1 {
                    longest_run = longest_run.max(ones + current_run);
                } else {
                    longest_run = longest_run.max(ones);
                }

                current_run = ones;
                if ones == 64 {
                    break;
                }
                value >>= ones + value.trailing_zeros();
            }
            last_bit = sample >> 63;
        }
        if longest_run <= 10 {
            bins[0] += 1.0;
        } else if longest_run >= 15 {
            bins[5] += 1.0;
        } else {
            bins[(longest_run - 10) as usize] += 1.0;
        }
    }
    let mut chi_squared: f64 = 0.0;
    let n: f64 = bins.iter().sum();
    for i in 0..=BIN_COUNT {
        chi_squared += (bins[i] - (n * PI_TABLE[i])).powi(2) / (n * PI_TABLE[i])
    }
    if chi_squared == 0.0 {
        return 0.0;
    }
    statrs::function::gamma::gamma_ur(BIN_COUNT as f64 / 2.0, chi_squared / 2.0).clamp(0.0, 1.0)
}

/// Divides the bitstream into 32x32 bit binary matrices.
/// NIST Special Publication 800-22 Test 2.5
/// Each matrix is 1024 bits (128 bytes, 16 * u64).
/// Determines the rank of each matrix over GF(2)
/// and bins the results into three categories.
/// Determine p-value via the chi2 statistic.
/// Returns p value
pub fn matrix_ranks(test_data: &[u64]) -> f64 {
    // All matrices are square.
    const MATRIX_SIZE: usize = 32;
    // Matrix ranks are binned as follows:
    // Full rank, one less than full rank, any lower rank
    // Expected distributions for 32x32 matrix come from:
    // NIST Special Publication 800-22 Section 3.5
    const EXPECTED_DISTRIBUTION: [f64; 3] = [0.2888, 0.5776, 0.1336];
    if test_data.is_empty() {
        return 0.0;
    }
    let mut matrix_ranks: [f64; 3] = [0.0; 3];
    for chunks in test_data.chunks_exact((MATRIX_SIZE * MATRIX_SIZE) / 64) {
        let mut matrix: [u32; MATRIX_SIZE] = [0; MATRIX_SIZE];
        for (i, &block) in chunks.iter().enumerate() {
            matrix[2 * i] = (block >> 32) as u32;
            matrix[2 * i + 1] = block as u32;
        }
        let rank: usize = utils::rank_binary_matrix(matrix);
        if rank == MATRIX_SIZE {
            matrix_ranks[0] += 1.0;
        } else if rank == MATRIX_SIZE - 1 {
            matrix_ranks[1] += 1.0;
        } else {
            matrix_ranks[2] += 1.0;
        }
    }
    let n: f64 = matrix_ranks.iter().fold(0.0, |acc, x| acc + { *x });
    let mut chi_squared: f64 = 0.0;
    for (i, bin) in matrix_ranks.iter().enumerate() {
        chi_squared += (bin - EXPECTED_DISTRIBUTION[i] * n).powi(2) / (EXPECTED_DISTRIBUTION[i] * n)
    }
    ((-1.0 * chi_squared) / 2.0).exp().clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    // Specified in number of u64 blocks.
    const TEST_DATA_LENGTH: f64 = 512.0;
    const TEST_DATA_BITS: f64 = TEST_DATA_LENGTH * 64.0;
    const DEFAULT_PMAX: f64 = 1.0;
    const DEFAULT_PMIN: f64 = 0.0;
    use super::*;
    use crate::rngs;

    fn rng_test_verification(
        test_rng: &mut impl RNG,
        max_p: f64,
        min_p: f64,
        test_func: fn(&[u64]) -> f64,
    ) {
        let (test_data, _) = generate_test_data(test_rng, TEST_DATA_LENGTH as usize);
        let p = test_func(&test_data);
        assert!(
            (min_p..=max_p).contains(&p),
            "p-value out of range: expected [{}, {}], got {}",
            min_p,
            max_p,
            p
        );
    }

    #[test]
    fn monobit_verification_onlyone() {
        rng_test_verification(
            &mut rngs::testgens::OnlyOne::new(0),
            DEFAULT_PMIN,
            DEFAULT_PMIN,
            monobit_test,
        );
    }

    #[test]
    fn monobit_verification_onlyzero() {
        rng_test_verification(
            &mut rngs::testgens::OnlyZero::new(0),
            DEFAULT_PMIN,
            DEFAULT_PMIN,
            monobit_test,
        );
    }

    #[test]
    fn monobit_verification_alternating_bytes() {
        rng_test_verification(
            &mut rngs::testgens::AlternatingBytes::new(0),
            DEFAULT_PMAX,
            DEFAULT_PMAX,
            monobit_test,
        );
    }
    #[test]
    fn monobit_verification_alternating_bits() {
        rng_test_verification(
            &mut rngs::testgens::AlternatingBits::new(0),
            DEFAULT_PMAX,
            DEFAULT_PMAX,
            monobit_test,
        );
    }
    #[test]
    fn monobit_verification_alternating_blocks() {
        rng_test_verification(
            &mut rngs::testgens::AlternatingBlocks::new(0),
            DEFAULT_PMAX,
            DEFAULT_PMAX,
            monobit_test,
        );
    }
    #[test]
    fn monobit_verification_random() {
        rng_test_verification(&mut rngs::ReferenceRand::new(0), 0.999, 0.001, monobit_test);
    }
}
