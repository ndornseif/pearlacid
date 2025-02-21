// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of methods for statistical analysis.

//TODO:
// Interesing tests:
// - Seed output difference
// - Runs Test (Wald-Wolfowitz)
// - Birthday spacings test

use core::f64;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::{rngs::RNG, utils};
use statrs::distribution::{ChiSquared, ContinuousCDF};

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
pub fn generate_test_data(test_rng: &mut impl RNG, sample_size: usize) -> Vec<u64> {
    let mut testdata: Vec<u64> = vec![];
    for _ in 0..sample_size {
        testdata.push(test_rng.next());
    }
    testdata
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

/// Get p value for given degrees of freedom and chi squared value.
fn chi_squared_p_value(df: u32, chi_squared: f64) -> f64 {
    let chi_squared_dist = ChiSquared::new(df as f64).unwrap();
    chi_squared_dist.cdf(chi_squared)
}

/// Measures the distribution among the bytes.
/// Returns chi2 statistic, p value
pub fn byte_distribution_test(test_data: &[u64]) -> (f64, f64) {
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
    let p = 1.0 - chi_squared_p_value(255, chi_squared);
    (chi_squared, p)
}

/// Examines the average distance between u64 values with 'zero_count' leading zeroes.
/// Returns the average distance.
pub fn leading_zeros_frequency_test(test_data: &[u64], zero_count: usize) -> f64 {
    let mut distances: Vec<usize> = vec![];
    let mask: u64 = u64::MAX >> (64 - zero_count);
    let mut current_distance: usize = 0;
    for sample in test_data.iter() {
        if mask & sample == 0 {
            distances.push(current_distance);
            current_distance = 0;
        } else {
            current_distance += 1;
        }
    }
    let sum: f64 = distances.iter().fold(0.0, |acc, x| acc + *x as f64);
    sum / (distances.len() as f64)
}

/// Measures the difference between the number of ones and zeros generated.
/// NIST Special Publication 800-22 Test 2.1
/// Returns the cummulative difference, p value.
pub fn monobit_test(test_data: &[u64]) -> (i64, f64) {
    let mut difference: i64 = 0;
    for sample in test_data.iter() {
        difference += (sample.count_ones() as i64) - 32;
    }
    let p: f64 = statrs::function::erf::erfc(
        (difference.abs() as f64 / f64::sqrt(test_data.len() as f64 * 64.0)) * utils::INV_ROOT2,
    );
    (difference, p)
}

/// Measures the ratio of ones and zeroes in each u64
/// NIST Special Publication 800-22 Test 2.2
/// Returns chi2 statistic, p value
pub fn u64_block_bit_frequency_test(test_data: &[u64]) -> (f64, f64) {
    let mut chi_squared: f64 = 0.0;
    let expected: f64 = 0.5;
    for sample in test_data.iter() {
        chi_squared += ((sample.count_ones() as f64) / 64.0 - expected).powi(2);
    }
    chi_squared *= 4.0 * 64.0;
    let p: f64 = statrs::function::gamma::checked_gamma_lr(
        (test_data.len() as f64) / 2.0,
        chi_squared / 2.0,
    )
    .unwrap();
    (chi_squared, p)
}

/// Meansures the number of unintterupted sequence of ones/zeroes.
/// NIST Special Publication 800-22 Test 2.3
///     -> generates 'sample_size' * 8 bytes.
/// Returns number of runs, p value
pub fn runs_test(test_data: &[u64], excess_ones: i64) -> (u64, f64) {
    let mut runs: u64 = 0;
    // This sometimes introduces a off by one error
    // If the first bit is a 1.
    // Considerd acceptable error to save additional complexitiy and execution time.
    let mut last_bit = (test_data[0] >> 63) & 1; // Extract the MSB of the first word

    for &sample in test_data.iter() {
        let transitions = sample ^ (sample >> 1); // Transitions within the word
        runs += transitions.count_ones() as u64; // Count them
                                                 // Check transition between words
        let first_bit = sample & 1;
        if first_bit != last_bit {
            runs += 1; // Count transition across words
        }

        last_bit = (sample >> 63) & 1; // Store last bit for next iteration
        if last_bit != 0 {
            runs -= 1;
        }
    }
    let num_bits: f64 = test_data.len() as f64 * 64.0;
    let ones_ratio: f64 = ((num_bits / 2.0) + excess_ones as f64) / num_bits;
    let p: f64 = statrs::function::erf::erfc(
        ((runs as f64) - (2.0 * ones_ratio * num_bits * (1.0 - ones_ratio))).abs()
            / (2.0 * f64::sqrt(2.0 * num_bits) * ones_ratio * (1.0 - ones_ratio)),
    );
    (runs, p)
}

/// Measures the time taken to generate the specified amount of samples.
/// Returns RNG speed in bytes per second.
pub fn speed_test(test_rng: &mut impl RNG, sample_size: usize) -> f64 {
    let start = std::time::Instant::now();
    for _ in 0..sample_size {
        let sample = test_rng.next();
        std::hint::black_box(sample);
    }
    let timer = start.elapsed();
    ((sample_size as f64) * 8.0) / ((timer.as_nanos() as f64) / 1e9)
}
