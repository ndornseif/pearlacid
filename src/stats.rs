// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of methods for statistical analysis.

//TODO:
// Interesing tests:
// - Seed output difference
// - Runs Test (Wald-Wolfowitz)
// - Birthday spacings test

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

/// Generate 'sample size' u64s using the supplied rng.
///     -> generates 'sample_size' * 8 bytes.
/// Return chi squared value of the byte distribution.
pub fn byte_distribution_test(test_rng: &mut impl RNG, sample_size: usize) -> (f64, f64) {
    let mut counts: [usize; 256] = [0; 256];
    for _ in 0..sample_size {
        let sample = test_rng.next().to_le_bytes();
        for by in sample {
            counts[by as usize] += 1;
        }
    }
    let expected: f64 = (sample_size as f64 * 8.0) / 256.0;
    let mut chi_squared: f64 = 0.0;
    for value in counts {
        chi_squared += (value as f64 - expected).powi(2) / expected;
    }
    let p = 1.0 - chi_squared_p_value(255, chi_squared);
    (chi_squared, p)
}

/// Generate 'sample size' u64s using rand crate as a reference.
///     -> generates 'sample_size' * 8 bytes.
/// Return chi squared value and p of the byte distribution.
pub fn byte_distribution_test_reference(sample_size: usize) -> (f64, f64) {
    let mut counts: [usize; 256] = [0; 256];
    for _ in 0..sample_size {
        let sample = rand::random::<u64>().to_le_bytes();
        for by in sample {
            counts[by as usize] += 1;
        }
    }
    let expected: f64 = (sample_size as f64 * 8.0) / 256.0;
    let mut chi_squared: f64 = 0.0;
    for value in counts {
        chi_squared += (value as f64 - expected).powi(2) / expected;
    }
    let p = 1.0 - chi_squared_p_value(255, chi_squared);
    (chi_squared, p)
}

/// Get p value for given degrees of freedom and chi squared value.
fn chi_squared_p_value(df: u32, chi_squared: f64) -> f64 {
    let chi_squared_dist = ChiSquared::new(df as f64).unwrap();
    chi_squared_dist.cdf(chi_squared)
}

/// Examines the average distance between u64 values with 'zero_count' leading zeroes.
/// This is the reference test using the rand crate.
///     -> generates 'sample_size' * 8 bytes.
/// Returns the average distance.
pub fn leading_zeros_spacing_test_reference(sample_size: usize, zero_count: usize) -> f64 {
    let mut distances: Vec<usize> = vec![];
    let mask: u64 = u64::MAX >> (64 - zero_count);
    let mut current_distance: usize = 0;
    for _ in 0..sample_size {
        let sample = rand::random::<u64>();
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

/// Examines the average distance between u64 values with 'zero_count' leading zeroes.
/// Using the supplied RNG.
///     -> generates 'sample_size' * 8 bytes.
/// Returns the average distance.
pub fn leading_zeros_frequency_test(
    test_rng: &mut impl RNG,
    sample_size: usize,
    zero_count: usize,
) -> f64 {
    let mut distances: Vec<usize> = vec![];
    let mask: u64 = u64::MAX >> (64 - zero_count);
    let mut current_distance: usize = 0;
    for _ in 0..sample_size {
        let sample = test_rng.next();
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
/// This is the reference test using the rand crate.
///     -> generates 'sample_size' * 8 bytes.
/// Returns the cummulative difference, test statistic.
pub fn monobit_test_reference(sample_size: usize) -> (i64, f64) {
    let mut difference: i64 = 0;
    for _ in 0..sample_size {
        let sample = rand::random::<u64>();
        difference += (sample.count_ones() as i64) - 32;
    }
    let p: f64 = statrs::function::erf::erfc(
        (difference.abs() as f64 / f64::sqrt(sample_size as f64 * 64.0)) * utils::INV_ROOT2,
    );
    (difference, p)
}

/// Measures the difference between the number of ones and zeros generated.
/// Using the supplied RNG.
///     -> generates 'sample_size' * 8 bytes.
/// Returns the cummulative difference, test statistic.
pub fn monobit_test(test_rng: &mut impl RNG, sample_size: usize) -> (i64, f64) {
    let mut difference: i64 = 0;
    for _ in 0..sample_size {
        let sample = test_rng.next();
        difference += (sample.count_ones() as i64) - 32;
    }
    let p: f64 = statrs::function::erf::erfc(
        (difference.abs() as f64 / f64::sqrt(sample_size as f64 * 64.0)) * utils::INV_ROOT2,
    );
    (difference, p)
}
