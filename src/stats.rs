// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of methods for statistical analysis.

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::rngs::RNG;

/// Count the number of each possible byte in a Vec<u8>.
pub fn count_bytes(sample: &Vec<u8>) -> [usize; 256] {
    let mut counts: [usize; 256] = [0; 256];
    for by in sample {
        counts[*by as usize] += 1;
    }
    counts
}

/// Generate 'sample size' u64s using the supplied rng.
///     -> generates 'sample_size' * 8 bytes.
/// Return count of all possible bytes over the entire sample.
pub fn check_byte_distribution(test_rng: &mut impl RNG, sample_size: usize) -> [usize; 256] {
    let mut counts: [usize; 256] = [0; 256];
    for _ in 0..sample_size {
        let sample = test_rng.next().to_le_bytes();
        for by in sample {
            counts[by as usize] += 1;
        }
    }
    counts
}

/// Generate 'sample size' u64s using the rand crate as a reference.
///     -> generates 'sample_size' * 8 bytes.
/// Return count of all possible bytes over the entire sample.
pub fn check_byte_distribution_reference(sample_size: usize) -> [usize; 256] {
    let mut counts: [usize; 256] = [0; 256];
    for _ in 0..sample_size {
        let sample = rand::random::<u64>().to_le_bytes();
        for by in sample {
            counts[by as usize] += 1;
        }
    }
    counts
}

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
pub fn bytes_chi_squared_test(test_rng: &mut impl RNG, sample_size: usize) -> (f64, f64) {
    let distribution: [usize; 256] = check_byte_distribution(test_rng, sample_size);
    let expected: f64 = (sample_size as f64 * 8.0) / 256.0;
    let mut chi_squared: f64 = 0.0;
    for value in distribution {
        chi_squared += (value as f64 - expected).powi(2) / expected;
    }
    let p = 1.0 - chi_squared_p_value(255, chi_squared);
    (chi_squared, p)
}

/// Generate 'sample size' u64s using rand crate as a reference.
///     -> generates 'sample_size' * 8 bytes.
/// Return chi squared value of the byte distribution.
pub fn bytes_chi_squared_test_reference(sample_size: usize) -> (f64, f64) {
    let distribution: [usize; 256] = check_byte_distribution_reference(sample_size);
    let expected: f64 = (sample_size as f64 * 8.0) / 256.0;
    let mut chi_squared: f64 = 0.0;
    for value in distribution {
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
