// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of PRNGS and methods for statistical analysis.

#![allow(dead_code, unused_macros)]

pub mod rngs;
pub mod stats;
mod utils;

use rngs::RNG;

macro_rules! time_it {
    ($tip:literal, $func:stmt) => {
        let start = std::time::Instant::now();
        $func
        println!("{}: {:?}", $tip, start.elapsed());
    };
}

/// Perform performance tests for supplied RNGs.
/// Uses 0,1 and all ones as fixed seeds.
/// Also generates 'randomseeds' additional random test seeds using the rand crate.
/// The rand crate is also used to run a reference test.
fn test_suite(test_rng: &mut impl RNG, sample_exponent: usize, randomseeds: usize) {
    let sample_size: usize = 1 << sample_exponent;
    let leading_zeroes: usize = if sample_exponent > 14 {
        sample_exponent - 14
    } else {
        1
    };

    println!(
        "Generating {} per test.",
        utils::format_byte_count(sample_size * 8)
    );
    println!("Running rand crate reference test.");
    let start = std::time::Instant::now();
    let (chi_squared, p) = stats::bytes_chi_squared_test_reference(sample_size);
    println!(
        "Time: {:?}    Chi2: {:.2}  p: {:.4}",
        start.elapsed(),
        chi_squared,
        p
    );
    let start = std::time::Instant::now();
    let avg_distance = stats::leading_zeros_frequency_test_reference(sample_size, leading_zeroes);
    println!(
        "Time: {:?}    Leading zeros: {:.2}   Dist:  Expected: {:.4}    Measured: {:.0}",
        start.elapsed(),
        leading_zeroes,
        1 << leading_zeroes,
        avg_distance
    );
    let mut seeds: Vec<u64> = vec![0, 1, u64::MAX];
    for _ in 0..randomseeds {
        seeds.push(rand::random::<u64>());
    }
    for seed in seeds.iter() {
        test_rng.reseed(*seed);
        println!("Testing for seed: {:#01x}", seed);
        let start = std::time::Instant::now();
        let (chi_squared, p) = stats::bytes_chi_squared_test(test_rng, sample_size);
        println!(
            "Time: {:?}    Chi2: {:.2}  p: {:.4}",
            start.elapsed(),
            chi_squared,
            p
        );
        let start = std::time::Instant::now();
        let avg_distance =
            stats::leading_zeros_frequency_test(test_rng, sample_size, leading_zeroes);
        println!(
            "Time: {:?}    Leading zeros: {:.2}   Dist:  Expected: {:.4}    Measured: {:.0}",
            start.elapsed(),
            leading_zeroes,
            1 << leading_zeroes,
            avg_distance
        );
    }
}

fn main() {
    const TEST_SIZE_EXPONENT: usize = 32;
    const RANDOMSEEDS: usize = 4;
    let mut r = rngs::lcg::Lehmer64::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
    let mut r = rngs::lcg::Randu::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
    let mut r = rngs::lcg::Mmix::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
    let mut r = rngs::lcg::UlsLcg512::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
    let mut r = rngs::lcg::UlsLcg512H::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
    let mut r = rngs::xorshift::XORShift128::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
    let mut r = rngs::stream_nlarx::StreamNLARXu128::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, RANDOMSEEDS);
}
