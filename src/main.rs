// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of PRNGS and methods for statistical analysis.

#![allow(unused_macros)]

pub mod conditioning;
pub mod rngs;
pub mod stats;
pub mod utils;

use rngs::RNG;

macro_rules! time_it {
    ($tip:literal, $func:stmt) => {
        let start = std::time::Instant::now();
        $func
        println!("{}: {:?}", $tip, start.elapsed());
    };
}

/// Perform performance tests for supplied RNGs.
/// Performs all tests using any of the supplied seeds.
/// Runs: Byte distribution, LZ-Distance, Monobit, U64 blocks.
fn test_suite(test_rng: &mut impl RNG, sample_exponent: usize, seeds: &[u64]) {
    let sample_size: usize = 1 << sample_exponent;
    let leading_zeroes: usize = if sample_exponent > 14 {
        sample_exponent - 14
    } else {
        1
    };
    for seed in seeds.iter() {
        test_rng.reseed(*seed);

        println!("Testing for seed: {:#01x}", seed);
        let start = std::time::Instant::now();
        let testdata = stats::generate_test_data(test_rng, sample_size);
        println!(
            "Generated {} test data in {:?}.",
            utils::format_byte_count(sample_size * 8),
            start.elapsed()
        );
        let start = std::time::Instant::now();
        let (chi_squared, p) = stats::byte_distribution_test(&testdata);
        println!(
            "Bytes: Time: {:?}    Chi2: {:.2}  p: {:.4}",
            start.elapsed(),
            chi_squared,
            p
        );
        let start = std::time::Instant::now();
        let avg_distance = stats::leading_zeros_frequency_test(&testdata, leading_zeroes);
        println!(
            "LZ-Space: Time: {:?}    Leading zeros: {:.2}   Dist:  Expected: {:.4}    Measured: {:.0}",
            start.elapsed(),
            leading_zeroes,
            1 << leading_zeroes,
            avg_distance
        );
        let start = std::time::Instant::now();
        let (bit_difference, p) = stats::monobit_test(&testdata);
        println!(
            "Mono: Time: {:?}    Bit difference: {:.0}   p: {:.4}",
            start.elapsed(),
            bit_difference,
            p
        );
        let start = std::time::Instant::now();
        let (runs, p) = stats::runs_test(&testdata, bit_difference);
        println!(
            "Runs: Time: {:?}    Runs count: {:.0}   p: {:.4}",
            start.elapsed(),
            runs,
            p
        );
        let start = std::time::Instant::now();
        let (chi_squared, p) = stats::u64_block_bit_frequency_test(&testdata);
        println!(
            "Blocks: Time: {:?}    Chi2: {:.0}   p: {:.4}",
            start.elapsed(),
            chi_squared,
            p
        );
        test_rng.reseed(*seed);
        let speed: f64 = stats::speed_test(test_rng, sample_size);
        // Relative speed compared to a reference speed of 3.78 GiB/s
        // The reference speed is the speed the rand crate generator runs
        // on a AMD Ryzen 7 5800X
        let rel_speed: f64 = (speed / (3.78 * ((1 << 30) as f64))) * 100.0;
        println!(
            "Speed: {}/s  ({:.4}%)",
            utils::format_byte_count(speed as usize),
            rel_speed
        );
    }
}

fn main() {
    const TEST_SIZE_EXPONENT: usize = 32;
    const RANDOMSEEDS: usize = 2;
    let mut seeds: Vec<u64> = vec![0, 1, u64::MAX];
    for _ in 0..RANDOMSEEDS {
        seeds.push(rand::random::<u64>());
    }
    println!("\nTesting Reference");
    let mut r = rngs::RefefenceRand::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting RijndaelStream");
    let mut r = rngs::spn::RijndaelStream::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting Lehmer64");
    let mut r = rngs::lcg::Lehmer64::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting RANDU");
    let mut r = rngs::lcg::Randu::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting MMIX");
    let mut r = rngs::lcg::Mmix::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting UlsLcg512");
    let mut r = rngs::lcg::UlsLcg512::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting UlsLcg512H");
    let mut r = rngs::lcg::UlsLcg512H::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting XORShift128");
    let mut r = rngs::xorshift::XORShift128::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting StreamNLARXu128");
    let mut r = rngs::stream_nlarx::StreamNLARXu128::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
}
