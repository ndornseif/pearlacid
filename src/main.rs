// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of PRNGS and methods for statistical analysis.

pub mod conditioning;

pub mod rngs;
pub mod stats;
mod testdata;
pub mod utils;

use rngs::RNG;

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
        let (test_data, speed) = stats::generate_test_data(test_rng, sample_size);
        // Relative speed compared to a reference speed of 3.78 GiB/s
        // The reference speed is the speed the rand crate generator runs
        // on a AMD Ryzen 7 5800X
        let rel_speed: f64 = (speed / (2.66 * 1073741824.0)) * 100.0;
        println!(
            "Generated {} test data in {:?}. (Speed: {}/s  ({:.4}%))",
            utils::format_byte_count(sample_size * 8),
            start.elapsed(),
            utils::format_byte_count(speed as usize),
            rel_speed
        );
        let start = std::time::Instant::now();
        let p = stats::byte_distribution_test(&test_data);
        println!(
            "Bytes: Time: {:?} p: {:.6}",
            start.elapsed(),
            p
        );
        let start = std::time::Instant::now();
        let avg_distance = stats::leading_zeros_frequency_test(&test_data, leading_zeroes);
        println!(
            "LZ-Space: Time: {:?}    Leading zeros: {:.2}   Dist:  Expected: {:.4}    Measured: {:.0}",
            start.elapsed(),
            leading_zeroes,
            1 << leading_zeroes,
            avg_distance
        );
        let start = std::time::Instant::now();
        let p = stats::monobit_test(&test_data);
        println!(
            "Mono: Time: {:?} p: {:.6}",
            start.elapsed(),
            p
        );
        let start = std::time::Instant::now();
        let p = stats::runs_test(&test_data);
        println!(
            "Runs: Time: {:?}  p: {:.6}",
            start.elapsed(),
            p
        );
        let start = std::time::Instant::now();
        let p = stats::u64_block_bit_frequency_test(&test_data);
        println!(
            "Blocks: Time: {:?} p: {:.6}",
            start.elapsed(),
            p
        );
        let start = std::time::Instant::now();
        let p = stats::longest_ones_run(&test_data);
        println!(
            "MaxOnes: Time: {:?} p: {:.6}",
            start.elapsed(),
            p
        );
        let start = std::time::Instant::now();
        let p = stats::matrix_ranks(&test_data);
        println!(
            "Matrix: Time: {:?} p: {:.6}",
            start.elapsed(),
            p
        );
    }
}

fn main() {
    let start = std::time::Instant::now();
    const TEST_SIZE_EXPONENT: usize = 24;
    const RANDOMSEEDS: usize = 2;
    let mut seeds: Vec<u64> = vec![0, 1, u64::MAX];
    for _ in 0..RANDOMSEEDS {
        seeds.push(rand::random::<u64>());
    }
    println!("\nTesting Reference");
    let mut r = rngs::ReferenceRand::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &seeds);
    println!("\nTesting OnlyOnes");
    let mut r = rngs::testgens::OnlyOne::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &[0]);
    println!("\nTesting OnlyZero");
    let mut r = rngs::testgens::OnlyZero::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &[0]);
    println!("\nTesting AlternatingBlocks");
    let mut r = rngs::testgens::AlternatingBlocks::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &[0]);
    println!("\nTesting AlternatingBytes");
    let mut r = rngs::testgens::AlternatingBytes::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &[0]);
    println!("\nTesting AlternatingBits");
    let mut r = rngs::testgens::AlternatingBits::new(0);
    test_suite(&mut r, TEST_SIZE_EXPONENT, &[0]);
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
    println!("Total runtime: {:?}", start.elapsed());
}
