// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Statistical testing of an RNGs output.

use std::{ops::Mul, time::Instant};

use crate::{
    rngs::{self, RNG},
    stats, testdata, utils,
};

const P_LOG_STAT_LIMIT: f64 = 3.0;
const FAIL_STR: &str = "FAILED!!";
const PASS_STR: &str = "PASSED";

const TEST_SEED_COUNT: usize = 8;

macro_rules! run_test {
    ($test_fn:expr, $test_name:expr, $test_data:expr, $p_log_stat_values:expr, $test_results:expr) => {{
        let start: Instant = Instant::now();
        let p: f64 = $test_fn(&$test_data);
        let p_log_stat: f64 = p_log_stat(p);
        $p_log_stat_values[p_log_stat.floor() as usize] += 1;
        print!(
            "{:<10}: Time: {}     p: {:.6}     pls: {:.4}",
            $test_name,
            utils::format_elapsed_time(start.elapsed()),
            p,
            p_log_stat
        );
        if p_log_stat < P_LOG_STAT_LIMIT {
            println!("   - {}", PASS_STR);
            $test_results[0] += 1;
        } else {
            println!("   - {}", FAIL_STR);
            $test_results[1] += 1;
        }
    }};
}

/// Measure the speed of the rand crates default RNG.
/// Return in bytes per second.
fn measure_reference_speed(sample_size: usize) -> f64 {
    let mut ref_rng = rngs::ReferenceRand::new(0);
    let (_, speed) = stats::generate_test_data(&mut ref_rng, sample_size);
    speed
}

/// Logarithmic quantity to specify how close to 1.0 or 0.0 a p-value is.
/// Has a range of 0-9.9999.
/// -0.2 * log2(min(p, 1-p)) clamped to 9.9999
fn p_log_stat(p: f64) -> f64 {
    (p.min(1.0 - p).log2()).mul(-0.2).min(9.9999)
}

/// Perform performance tests for supplied RNG.
pub fn test_suite(test_rng: &mut impl RNG, sample_size: usize, rng_name: &str) {
    test_suite_with_seeds(
        test_rng,
        sample_size,
        &testdata::rng_test::STATIC_TEST_SEEDS[0..TEST_SEED_COUNT],
        rng_name,
    );
}
/// Perform performance tests for supplied RNG.
/// Allows supplying a custom list of seeds for testing.
pub fn test_suite_with_seeds(
    test_rng: &mut impl RNG,
    sample_size: usize,
    seeds: &[u64],
    rng_name: &str,
) {
    println!("\nTesting {}", rng_name);
    // Index zero is passed, one is failed.
    let mut test_results = [0u32; 2];
    let mut p_log_stat_values = [0u32; 10];
    test_rng.reseed(testdata::rng_test::STATIC_TEST_SEEDS[0]);
    let start = std::time::Instant::now();
    let pre_clock: u64 = unsafe { core::arch::x86_64::_rdtsc() };
    let (_, speed) = stats::generate_test_data(test_rng, sample_size);
    let cycle_count: f64 = unsafe { core::arch::x86_64::_rdtsc() - pre_clock } as f64;
    let ref_speed: f64 = measure_reference_speed(sample_size);
    let rel_speed: f64 = (speed / ref_speed) * 100.0;
    println!(
        "Generated {} test data in {:?}. (Speed: {}/s  ({:.4}%)) ({} cycles ({:.4} cycles/byte))",
        utils::format_byte_count(sample_size * 8),
        start.elapsed(),
        utils::format_byte_count(speed as usize),
        rel_speed,
        cycle_count,
        cycle_count / (sample_size as f64 * 8.0)
    );
    for &seed in seeds.iter() {
        test_rng.reseed(seed);
        println!("Testing for seed: {:#018x}", seed);
        let (test_data, _) = stats::generate_test_data(test_rng, sample_size);
        run_test!(
            stats::byte_distribution_test,
            "Bytes",
            test_data,
            p_log_stat_values,
            test_results
        );
        run_test!(
            stats::leading_zeros_frequency_test,
            "LZ-Space",
            test_data,
            p_log_stat_values,
            test_results
        );
        run_test!(
            stats::monobit_test,
            "Mono",
            test_data,
            p_log_stat_values,
            test_results
        );
        run_test!(
            stats::runs_test,
            "Runs",
            test_data,
            p_log_stat_values,
            test_results
        );
        run_test!(
            stats::u64_block_bit_frequency_test,
            "Blocks",
            test_data,
            p_log_stat_values,
            test_results
        );
        run_test!(
            stats::longest_ones_run,
            "MaxOnes",
            test_data,
            p_log_stat_values,
            test_results
        );
        run_test!(
            stats::matrix_ranks,
            "Matrix",
            test_data,
            p_log_stat_values,
            test_results
        );
    }
    let total_tests: u32 = test_results.iter().sum();
    println!();
    println!("Summary for: {}", rng_name);
    println!("P log stats:");
    for (val, count) in p_log_stat_values.iter().enumerate() {
        if val == 9 {
            print!(" {:>1}+ : {:<4}", val, count);
        } else {
            print!(" {:>2} : {:<4}", val, count);
        }
    }
    println!();
    println!(
        "Overall result: {}          ( {} / {} passed)",
        if total_tests == test_results[0] {
            PASS_STR
        } else {
            FAIL_STR
        },
        test_results[0],
        total_tests
    );
    println!();
}
