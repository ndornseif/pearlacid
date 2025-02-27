// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Statistical testing of an RNGs output.

use std::time::Instant;

use crate::{
    rngs::{self, RNG},
    stats, testdata, utils,
};

const P_TOLERANCE: f64 = 0.001;
const MAX_P: f64 = 1.0 - P_TOLERANCE;
const MIN_P: f64 = P_TOLERANCE;
const FAIL_STR: &str = "FAILED!!";
const PASS_STR: &str = "PASSED";

const TEST_SEED_COUNT: usize = 8;

const fn construct_seed_array() -> [u64; TEST_SEED_COUNT + 3] {
    let mut arr = [0u64; TEST_SEED_COUNT + 3];
    // TODO: Factor these out into weak seed test.
    arr[0] = u64::MIN;
    arr[1] = u64::MIN + 1;
    arr[2] = u64::MAX;

    let mut i = 0;
    while i < TEST_SEED_COUNT {
        arr[i + 3] = testdata::rng_test::STATIC_TEST_SEEDS[i];
        i += 1;
    }

    arr
}
pub const TEST_SEEDS: [u64; TEST_SEED_COUNT + 3] = construct_seed_array();

macro_rules! run_test {
    ($test_fn:expr, $test_name:expr, $test_data:expr, $p_values:expr, $test_results:expr) => {{
        let start: Instant = Instant::now();
        let p: f64 = $test_fn(&$test_data);
        $p_values.push(p);
        print!(
            "{:<10}: Time: {}     p: {:.6}     flq: {:<4}",
            $test_name,
            utils::format_elapsed_time(start.elapsed()),
            p,
            fail_quantity(p)
        );
        if (MIN_P..=MAX_P).contains(&p) {
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
fn measure_reference_speed() -> f64 {
    const SPEED_MEASURE_SAMPLE_SIZE: usize = 1 << 26; //  512 MiB
    let mut ref_rng = rngs::ReferenceRand::new(0);
    let (_, speed) = stats::generate_test_data(&mut ref_rng, SPEED_MEASURE_SAMPLE_SIZE);
    speed
}

/// Logarithmic quantity to specify how far away from passing a failed test was.
fn fail_quantity(p: f64) -> u32 {

    if p < P_TOLERANCE {
        if p == 0.0 {
            9999
        }else {
            utils::fast_log2((P_TOLERANCE / p) as u64)  
        }
    } else if p > 1.0 - P_TOLERANCE {
        if 1.0 - p == 0.0 {
            9999
        } else{
            utils::fast_log2((P_TOLERANCE / (1.0 - p)) as u64)
        }
    } else {
        0
    }
}

/// Perform performance tests for supplied RNG.
pub fn test_suite(test_rng: &mut impl RNG, sample_size: usize) {
    test_suite_with_seeds(test_rng, sample_size, &TEST_SEEDS);
}
/// Perform performance tests for supplied RNG.
/// Allows supplying a custom list of seeds for testing.
pub fn test_suite_with_seeds(test_rng: &mut impl RNG, sample_size: usize, seeds: &[u64]) {
    // Index zero is passed, one is failed.
    let mut test_results = [0u32; 2];
    let mut p_values: Vec<f64> = vec![];
    test_rng.reseed(testdata::rng_test::STATIC_TEST_SEEDS[0]);
    let start = std::time::Instant::now();
    let pre_clock: u64 = unsafe { core::arch::x86_64::_rdtsc() };
    let (_, speed) = stats::generate_test_data(test_rng, sample_size);
    let cycle_count: f64 = unsafe { core::arch::x86_64::_rdtsc() - pre_clock } as f64;
    let ref_speed: f64 = measure_reference_speed();
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
            p_values,
            test_results
        );
        run_test!(
            stats::leading_zeros_frequency_test,
            "LZ-Space",
            test_data,
            p_values,
            test_results
        );
        run_test!(
            stats::monobit_test,
            "Mono",
            test_data,
            p_values,
            test_results
        );
        run_test!(stats::runs_test, "Runs", test_data, p_values, test_results);
        run_test!(
            stats::u64_block_bit_frequency_test,
            "Blocks",
            test_data,
            p_values,
            test_results
        );
        run_test!(
            stats::longest_ones_run,
            "MaxOnes",
            test_data,
            p_values,
            test_results
        );
        run_test!(
            stats::matrix_ranks,
            "Matrix",
            test_data,
            p_values,
            test_results
        );
    }
    let total_tests: u32 = test_results.iter().sum();
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
}
