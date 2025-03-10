// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Statistical testing of an RNGs output.

use std::{ops::Mul, time::Duration, time::Instant};

use crate::utils::write_and_print;
use crate::{
    rngs::{self, RNG},
    stats, strings, testdata, utils,
};

const P_LOG_STAT_LIMIT_MARGINAL: f64 = 2.0;
const P_LOG_STAT_LIMIT_FAIL: f64 = 4.0;
/// The fraction of all tests that can be marginal
/// while returning a passed overall result.
const MAX_MARGINAL_FRACTION: f64 = 0.05;
const TEST_SEED_COUNT: usize = 16;

const TEST_F_POINTERS: [fn(&[u64]) -> f64; 7] = [
    stats::byte_distribution_test,
    stats::leading_zeros_frequency_test,
    stats::monobit_test,
    stats::runs_test,
    stats::u64_block_bit_frequency_test,
    stats::longest_ones_run,
    stats::matrix_ranks,
];

#[derive(Debug, Copy, Clone)]
struct TestResult {
    test_id: usize,
    p: f64,
    time_used: Duration,
}

impl TestResult {
    pub fn logstat(&self) -> f64 {
        p_log_stat(self.p)
    }
    pub fn passed(&self) -> bool {
        self.logstat() < P_LOG_STAT_LIMIT_MARGINAL
    }
    pub fn marginal(&self) -> bool {
        (P_LOG_STAT_LIMIT_MARGINAL..=P_LOG_STAT_LIMIT_FAIL).contains(&self.logstat())
    }
    pub fn failed(&self) -> bool {
        self.logstat() > P_LOG_STAT_LIMIT_FAIL
    }
    pub fn format(&self) -> String {
        format!(
            "{:<10}: Time: {}     p: {:.6}     pls: {:.4}   - {}",
            strings::TEST_NAMES[self.test_id],
            utils::format_elapsed_time(self.time_used),
            self.p,
            self.logstat(),
            if self.passed() {
                strings::PASS_STR
            } else if self.marginal() {
                strings::MARGINAL_STR
            } else {
                strings::FAIL_STR
            }
        )
    }
}

/// Get the file path used for saving test results.
fn get_result_file_path(rng_name: &str) -> String {
    let mut strvec: Vec<String> = vec![chrono::Local::now()
        .format("pearlacid-%Y-%m-%dT%H:%M:%S-")
        .to_string()];
    strvec.push(rng_name.to_string());
    strvec.push(".txt".to_string());
    strvec.join("")
}

/// Run a test function located at `TEST_F_POINTERS[test_id]`
/// and return the result and excution time.
fn run_single_test(test_data: &[u64], test_id: usize) -> TestResult {
    let start: Instant = Instant::now();
    let p: f64 = TEST_F_POINTERS[test_id](test_data);
    let time_used: Duration = start.elapsed();
    TestResult {
        test_id,
        p,
        time_used,
    }
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
/// -0.2 * (log2(min(p, 1-p)) - 1) clamped to 9.9999
fn p_log_stat(p: f64) -> f64 {
    (p.min(1.0 - p).log2() - 1.0).mul(-0.2).min(9.9999)
}

/// Measure rng speed over sample size and report in bytes/s and cycles/bytes.
/// Also reports speed relative to reference speed.
fn speed_test(test_rng: &mut impl RNG, sample_size: usize) -> String {
    test_rng.reseed(testdata::rng_test::STATIC_TEST_SEEDS[0]);
    let pre_clock: u64 = unsafe { core::arch::x86_64::_rdtsc() };
    let (_, speed) = stats::generate_test_data(test_rng, sample_size);
    let cycle_count: f64 = unsafe { core::arch::x86_64::_rdtsc() - pre_clock } as f64;
    let ref_speed: f64 = measure_reference_speed(sample_size);
    let rel_speed: f64 = (speed / ref_speed) * 100.0;
    format!(
        "Generated {} test data. (Speed: {}/s  ({:.4}%)) ({} cycles ({:.4} cycles/byte))",
        utils::format_byte_count(sample_size * 8),
        utils::format_byte_count(speed as usize),
        rel_speed,
        cycle_count,
        cycle_count / (sample_size as f64 * 8.0)
    )
}

/// Peform all tests listed in `TEST_F_POINTERS` and add the results to `test_results`.
fn test_single_seed(
    test_rng: &mut impl RNG,
    sample_size: usize,
    seed: u64,
    test_results: &mut Vec<TestResult>,
    result_file_path: &str,
) {
    test_rng.reseed(seed);
    write_and_print(
        format!("Testing for seed: {:#018x}", seed),
        result_file_path,
    );
    let (test_data, _) = stats::generate_test_data(test_rng, sample_size);
    for test_id in 0..TEST_F_POINTERS.len() {
        let rslt = run_single_test(&test_data, test_id);
        write_and_print(rslt.format(), result_file_path);
        test_results.push(rslt);
    }
}

fn weak_seeds_tests(
    test_rng: &mut impl RNG,
    sample_size: usize,
    result_file_path: &str,
) -> Vec<u64> {
    let mut found_weak_seeds: Vec<u64> = vec![];
    for seed in testdata::rng_test::WEAK_SEEDS {
        write_and_print(
            format!("Testing weak seed: {:#018x}", seed),
            result_file_path,
        );
        test_rng.reseed(seed);
        let (test_data, _) = stats::generate_test_data(test_rng, sample_size);
        let mut seed_test_results: Vec<TestResult> = vec![];
        for test_id in 0..TEST_F_POINTERS.len() {
            let rslt = run_single_test(&test_data, test_id);
            write_and_print(rslt.format(), result_file_path);
            seed_test_results.push(rslt);
        }
        for rslt in seed_test_results {
            if rslt.failed() {
                found_weak_seeds.push(seed);
                break;
            }
        }
    }
    found_weak_seeds
}

/// Format a vec of `TestResults` and print a summary of the results.
fn format_test_results_summary(test_results: &Vec<TestResult>) -> String {
    const P_LOG_STAT_BINS: usize = 10;
    let mut p_logstat_bins = [0u32; P_LOG_STAT_BINS];
    let mut failed_tests = 0usize;
    let mut marginal_tests = 0usize;
    for rslt in test_results {
        p_logstat_bins[rslt.logstat().floor() as usize] += 1;
        if rslt.marginal() {
            marginal_tests += 1;
        } else if rslt.failed() {
            failed_tests += 1;
        }
    }
    let logstat_summary: String = p_logstat_bins
        .iter()
        .enumerate()
        .map(|(bin, &value)| {
            if bin == P_LOG_STAT_BINS - 1 {
                format!("{:>2}+ : {:04}", bin, value) // Handle last bin with '+'
            } else {
                format!("{:>2} : {:04}|", bin, value)
            }
        })
        .collect::<Vec<String>>()
        .join("");
    let total_tests: usize = test_results.len();
    format!(
        "P log stats: \n{}\nOverall result: {}          ( {} passed; {} marginal; {} failed; {} total)",
        logstat_summary,
        if failed_tests > 0 || marginal_tests as f64 > MAX_MARGINAL_FRACTION * total_tests as f64{
            strings::FAIL_STR
        } else {
            strings::PASS_STR
        },
        total_tests - failed_tests - marginal_tests,
        marginal_tests,
        failed_tests,
        total_tests
    )
}
/// Perform performance tests for supplied RNG.
pub fn test_suite(test_rng: &mut impl RNG, sample_size: usize, rng_name: &str) {
    test_suite_with_seeds(
        test_rng,
        sample_size,
        &testdata::rng_test::STATIC_TEST_SEEDS[0..TEST_SEED_COUNT],
        rng_name,
        true,
    );
}
/// Perform performance tests for supplied RNG.
/// Allows supplying a custom list of seeds for testing.
pub fn test_suite_with_seeds(
    test_rng: &mut impl RNG,
    sample_size: usize,
    seeds: &[u64],
    rng_name: &str,
    test_weak_seeds: bool,
) {
    let full_start = std::time::Instant::now();
    let result_file_path = get_result_file_path(rng_name);
    utils::write_and_print(format!("\nTesting: {}", rng_name), &result_file_path);
    let mut test_results: Vec<TestResult> = vec![];
    utils::write_and_print(speed_test(test_rng, sample_size), &result_file_path);
    for &seed in seeds.iter() {
        test_single_seed(
            test_rng,
            sample_size,
            seed,
            &mut test_results,
            &result_file_path,
        );
    }
    if test_weak_seeds {
        utils::write_and_print(
            format!(
                "Found weak seeds: {:?}",
                weak_seeds_tests(test_rng, sample_size, &result_file_path)
            ),
            &result_file_path,
        );
    }
    utils::write_and_print(format!("\nSummary for: {}", rng_name), &result_file_path);
    utils::write_and_print(
        format_test_results_summary(&test_results),
        &result_file_path,
    );
    write_and_print(
        format!("Total runtime: {:?}", full_start.elapsed()),
        &result_file_path,
    );
}
