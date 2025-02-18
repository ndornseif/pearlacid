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

/// Perform chi squared test for supplied RNG.
/// Uses 0,1 and all ones as fixed seeds.
/// Also generates 'randomseeds' additional random testseeds using the rand crate.
/// The rand crate is also used to run a reference chi squared test.
fn test_suite(test_rng: &mut impl RNG, sample_size: usize, randomseeds: usize) {
    println!(
        "Generating {} per test.",
        utils::format_byte_count(sample_size * 8)
    );
    println!("Running reference RNG byte chi2 test.");
    let start = std::time::Instant::now();
    let (chi_squared, p) = stats::bytes_chi_squared_test_reference(sample_size);
    println!(
        "Time: {:?}    Chi2: {:.2}  p: {:.4}",
        start.elapsed(),
        chi_squared,
        p
    );
    let mut seeds: Vec<u64> = vec![0, 1, 0xffffffffffffffff];
    for _ in 0..randomseeds {
        seeds.push(rand::random::<u64>());
    }
    for seed in seeds.iter() {
        test_rng.reseed(*seed);
        println!("Testing byte chi2 for seed: {:#01x}", seed);
        let start = std::time::Instant::now();
        let (chi_squared, p) = stats::bytes_chi_squared_test(test_rng, sample_size);
        println!(
            "Time: {:?}    Chi2: {:.2}  p: {:.4}",
            start.elapsed(),
            chi_squared,
            p
        );
    }
}

fn main() {
    let sample_exponent: usize = 30;
    let mut r = rngs::lcg::Lehmer64::new(0);
    test_suite(&mut r, 1 << sample_exponent, 3);
    //let _ = stats::fill_test_image("testfiles/d.ppm", &mut r, 7680, 4320);
}
