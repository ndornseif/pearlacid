// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of PRNGS and methods for statistical analysis.

pub mod conditioning;

pub mod rng_testing;
pub mod rngs;
pub mod stats;
pub mod testdata;
pub mod utils;

use rng_testing::{test_suite, test_suite_with_seeds};
use rngs::RNG;

fn main() {
    let start = std::time::Instant::now();
    const TEST_SIZE_EXPONENT: usize = 22;
    const TEST_SIZE: usize = 1 << TEST_SIZE_EXPONENT;
    println!("\nTesting Reference");
    let mut r = rngs::ReferenceRand::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting OnlyOnes");
    let mut r = rngs::testgens::OnlyOne::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0]);
    println!("\nTesting OnlyZero");
    let mut r = rngs::testgens::OnlyZero::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0]);
    println!("\nTesting AlternatingBlocks");
    let mut r = rngs::testgens::AlternatingBlocks::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0]);
    println!("\nTesting AlternatingBytes");
    let mut r = rngs::testgens::AlternatingBytes::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0]);
    println!("\nTesting AlternatingBits");
    let mut r = rngs::testgens::AlternatingBits::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0]);
    println!("\nTesting RijndaelStream");
    let mut r = rngs::spn::RijndaelStream::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting Lehmer64");
    let mut r = rngs::lcg::Lehmer64::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting RANDU");
    let mut r = rngs::lcg::Randu::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting MMIX");
    let mut r = rngs::lcg::Mmix::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting UlsLcg512");
    let mut r = rngs::lcg::UlsLcg512::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting UlsLcg512H");
    let mut r = rngs::lcg::UlsLcg512H::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting XORShift128");
    let mut r = rngs::xorshift::XORShift128::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("\nTesting StreamNLARXu128");
    let mut r = rngs::stream_nlarx::StreamNLARXu128::new(0);
    test_suite(&mut r, TEST_SIZE);
    println!("Total runtime: {:?}", start.elapsed());
}
