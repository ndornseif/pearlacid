// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Collection of PRNGS and methods for statistical analysis.

pub mod conditioning;

pub mod rng_testing;
pub mod rngs;
pub mod stats;
mod strings;
pub mod testdata;
pub mod utils;

use rng_testing::{test_suite, test_suite_with_seeds};
use rngs::RNG;

fn main() {
    let start = std::time::Instant::now();
    const TEST_SIZE_EXPONENT: usize = 22;
    const TEST_SIZE: usize = 1 << TEST_SIZE_EXPONENT;
    let mut r = rngs::ReferenceRand::new(0);
    test_suite(&mut r, TEST_SIZE, "Reference");
    let mut r = rngs::testgens::OnlyOne::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0], "OnlyOnes");
    let mut r = rngs::testgens::OnlyZero::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0], "OnlyZero");
    let mut r = rngs::testgens::AlternatingBlocks::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0], "AlternatingBlocks");
    let mut r = rngs::testgens::AlternatingBytes::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0], "AlternatingBytes");
    let mut r = rngs::testgens::AlternatingBits::new(0);
    test_suite_with_seeds(&mut r, TEST_SIZE, &[0], "AlternatingBits");
    let mut r = rngs::spn::RijndaelStream::new(0);
    test_suite(&mut r, TEST_SIZE, "RijndaelStream");
    let mut r = rngs::lcg::Lehmer64::new(0);
    test_suite(&mut r, TEST_SIZE, "Lehmer64");
    let mut r = rngs::lcg::Randu::new(0);
    test_suite(&mut r, TEST_SIZE, "RANDU");
    let mut r = rngs::lcg::Mmix::new(0);
    test_suite(&mut r, TEST_SIZE, "MMIX");
    let mut r = rngs::lcg::UlsLcg512::new(0);
    test_suite(&mut r, TEST_SIZE, "UlsLcg512");
    let mut r = rngs::lcg::UlsLcg512H::new(0);
    test_suite(&mut r, TEST_SIZE, "UlsLcg512H");
    let mut r = rngs::xorshift::XORShift128::new(0);
    test_suite(&mut r, TEST_SIZE, "XORShift128");
    let mut r = rngs::stream_nlarx::StreamNLARXu128::new(0);
    test_suite(&mut r, TEST_SIZE, "StreamNLARXu128");
    println!("Full program runtime: {:?}", start.elapsed());
}
