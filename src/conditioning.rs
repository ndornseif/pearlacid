// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Methods to turn random bits into more constrained data types.

use crate::rngs::RNG;

/// Maps a u64 to the 0..1 range in f64.
/// The destribution is uniform but only uses
/// the lower 52 bits of the u64.
/// Not all possible f64 in the output range are produced by this function.
pub fn u64_to_double(int: u64) -> f64 {
    let return_float = (int & 0x000fffffffffffff) | 0x3ff0000000000000;
    f64::from_bits(return_float) - 1.0
}

/// Generate integer between 'lower' (inclusive) and 'upper' (exclusive).
/// Uses rejection sampling so the number of rng calls required is theoretically unbounded.
pub fn rs_random_int(test_rng: &mut impl RNG, lower: u64, upper: u64) -> u64 {
    assert!(upper > lower);
    let range: u64 = upper - lower;
    let mut rn;
    loop {
        rn = test_rng.next() & (u64::MAX >> (range - 1).leading_zeros());
        if rn < range {
            break;
        }
    }
    rn + lower
}
