// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Implementation of various rngs.
//! All implement the RNG interface, some feature additional methods like:
//! seek(delta: usize)

use rand::{RngCore, SeedableRng};

/// General trait for PRNGs
pub trait RNG {
    /// Initialize with specified seed.
    fn new(seed: u64) -> Self;
    /// Generate u32 and advance the state one step.
    fn next_u32(&mut self) -> u32;
    /// Generate u64 and advance the state one step.
    /// For generators that dont support full u64 might advance
    /// state more than one step.
    fn next(&mut self) -> u64;
    /// Advance the generator state by the specified amount of steps.
    /// For generators that dont support seek this takes a similar
    /// amount of time to generating (delta) outputs.
    fn advance(&mut self, delta: usize);
    /// Reset to inital state, equivalent to repalcing with ::new(seed).
    fn reseed(&mut self, seed: u64);
}

pub struct RefefenceRand {
    rng: rand::rngs::StdRng,
}

impl RNG for RefefenceRand {
    fn new(seed: u64) -> Self {
        RefefenceRand {
            rng: rand::rngs::StdRng::seed_from_u64(seed),
        }
    }

    fn next_u32(&mut self) -> u32 {
        self.rng.next_u32()
    }

    fn next(&mut self) -> u64 {
        self.rng.next_u64()
    }

    fn advance(&mut self, delta: usize) {
        for _ in 0..delta {
            let _ = self.next();
        }
    }

    fn reseed(&mut self, seed: u64) {
        self.rng = rand::rngs::StdRng::seed_from_u64(seed);
    }
}

/// Steam cipher based, add–rotate–XOR PRNG with non linear step.
/// Allows seeking to any position in the output stream.
pub mod stream_nlarx {
    use super::RNG;
    const INITIAL_STATE: u64 = 0;
    const N_ROUNDS: usize = 6;

    #[derive(Debug, Copy, Clone)]
    pub struct StreamNLARXu128 {
        state: u128,
    }

    fn mix_u128(in_state: u128) -> u128 {
        let mut out_state = in_state;
        for _ in 0..N_ROUNDS {
            out_state = out_state.swap_bytes();
            out_state ^= out_state.rotate_left(17);
            if out_state & 1 != 0 {
                out_state = out_state.wrapping_add(out_state.rotate_left(23));
            } else {
                out_state = out_state.wrapping_add(out_state.rotate_left(41));
            }
            if out_state & 2 != 0 {
                out_state = out_state.wrapping_add(out_state.rotate_left(33));
            } else {
                out_state = out_state.wrapping_add(out_state.rotate_left(17));
            }
        }
        out_state
    }

    impl RNG for StreamNLARXu128 {
        fn new(seed: u64) -> StreamNLARXu128 {
            StreamNLARXu128 {
                state: (seed as u128) << 64 | INITIAL_STATE as u128,
            }
        }
        fn advance(&mut self, delta: usize) {
            self.state = (self.state & 0xffffffffffffffff0000000000000000)
                | (self.state.wrapping_add(delta as u128) & 0x0000000000000000ffffffffffffffff);
        }
        fn next(&mut self) -> u64 {
            self.advance(1);
            mix_u128(self.state) as u64
        }

        fn next_u32(&mut self) -> u32 {
            self.advance(1);
            mix_u128(self.state) as u32
        }

        fn reseed(&mut self, seed: u64) {
            self.state = (seed as u128) << 64 | INITIAL_STATE as u128;
        }
    }
    impl StreamNLARXu128 {
        pub fn seek(&mut self, counter: u64) {
            self.state = (self.state & 0xffffffffffffffff0000000000000000) | counter as u128;
        }
    }
}

// Xorshift PRNGs
pub mod xorshift {
    use super::RNG;
    #[derive(Debug, Copy, Clone)]
    pub struct XORShift128 {
        state: [u32; 4],
    }

    impl RNG for XORShift128 {
        fn new(seed: u64) -> Self {
            XORShift128 {
                state: [
                    seed as u32,
                    (seed >> 32) as u32,
                    seed as u32,
                    (seed >> 32) as u32,
                ],
            }
        }

        fn next_u32(&mut self) -> u32 {
            let mut t: u32 = self.state[3];
            let s: u32 = self.state[0];
            self.state[3] = self.state[2];
            self.state[2] = self.state[1];
            self.state[1] = s;
            t ^= t << 11;
            t ^= t >> 8;
            self.state[0] = t ^ s ^ (s >> 19);
            self.state[0]
        }

        fn next(&mut self) -> u64 {
            let a: u64 = self.next_u32() as u64;
            let b: u64 = self.next_u32() as u64;
            (a << 32) | b
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next_u32();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = [
                seed as u32,
                (seed >> 32) as u32,
                seed as u32,
                (seed >> 32) as u32,
            ];
        }
    }
}

// Linear congruential generators
pub mod lcg {
    use super::RNG;
    /// Ill concieved early LCG, that fails the spectral test badly.
    /// Only has output space of 0-2**31-1.
    /// The .next() method uses three RANDU calls to fill the 64 bit output space,
    /// The .next_u32() method uses two RANDU calls.
    /// the .next_small() method returns the reduced original output space.
    #[derive(Debug, Copy, Clone)]
    pub struct Randu {
        state: u32,
    }

    impl RNG for Randu {
        fn new(seed: u64) -> Self {
            Randu { state: seed as u32 }
        }

        fn next_u32(&mut self) -> u32 {
            let a: u32 = self.next_small();
            let b: u32 = self.next_small();
            a << 15 | (b & 0xffff)
        }

        fn next(&mut self) -> u64 {
            let a: u64 = self.next_small() as u64;
            let b: u64 = self.next_small() as u64;
            let c: u64 = self.next_small() as u64;
            (a << 42) | ((b & 0x3fffff) << 20) | (c & 0xfffff)
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next_small();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = seed as u32;
        }
    }
    impl Randu {
        /// Generate a number in the original reduced output space of 0 to 2**31 - 1.
        fn next_small(&mut self) -> u32 {
            self.state = (self.state * 65539) & 0x7fffffff;
            self.state
        }
    }
    /// Originaly designed by Donald Knuth
    #[derive(Debug, Copy, Clone)]
    pub struct Mmix {
        state: u64,
    }

    impl RNG for Mmix {
        fn new(seed: u64) -> Self {
            Mmix { state: seed }
        }

        fn next_u32(&mut self) -> u32 {
            self.next() as u32
        }

        fn next(&mut self) -> u64 {
            self.state = self.state.wrapping_mul(0x5851f42d4c957f2d);
            self.state = self.state.wrapping_add(0x14057b7ef767814f);
            self.state
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = seed;
        }
    }
    #[derive(Debug, Copy, Clone)]
    pub struct UlsLcg512 {
        state: [u128; 4],
    }

    impl RNG for UlsLcg512 {
        fn new(seed: u64) -> Self {
            UlsLcg512 {
                state: [
                    (!seed as u128) << 64 | !seed as u128,
                    (seed as u128) << 64 | seed as u128,
                    (seed as u128) << 64 | !seed as u128,
                    (!seed as u128) << 64 | seed as u128,
                ],
            }
        }

        fn next_u32(&mut self) -> u32 {
            self.next() as u32
        }

        fn next(&mut self) -> u64 {
            self.state[0] = self.state[0].wrapping_mul(0x59ca1b2888a0a80fc054cd25b1fde311);
            self.state[0] = self.state[0].wrapping_add(0xa53a3854d740d22b4802f2e6ea01e350);
            self.state[1] = self.state[1].wrapping_mul(0xade47f9859546ba094573e7c2194a93c);
            self.state[1] = self.state[1].wrapping_add(0xc77a0728309148b95143795d657a29f2);
            self.state[2] = self.state[2].wrapping_mul(0x85fec39e4833d57dd07f903f191ecfd3);
            self.state[2] = self.state[2].wrapping_add(0x77421f2a59df2305739f337afcad9edb);
            self.state[3] = self.state[3].wrapping_mul(0xcdf30907584f7e1551c0667353108b63);
            self.state[3] = self.state[3].wrapping_add(0x935fec88eaba8c39e94503587c22ce99);
            ((self.state[0] >> 64) as u64)
                ^ ((self.state[1] >> 64) as u64)
                ^ ((self.state[2] >> 64) as u64)
                ^ ((self.state[3] >> 64) as u64)
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = [
                (!seed as u128) << 64 | !seed as u128,
                (seed as u128) << 64 | seed as u128,
                (seed as u128) << 64 | !seed as u128,
                (!seed as u128) << 64 | seed as u128,
            ];
        }
    }
    #[derive(Debug, Copy, Clone)]
    pub struct UlsLcg512H {
        state: [u128; 4],
    }

    impl RNG for UlsLcg512H {
        fn new(seed: u64) -> Self {
            UlsLcg512H {
                state: [
                    (!seed as u128) << 64 | !seed as u128,
                    (seed as u128) << 64 | seed as u128,
                    (seed as u128) << 64 | !seed as u128,
                    (!seed as u128) << 64 | seed as u128,
                ],
            }
        }

        fn next_u32(&mut self) -> u32 {
            self.next() as u32
        }

        fn next(&mut self) -> u64 {
            self.state[0] = self.state[0].wrapping_mul(0xe7513927bf96492135e503ed7f5b837e);
            self.state[0] = self.state[0].wrapping_add(0x126b06c2bfe2dac7725ee66c0e1efe69);
            self.state[1] = self.state[1].wrapping_mul(0x6420fafa38bd7d81fc02e8cbfac57698);
            self.state[1] = self.state[1].wrapping_add(0xd2a884d8ed65a425999f67abfa901eba);
            self.state[2] = self.state[2].wrapping_mul(0x3072f956f9d4a9531efd7c4bd3f684f5);
            self.state[2] = self.state[2].wrapping_add(0x2f18c679c54a581aef3f88efa973d2c9);
            self.state[3] = self.state[3].wrapping_mul(0xa7b5b12dc766a03cfdbaf54bacac8382);
            self.state[3] = self.state[3].wrapping_add(0xb12c82d5df1c4e33fd207ba107b9c620);
            (self.state[0].wrapping_add(
                self.state[1].wrapping_add(self.state[2].wrapping_add(self.state[3])),
            ) >> 64) as u64
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = [
                (!seed as u128) << 64 | !seed as u128,
                (seed as u128) << 64 | seed as u128,
                (seed as u128) << 64 | !seed as u128,
                (!seed as u128) << 64 | seed as u128,
            ];
        }
    }

    #[derive(Debug, Copy, Clone)]
    pub struct Lehmer64 {
        state: u128,
    }
    impl RNG for Lehmer64 {
        fn new(seed: u64) -> Self {
            Lehmer64 {
                state: (seed as u128) << 64 | seed as u128,
            }
        }

        fn next_u32(&mut self) -> u32 {
            self.next() as u32
        }

        fn next(&mut self) -> u64 {
            self.state = self.state.wrapping_mul(0xda942042e4dd58b5);
            (self.state >> 64) as u64
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = (seed as u128) << 64 | seed as u128;
        }
    }
}

/// RNGs based on permutation substitution networks.
pub mod spn {
    use std::arch::x86_64::*;

    use super::RNG;

    /// Implementation is x86 architecture specific.
    /// Will crash if x86 AES instruction set is not available.
    pub struct RijndaelStream {
        counter: u128,
        key: [u8; 16],
    }
    impl RNG for RijndaelStream {
        fn new(seed: u64) -> Self {
            let mut key: [u8; 16] = [0; 16];
            key[0..8].clone_from_slice(&seed.to_le_bytes());
            key[8..16].clone_from_slice(&(!seed).to_le_bytes());
            RijndaelStream { counter: 0, key }
        }

        fn next_u32(&mut self) -> u32 {
            self.next() as u32
        }

        fn next(&mut self) -> u64 {
            #![feature(stdarch)]
            self.advance(1);

            let mut encrypted = [0u8; 16];
            unsafe {
                // Load key and block into SIMD registers
                let key = _mm_loadu_si128(self.key.as_ptr() as *const __m128i);
                let mut block =
                    _mm_loadu_si128(self.counter.to_le_bytes().as_ptr() as *const __m128i);

                for _ in 0..4 {
                    block = _mm_aesenc_si128(block, key);
                }
                _mm_storeu_si128(encrypted.as_mut_ptr() as *mut __m128i, block);
            }
            u128::from_le_bytes(encrypted) as u64
        }

        fn advance(&mut self, delta: usize) {
            self.counter += delta as u128;
        }

        fn reseed(&mut self, seed: u64) {
            let mut key: [u8; 16] = [0; 16];
            key[0..8].clone_from_slice(&seed.to_le_bytes());
            key[8..16].clone_from_slice(&(!seed).to_le_bytes());
            self.key = key;
        }
    }
    impl RijndaelStream {
        pub fn seek(&mut self, counter: u64) {
            self.counter = counter as u128;
        }
    }
}

pub mod testgens {
    use super::RNG;

    pub struct OnlyOne {}
    impl RNG for OnlyOne {
        fn new(_seed: u64) -> Self {
            OnlyOne {}
        }

        fn next_u32(&mut self) -> u32 {
            u32::MAX
        }

        fn next(&mut self) -> u64 {
            u64::MAX
        }

        fn advance(&mut self, _delta: usize) {}

        fn reseed(&mut self, _seed: u64) {}
    }

    pub struct OnlyZero {}
    impl RNG for OnlyZero {
        fn new(_seed: u64) -> Self {
            OnlyZero {}
        }

        fn next_u32(&mut self) -> u32 {
            0
        }

        fn next(&mut self) -> u64 {
            0
        }

        fn advance(&mut self, _delta: usize) {}

        fn reseed(&mut self, _seed: u64) {}
    }

    pub struct AlternatingBlocks {
        state: u64,
    }
    impl RNG for AlternatingBlocks {
        fn new(_seed: u64) -> Self {
            AlternatingBlocks { state: 0 }
        }

        fn next_u32(&mut self) -> u32 {
            self.next() as u32
        }

        fn next(&mut self) -> u64 {
            self.advance(1);
            self.state
        }

        fn advance(&mut self, delta: usize) {
            if delta & 1 == 1 {
                self.state = !self.state;
            }
        }

        fn reseed(&mut self, _seed: u64) {}
    }

    pub struct AlternatingBytes {}
    impl RNG for AlternatingBytes {
        fn new(_seed: u64) -> Self {
            AlternatingBytes {}
        }

        fn next_u32(&mut self) -> u32 {
            0xff00ff00
        }

        fn next(&mut self) -> u64 {
            0xff00ff00ff00ff00
        }

        fn advance(&mut self, _delta: usize) {}

        fn reseed(&mut self, _seed: u64) {}
    }
}
