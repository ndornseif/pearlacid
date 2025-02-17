// Copyright 2025 N. Dornseif
//
// Dual-licensed under Apache 2.0 and MIT terms.

//! Implementation of various rngs.
//! All implement the RNG interface, some feature additional methods like:
//! seek(delta: usize)

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

/// Steam cipher based, add–rotate–XOR PRNG with non linear step.
/// Allows seeking to any position in the output stream.
/// Available: StreamNLARXu128
pub mod stream_nlarx {
    use super::RNG;
    const INITIAL_STATE: u64 = 0;
    const N_ROUNDS: usize = 6;
    const XOR_CONST: u128 = 0x65dcfc916d8e80c9c3cdd6d59b50c964;

    #[derive(Debug, Copy, Clone)]
    ///
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
            let nrng = StreamNLARXu128 {
                state: (seed as u128) << 64 | INITIAL_STATE as u128,
            };
            nrng
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

    pub struct XORShift128 {
        state: [u32; 4],
    }

    impl RNG for XORShift128 {
        fn new(seed: u64) -> Self {
            XORShift128 {
                state: [
                    (seed & 0xFFFFFFFF) as u32,
                    (seed >> 32) as u32,
                    (seed & 0xFFFFFFFF) as u32,
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
            self.state =  [
                    (seed & 0xFFFFFFFF) as u32,
                    (seed >> 32) as u32,
                    (seed & 0xFFFFFFFF) as u32,
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
    /// The .next() method used three RANDU calls to fill the 64 bit output space,
    /// the .next_u32() method returns the reduced original output space.
    pub struct RANDU {
        state: u32,
    }
    impl RNG for RANDU {
        fn new(seed: u64) -> Self {
            RANDU { state: seed as u32 }
        }

        fn next_u32(&mut self) -> u32 {
            self.state = (self.state * 65539) & 0x7fffffff;
            self.state
        }

        fn next(&mut self) -> u64 {
            let a: u64 = self.next_u32() as u64;
            let b: u64 = self.next_u32() as u64;
            let c: u64 = self.next_u32() as u64;
            (a << 42) | ((b & 0x3fffff) << 20) | (c & 0xfffff)
        }

        fn advance(&mut self, delta: usize) {
            for _ in 0..delta {
                let _ = self.next_u32();
            }
        }

        fn reseed(&mut self, seed: u64) {
            self.state = seed as u32;
        }
    }
}
