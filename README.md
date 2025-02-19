# pearlacid

A collection of PRNGs and methods for statistical analysis.

## Description
A collection of various PRNGs implemented in Rust, including both pre-existing designs and custom implementations.
Also includes statistical analysis tools to evaluate RNG performance.

## RNG Trait
All PRNGs implement the `RNG` trait, which includes the following methods:

##### `new(seed: u64) -> Self`
Initializes and returns a new RNG with the provided seed.
If the internal state or seed handling does not support 64-bit seeds, the seed is truncated.
This means that not all seeds generate unique streams for every RNG.
Also note that some generators produce low quality output if the seed is zero.

##### `next_u32(&mut self) -> u32`
Generates a `u32` and advances the internal state by one step.
Some RNGs generate fewer than 32 bits per step, in which case the internal state may advance more than one step.
Refer to individual RNG implementations for details.

##### `next(&mut self) -> u64`
Generates a `u64` and advances the internal state by one step.
Some RNGs generate fewer than 64 bits per step, in which case the internal state may advance more than one step.
Refer to individual RNG implementations for details.

##### `advance(&mut self, delta: usize)`
Advances the RNG's internal state by `delta` steps.
For generators that do not support `seek`, this takes a similar amount of time as generating `delta` random numbers.

##### `reseed(&mut self, seed: u64)`
Resets the internal state as if initialized with the provided seed.
`r.reseed(x)` is equivalent to `let mut r = RNGCALLHERE::new(x)`.

### Additional Methods
These methods are not universally implemented.

##### `seek(&mut self, counter: u64)`
Sets the internal state as if, after initialization, it had been advanced `counter` times.
If implemented, this operation is usually very fast.

##### `next_small(&mut self) -> some uint`
For RNGs where advancing the internal state produces fewer than 32 bits, this method is implemented.
It returns only the number of bits generated in one step.

## RNGs

### stream_nlarx
A stream cipher-based add–rotate–XOR PRNG with a non-linear step.
Allows seeking to any position in the output stream with the `seek` method.

| StreamNLARXu128 |   |
|---|---|
| Speed | 12% |
| Fails Tests | None |
| Bits per Step | 64 |
| State Size | 128 |
| Supports | `seek` |

### xorshift
Based on the well-established xorshift architecture.

| XORShift128 |   |
|---|---|
| Speed | 170% |
| Fails Tests | None |
| Output per Step | 32 bits |
| State Size | 128 bits |
| Supports | |

### lcg
Linear congruential generators.

| Randu |   |
|---|---|
| Speed | 170% |
| Fails Tests | Chi², Spectral |
| Output per Step | 31 bits |
| State Size | 32 bits |
| Supports | `next_small` |



| Mmix |   |
|---|---|
| Speed | 220% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 64 bits |
| Supports | |


| UlsLcg512 |   |
|---|---|
| Speed | 100% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 512 bits |
| Supports | |


| UlsLcg512H |   |
|---|---|
| Speed | 100% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 512 bits |
| Supports | |


| Lehmer64 |   |
|---|---|
| Speed | 250% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 128 bits |
| Supports | |

## License

Licensed under either of the following, at your option:

 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT License ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

