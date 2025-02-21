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
| Speed | 10% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 128 bits |
| Supports | `seek` |

### xorshift
Based on the well-established xorshift architecture.

| XORShift128 |   |
|---|---|
| Speed | 160% |
| Fails Tests | None |
| Output per Step | 32 bits |
| State Size | 128 bits |
| Supports | |


### spn
Substitution–permutation networks.

| RijndaelStream |   |
|---|---|
| Speed | 10% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 128 bits |
| Supports | `seek` |


### lcg
Linear congruential generators.

| Randu |   |
|---|---|
| Speed | 160% |
| Fails Tests | Bytes, Spectral, LZ-Space, Blocks, Runs, Mono |
| Output per Step | 31 bits |
| State Size | 32 bits |
| Supports | `next_small` |



| Mmix |   |
|---|---|
| Speed | 240% |
| Fails Tests | Bytes |
| Output per Step | 64 bits |
| State Size | 64 bits |
| Supports | |


| UlsLcg512 |   |
|---|---|
| Speed | 50% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 512 bits |
| Supports | |


| UlsLcg512H |   |
|---|---|
| Speed | 45% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 512 bits |
| Supports | |


| Lehmer64 |   |
|---|---|
| Speed | 240% |
| Fails Tests | None |
| Output per Step | 64 bits |
| State Size | 128 bits |
| Supports | |

## Tests

### Speed
Measures the absolute speed in bytes/s and the relative speed compared to the reference speed.
The reference speed is the speed at which the rand crate generator runs on an AMD Ryzen 7 5800X.

### Monobit
Shorthand: Mono   
Measures the cumulative difference between the number of ones and zeroes generated over the entire bitstream.
A positive value indicates an excess of ones.
Based on NIST Special Publication 800-22 Test 2.1.

### Block bit frequency
Shorthand: Blocks    
Evaluates the ratio of ones and zeroes in every 64-bit block produced by the generator.
Based on NIST Special Publication 800-22 Test 2.2.
Calculates the p-value based on the chi² statistic.

### Runs
Shorthand: Runs   
Measures the number of uninterrupted 'runs' of ones or zeroes in the bitstream.
Based on NIST Special Publication 800-22 Test 2.3.

### Leading zeroes spacing
Shorthand: LZ-Space    
Measures the average distance between blocks of 64 bits that contain at least n leading zeroes.
n is chosen by default so that the expected value for the number of blocks found is 4096.
TODO: Calculate a proper p-value for this test.
Currently, it only returns the average distance. Comparison of measured distances to the expected exponential distribution is possible.

### Byte frequency
Shorthand: Bytes   
Measures the occurrences of the 256 possible byte values in the output stream, splitting each 64-bit output block into 8 bytes.
Calculates the p-value based on the chi² statistic.

## License

Licensed under either of the following, at your option:

 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT License ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

