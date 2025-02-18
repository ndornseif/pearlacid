# pearlacid

Rust PRNG toolkit.

## Description
All collection of different PRNGs in rust. Some prexisting designs, some custon.
Also includes statistical analysis tools to check RNG performance.

## RNG trait
All PRNGs implement the `RNG` trait that includes the following methods:

### new(seed: u64) -> Self
Initializes and returns a new RNG with the supplied seed.
If the internal state or support for seeds does not allow for 64 seed bits,
the seed is truncated, this means that not all seeds generate uniqe streams for every rng.
### next_u32(&mut self) -> u32
Generate a u32 and andvance the internal state one step. 
Some RNGs generate less than 32 bits per step. In theses cases the internal state might be advanced more than one step.
See individual RNGs for details.
### next(&mut self) -> u64
Generate a u64 and andvance the internal state one step. 
Some RNGs generate less than 64 bits per step. In theses cases the internal state might be advanced more than one step.
See individual RNGs for details.
### advance(&mut self, delta: usize);
Advance the RNGs internal state by 'delta' steps. 
This takes a similar amount of time to generating 'delta' random numbers.
### reseed(&mut self, seed: u64);
Reset internal state as if initialized with supplied seed.
`r.reseed(x)` is equivalent to `let mut r = RNGCALLHERE::new(x)`.

## RNGs


## License

Licensed under either of

 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.