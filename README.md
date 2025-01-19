# Ajtai

A no-std implementation of the Ajtai commitment scheme with compile-time parameter validation and target-specific optimizations.

## Overview

This library provides a high-performance implementation of cyclic/negacyclic convolution over finite fields using the Number Theoretic Transform (NTT), primarily targeting the Ajtai commitment scheme and other lattice-based cryptographic constructions.

## Features

- **Minimal Dependencies**: Core implementation is `no-std` compatible
- **Compile-time Parameter Validation**: Automatically validates mathematical correctness of NTT parameters
- **Type-safe Ring Operations**: Distinct types for standard and NTT basis prevent incorrect usage
- **Hardware Optimizations**: 
  - Efficient modular reduction tailored for NTT-friendly primes
  - AVX2 optimizations for x86_64 (coming soon)
  - NEON optimizations for ARM (coming soon)
  - optimizations for WASM (coming soon)

## Usage

Add to your `Cargo.toml`:
```toml
[dependencies]
ajtai = "0.1.0"
```

Basic example:
```rust, ignore
use ajtai::ring::{CyclotomicRing, StandardBasis};
use ff::PrimeField;

[derive(PrimeField)]
[PrimeFieldModulus = "17"]
[PrimeFieldGenerator = "3"]
[PrimeFieldReprEndianness = "little"]
pub struct ExampleField([u64; 1]);

// Create polynomials in Z[X]/(X^8 + 1)
let a = CyclotomicRing::<ExampleField, 8, 8, StandardBasis>::new([1, 1, 0, 0, 0, 0, 0, 0]);
let b = CyclotomicRing::<ExampleField, 8, 8, StandardBasis>::new([0, 0, 1, 0, 0, 0, 0, 0]);

// Multiply using NTT
let c = a * b;
```

## Implementation Details

The library implements fast polynomial multiplication in cyclotomic rings Z[X]/(X^n + 1) using the Number Theoretic Transform (NTT). Key components include:

- Type-safe representations of polynomial bases
- Optimized modular reduction for NTT-friendly primes
- Compile-time parameter validation
- Future support for vectorized implementations

## Performance Goals

- Minimal runtime overhead from safety checks
- Competitive with C/assembly implementations
- Platform-specific optimizations
- Constant-time operations where relevant for cryptographic usage

## Roadmap

- [x] Core ring implementation with NTT
- [x] Compile-time parameter validation
- [ ] AVX2 optimizations
- [ ] NEON optimizations
- [ ] Ajtai commitment scheme implementation
- [ ] Constant-time guarantees
- [ ] Benchmarking suite

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

- [Generating Hard Instances of Lattice Problems](https://dl.acm.org/doi/10.1145/237814.237838)
- [LatticeFold](https://eprint.iacr.org/2024/257.pdf)
- [Speeding Up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography](https://eprint.iacr.org/2016/504.pdf)
