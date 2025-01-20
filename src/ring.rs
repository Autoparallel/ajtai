//! Ring implementation for lattice-based cryptography using the Number Theoretic Transform (NTT).
//!
//! This module provides a type-safe implementation of cyclotomic rings commonly used in
//! lattice-based cryptography, particularly for schemes based on Ring-LWE. The implementation
//! uses compile-time checks to ensure mathematical correctness of parameters.
//!
//! The core type [`CyclotomicRing`] represents elements of the ring Z[X]/(X^D + 1) where:
//! - D is a power of two
//! - T is a divisor of D
//! - The modulus q satisfies q ≡ 1 + 2T (mod 4T)
//!
//! # Example
//! ```ignore
//! use ajtai::ring::{CyclotomicRing, StandardBasis};
//! use ff::PrimeField;
//!
//! #[derive(PrimeField)]
//! #[PrimeFieldModulus = "17"]
//! #[PrimeFieldGenerator = "3"]
//! #[PrimeFieldReprEndianness = "little"]
//! pub struct MockField([u64; 1]);
//!
//! // Create two polynomials in the ring and multiply them
//! let a = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([1, 1, 0, 0, 0, 0, 0, 0]);
//! let b = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([0, 0, 1, 0, 0, 0, 0, 0]);
//! let c = a * b;
//! ```

use core::ops::{Add, Mul};

use comptime::{unity_power, verify_modulus};

use super::*;

/// Private module to prevent external implementations of the [`Basis`] trait.
mod sealed {
  pub trait Sealed {}
}

/// A marker trait representing the basis in which ring elements are expressed.
///
/// This trait is sealed and cannot be implemented outside this crate.
pub trait Basis: sealed::Sealed {}

/// Represents elements in the standard polynomial basis.
///
/// In this basis, elements are represented as polynomials with coefficients modulo q.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct StandardBasis;

/// Represents elements in the Number Theoretic Transform (NTT) basis.
///
/// In this basis, elements are represented by their evaluations at powers of a primitive root of
/// unity.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct NTTBasis;

impl sealed::Sealed for StandardBasis {}
impl sealed::Sealed for NTTBasis {}

impl Basis for StandardBasis {}
impl Basis for NTTBasis {}

/// A type representing elements of a cyclotomic ring with compile-time parameter validation.
///
/// # Type Parameters
///
/// * `F` - The prime field used for coefficients
/// * `D` - The degree of the cyclotomic polynomial X^D + 1
/// * `T` - A parameter dividing D that determines properties of the modulus
/// * `B` - The basis in which the element is represented ([`StandardBasis`] or [`NTTBasis`])
///
/// # Mathematical Background
///
/// This type represents elements of the ring Z[X]/(X^D + 1) where arithmetic is performed
/// modulo a prime q. The parameters must satisfy:
/// - D is a power of two
/// - T divides D
/// - q ≡ 1 + 2T (mod 4T)
///
/// These conditions ensure the existence of appropriate roots of unity for the NTT.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CyclotomicRing<F: PrimeField, const D: usize, const T: usize, B: Basis> {
  /// The coefficients of the ring element
  pub coefficients: [F; D],
  /// Phantom data to track the basis
  _basis:           PhantomData<B>,
}

impl<F: PrimeField, const D: usize, const T: usize, B: Basis> CyclotomicRing<F, D, T, B>
where
  // D is a power of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:,
  // q = 1 + 2*T (mod 4*T)
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  pub const DEFAULT: Self = Self { coefficients: [F::ZERO; D], _basis: PhantomData::<B> };

  /// Creates a new ring element from an array of field element values.
  ///
  /// # Arguments
  ///
  /// * `values` - Array of length D containing the coefficients of the polynomial
  ///
  /// # Returns
  ///
  /// A new `CyclotomicRing` instance with the given coefficients in the specified basis.
  ///
  /// # Example
  ///
  /// ```ignore
  /// # use ajtai::ring::{CyclotomicRing, StandardBasis};
  /// # use ff::PrimeField;
  /// #
  /// # #[derive(PrimeField)]
  /// # #[PrimeFieldModulus = "17"]
  /// # #[PrimeFieldGenerator = "3"]
  /// # #[PrimeFieldReprEndianness = "little"]
  /// # pub struct ExampleField([u64; 1]);
  /// #
  /// let ring = CyclotomicRing::<ExampleField, 8, 8, StandardBasis>::new([1, 2, 3, 4, 5, 6, 7, 8]);
  /// ```
  pub fn new(values: [u64; D]) -> Self {
    Self { coefficients: values.map(F::from), _basis: PhantomData::<B> }
  }
}

impl<F: PrimeField, const D: usize, const T: usize, B: Basis> From<[F; D]>
  for CyclotomicRing<F, D, T, B>
{
  fn from(coefficients: [F; D]) -> Self {
    CyclotomicRing { coefficients, _basis: PhantomData::<B> }
  }
}

// TODO: I don't think we can always just square to get the omega we want
// TODO: This can be optimized heavily using Cooley-Tukey and other tricks
impl<F: PrimeField, const D: usize, const T: usize> CyclotomicRing<F, D, T, StandardBasis>
where
  // D is a power of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:,
  // q = 1 + 2*T (mod 4*T)
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  /// Converts this element to the NTT basis by evaluating at appropriate powers of a root of unity.
  ///
  /// This implementation uses a naive quadratic-time algorithm. For production use, it should
  /// be replaced with an FFT-based implementation.
  pub fn ntt(self) -> CyclotomicRing<F, D, T, NTTBasis> {
    // TODO: calling `.pow()` is going to be more inefficient than other methods probably. This can
    // likely all be done at compile time.
    let omega = F::MULTIPLICATIVE_GENERATOR.pow([unity_power::<F, D>()]);
    let mut result = [F::ZERO; D];

    #[allow(clippy::needless_range_loop)]
    for i in 0..D {
      let omega_i = omega.pow([i as u64]);

      for j in 0..D {
        let term = self.coefficients[j] * omega_i.pow([j as u64]);
        result[i] += term;
      }
    }

    CyclotomicRing { coefficients: result, _basis: PhantomData }
  }
}

impl<F: PrimeField, const D: usize, const T: usize> CyclotomicRing<F, D, T, NTTBasis>
where
  // D is a power
  // of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:,
  // q = 1 + 2*T (mod 4*T)
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  /// Converts this element from NTT basis back to standard basis.
  ///
  /// This involves evaluating the Lagrange interpolation formula, which is currently
  /// implemented in a naive quadratic-time algorithm.
  pub fn intt(self) -> CyclotomicRing<F, D, T, StandardBasis> {
    // TODO: calling `.pow()` is going to be more inefficient than other methods probably. This can
    // likely all be done at compile time.
    let omega = F::MULTIPLICATIVE_GENERATOR.pow([unity_power::<F, D>()]);
    // TODO: probably don't need to call invert explicitly since that is slow and we can probably
    // const compute the correct thing here
    let omega_inv = omega.invert().unwrap();
    let n_inv = F::from(D as u64).invert().unwrap();

    let mut result = [F::ZERO; D];

    #[allow(clippy::needless_range_loop)]
    for i in 0..D {
      for j in 0..D {
        let power = omega_inv.pow([(i * j) as u64]);
        let term = self.coefficients[j] * power;
        result[i] += term;
      }

      result[i] *= n_inv;
    }

    CyclotomicRing { coefficients: result, _basis: PhantomData }
  }
}

impl<F: PrimeField, const D: usize, const T: usize, B: Basis> Add for CyclotomicRing<F, D, T, B>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  type Output = Self;

  /// Adds two ring elements coefficient-wise.
  fn add(self, rhs: Self) -> Self::Output {
    Self {
      coefficients: core::array::from_fn(|i| self.coefficients[i] + rhs.coefficients[i]),
      _basis:       PhantomData,
    }
  }
}

impl<F: PrimeField, const D: usize, const T: usize> Mul for CyclotomicRing<F, D, T, NTTBasis>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  type Output = Self;

  /// Multiplies two elements in NTT basis coefficient-wise.
  fn mul(self, rhs: Self) -> Self::Output {
    Self {
      coefficients: core::array::from_fn(|i| self.coefficients[i] * rhs.coefficients[i]),
      _basis:       PhantomData,
    }
  }
}

impl<F: PrimeField, const D: usize, const T: usize> Mul for CyclotomicRing<F, D, T, StandardBasis>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  type Output = Self;

  /// Multiplies two elements in standard basis by converting to NTT basis,
  /// multiplying coefficient-wise, and converting back.
  fn mul(self, rhs: Self) -> Self::Output {
    let lhs_ntt = self.ntt();
    let rhs_ntt = rhs.ntt();
    (lhs_ntt * rhs_ntt).intt()
  }
}

#[cfg(test)]
mod tests {

  use super::*;

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_valid_dimensions() {
    let _ = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([1; 8]);
    let _ = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new([1; 16]);
  }

  #[rstest]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::basic_addition(
      [1; 16],
      [2; 16],
      [3; 16]
  )]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::with_zero(
      [0; 16],
      [5; 16],
      [5; 16]
  )]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::wrapping_around_modulus(
      [15; 16],
      [10; 16],
      [8; 16]  // (15+10)%17=8
  )]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::large_numbers(
      [16; 16],
      [16; 16],
      [15; 16]  // (16+16)%17=15
  )]
  fn test_add_dimension_sixteen(
    #[case] a: [u64; 16],
    #[case] b: [u64; 16],
    #[case] expected: [u64; 16],
  ) {
    let ring_a = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new(a);
    let ring_b = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new(b);
    let expected_ring = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new(expected);

    assert_eq!(ring_a + ring_b, expected_ring);
  }

  #[rstest]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::basic_addition(
      [1, 2, 3, 4, 5, 6, 7, 8],
      [2, 3, 4, 5, 6, 7, 8, 9],
      [3, 5, 7, 9, 11, 13, 15, 0]  // Last value wraps: (8+9)%17=0
  )]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::with_zero(
      [0, 0, 0, 0, 0, 0, 0, 0],
      [1, 2, 3, 4, 5, 6, 7, 8],
      [1, 2, 3, 4, 5, 6, 7, 8]
  )]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::wrapping_around_modulus(
      [15, 15, 15, 15, 15, 15, 15, 15],
      [10, 10, 10, 10, 10, 10, 10, 10],
      [8, 8, 8, 8, 8, 8, 8, 8]  // (15+10)%17=8
  )]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  #[case::alternating_pattern(
      [16, 0, 16, 0, 16, 0, 16, 0],
      [0, 16, 0, 16, 0, 16, 0, 16],
      [16, 16, 16, 16, 16, 16, 16, 16]
  )]
  fn test_add_dimension_eight(
    #[case] a: [u64; 8],
    #[case] b: [u64; 8],
    #[case] expected: [u64; 8],
  ) {
    let ring_a = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new(a);
    let ring_b = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new(b);
    let expected_ring = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new(expected);

    assert_eq!(ring_a + ring_b, expected_ring);
  }

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_ntt_simple() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([1, 1, 0, 0, 0, 0, 0, 0]);

    // Compute NTT
    let ntt_result = input.ntt();

    // Let's verify the result manually:
    // For polynomial 1 + X, evaluating at powers of omega (where omega = 3^2 = 9)
    // The evaluation points are: 9^0, 9^1, 9^2, 9^3, 9^4, 9^5, 9^6, 9^7
    // At X = 9^0:        1 + 1 = 2
    // At X = 9^1:        1 + 9 = 10
    // At X = 9^2 ≡ 13:   1 + 13 = 14
    // At X = 9^3 ≡ 15:   1 + 15 = 16
    // At X = 9^4 ≡ -1:   1 + -1 = 0
    // At X = 9^5 ≡ 8:    1 + 8 = 9
    // At X = 9^6 ≡ 4:    1 + 4 = 5
    // At X = 9^7 ≡ 2:    1 + 2 = 3
    let expected = CyclotomicRing::<MockField, 8, 8, NTTBasis>::new([2, 10, 14, 16, 0, 9, 5, 3]);

    assert_eq!(ntt_result, expected);
  }

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_inverse_ntt_simple() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([1, 1, 0, 0, 0, 0, 0, 0]);

    // Compute NTT
    let ntt_result = input.ntt();
    let intt_result = ntt_result.intt();

    assert_eq!(input, intt_result);
  }

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_mul_d_8() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input_1 = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([1, 1, 0, 0, 0, 0, 0, 0]);
    let input_2 = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([0, 0, 1, 0, 0, 0, 0, 0]);

    // Compute NTT
    let ntt_result_1 = input_1.ntt();
    let ntt_result_2 = input_2.ntt();
    let mul_result = ntt_result_1 * ntt_result_2;
    let intt_result = mul_result.intt();

    // Correct answer
    let expected = CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([0, 0, 1, 1, 0, 0, 0, 0]);

    assert_eq!(expected, intt_result);
  }

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_mul_d_16() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input_1 = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new([
      1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ]);
    let input_2 = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new([
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ]);

    let result = input_1 * input_2;

    // Correct answer
    let expected = CyclotomicRing::<MockField, 16, 8, StandardBasis>::new([
      0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ]);

    assert_eq!(expected, result);
  }
}
