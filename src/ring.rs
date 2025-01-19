use core::ops::{Add, Mul};

use comptime::{unity_power, verify_modulus};

use super::*;

// First define our basis tags as zero-sized types
#[derive(Clone, Copy)]
pub struct StandardBasis;
#[derive(Clone, Copy)]
pub struct NTTBasis;

trait Basis {}
impl Basis for StandardBasis {}
impl Basis for NTTBasis {}

// Modify Ring to take a basis parameter
#[derive(Debug, Clone, Copy)]
pub struct Ring<F: PrimeField, const D: usize, const T: usize, B: Basis>
where
  // D is a power of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:,
  // q = 1 + 2*T (mod 4*T)
  [(); verify_modulus::<F, T>() as usize - 1]:, {
  coefficients: [F; D],
  _basis:       core::marker::PhantomData<B>,
}

// TODO: I don't think we can always just square to get the omega we want
// TODO: This can be optimized heavily using Cooley-Tukey and other tricks
impl<F: PrimeField, const D: usize, const T: usize> Ring<F, D, T, StandardBasis>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  pub fn ntt(self) -> Ring<F, D, T, NTTBasis> {
    let omega = F::MULTIPLICATIVE_GENERATOR.pow([unity_power::<F, D>()]);
    let mut result = [F::ZERO; D];

    for i in 0..D {
      let omega_i = omega.pow([i as u64]);

      for j in 0..D {
        let term = self.coefficients[j] * omega_i.pow([j as u64]);
        result[i] += term;
      }
    }

    Ring { coefficients: result, _basis: core::marker::PhantomData }
  }
}

#[cfg(test)]
impl<F: PrimeField, const D: usize, const T: usize> Ring<F, D, T, NTTBasis>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  pub fn intt(self) -> Ring<F, D, T, StandardBasis> {
    // TODO: calling `.pow()` is going to be more inefficient than other methods probably. This can
    // likely all be done at compile time.
    let omega = F::MULTIPLICATIVE_GENERATOR.pow([unity_power::<F, D>()]);
    // TODO: probably don't need to call invert explicitly since that is slow and we can probably
    // const compute the correct thing here
    let omega_inv = omega.invert().unwrap();
    let n_inv = F::from(D as u64).invert().unwrap();

    let mut result = [F::ZERO; D];

    for i in 0..D {
      for j in 0..D {
        let power = omega_inv.pow([(i * j) as u64]);
        let term = self.coefficients[j] * power;
        result[i] += term;
      }

      result[i] *= n_inv;
    }

    Ring { coefficients: result, _basis: core::marker::PhantomData }
  }
}

impl<F: PrimeField, const D: usize, const T: usize, B: Basis> Add for Ring<F, D, T, B>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    Self {
      coefficients: core::array::from_fn(|i| self.coefficients[i] + rhs.coefficients[i]),
      _basis:       core::marker::PhantomData,
    }
  }
}

impl<F: PrimeField, const D: usize, const T: usize> Mul for Ring<F, D, T, NTTBasis>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  type Output = Self;

  fn mul(self, rhs: Self) -> Self::Output {
    Self {
      coefficients: core::array::from_fn(|i| self.coefficients[i] * rhs.coefficients[i]),
      _basis:       core::marker::PhantomData,
    }
  }
}

#[cfg(test)]
mod tests {
  use core::marker::PhantomData;

  use super::*;

  #[test]
  fn test_valid_dimensions() {
    let _: Ring<MockField, 8, 8, StandardBasis> = create_ring([1; 8]);
    let _: Ring<MockField, 16, 8, StandardBasis> = create_ring([1; 16]);
  }

  // Helper function to create Ring instances more easily
  fn create_ring<const D: usize, const T: usize, B: Basis>(
    values: [u64; D],
  ) -> Ring<MockField, D, T, B>
  where
    [(); D.is_power_of_two() as usize - 1]:,
    [(); (D % T == 0) as usize - 1]:,
    [(); verify_modulus::<MockField, T>() as usize - 1]:, {
    Ring { coefficients: values.map(MockField::from), _basis: PhantomData::<B> }
  }
  #[rstest]
  #[case::basic_addition(
      [1; 16],
      [2; 16],
      [3; 16]
  )]
  #[case::with_zero(
      [0; 16],
      [5; 16],
      [5; 16]
  )]
  #[case::wrapping_around_modulus(
      [15; 16],
      [10; 16],
      [8; 16]  // (15+10)%17=8
  )]
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
    let ring_a: Ring<MockField, 16, 8, StandardBasis> = create_ring(a);
    let ring_b: Ring<MockField, 16, 8, StandardBasis> = create_ring(b);
    let expected_ring: Ring<MockField, 16, 8, StandardBasis> = create_ring(expected);

    assert_eq!(
      (ring_a + ring_b).coefficients.map(|f| f.inner()[0]),
      expected_ring.coefficients.map(|f| f.inner()[0])
    );
  }

  #[rstest]
  #[case::basic_addition(
      [1, 2, 3, 4, 5, 6, 7, 8],
      [2, 3, 4, 5, 6, 7, 8, 9],
      [3, 5, 7, 9, 11, 13, 15, 0]  // Last value wraps: (8+9)%17=0
  )]
  #[case::with_zero(
      [0, 0, 0, 0, 0, 0, 0, 0],
      [1, 2, 3, 4, 5, 6, 7, 8],
      [1, 2, 3, 4, 5, 6, 7, 8]
  )]
  #[case::wrapping_around_modulus(
      [15, 15, 15, 15, 15, 15, 15, 15],
      [10, 10, 10, 10, 10, 10, 10, 10],
      [8, 8, 8, 8, 8, 8, 8, 8]  // (15+10)%17=8
  )]
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
    let ring_a: Ring<MockField, 8, 8, StandardBasis> = create_ring(a);
    let ring_b: Ring<MockField, 8, 8, StandardBasis> = create_ring(b);
    let expected_ring: Ring<MockField, 8, 8, StandardBasis> = create_ring(expected);

    assert_eq!(
      (ring_a + ring_b).coefficients.map(|f| f.inner()[0]),
      expected_ring.coefficients.map(|f| f.inner()[0])
    );
  }

  #[test]
  fn test_ntt_simple() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input: Ring<MockField, 8, 8, StandardBasis> = create_ring([1, 1, 0, 0, 0, 0, 0, 0]);

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
    let expected: Ring<MockField, 8, 8, NTTBasis> = create_ring([2, 10, 14, 16, 0, 9, 5, 3]);

    assert_eq!(
      ntt_result.coefficients.map(|f| f.inner()[0]),
      expected.coefficients.map(|f| f.inner()[0])
    );
  }

  #[test]
  fn test_inverse_ntt_simple() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input: Ring<MockField, 8, 8, StandardBasis> = create_ring([1, 1, 0, 0, 0, 0, 0, 0]);

    // Compute NTT
    let ntt_result = input.ntt();
    let intt_result = ntt_result.intt();

    assert_eq!(
      input.coefficients.map(|f| f.inner()[0]),
      intt_result.coefficients.map(|f| f.inner()[0])
    );
  }

  #[test]
  fn test_mul() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input_1: Ring<MockField, 8, 8, StandardBasis> = create_ring([1, 1, 0, 0, 0, 0, 0, 0]);
    let input_2: Ring<MockField, 8, 8, StandardBasis> = create_ring([0, 0, 1, 0, 0, 0, 0, 0]);

    // Compute NTT
    let ntt_result_1 = input_1.ntt();
    let ntt_result_2 = input_2.ntt();
    let mul_result = ntt_result_1 * ntt_result_2;
    let intt_result = mul_result.intt();

    // Correct answer
    let expected: Ring<MockField, 8, 8, StandardBasis> = create_ring([0, 0, 1, 1, 0, 0, 0, 0]);

    assert_eq!(
      expected.coefficients.map(|f| f.inner()[0]),
      intt_result.coefficients.map(|f| f.inner()[0])
    );
  }

  #[test]
  fn test_mul_D_16() {
    // Create polynomial 1 + X (coefficients: 1, 1, 0, 0, 0, 0, 0, 0)
    let input_1: Ring<MockField, 16, 8, StandardBasis> =
      create_ring([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    let input_2: Ring<MockField, 16, 8, StandardBasis> =
      create_ring([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    // Compute NTT
    let ntt_result_1 = input_1.ntt();
    let ntt_result_2 = input_2.ntt();
    let mul_result = ntt_result_1 * ntt_result_2;
    let intt_result = mul_result.intt();

    // Correct answer
    let expected: Ring<MockField, 16, 8, StandardBasis> =
      create_ring([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    assert_eq!(
      expected.coefficients.map(|f| f.inner()[0]),
      intt_result.coefficients.map(|f| f.inner()[0])
    );
  }
}
