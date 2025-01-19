use core::ops::{Add, Mul};

use comptime::verify_modulus;

use super::*;

// First define our basis tags as zero-sized types
pub struct StandardBasis;
pub struct NTTBasis;

// Modify Ring to take a basis parameter
#[derive(Debug, Clone, Copy)]
pub struct Ring<F: PrimeField, const D: usize, const T: usize, Basis = StandardBasis>
where
  // D is a power of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:,
  // q = 1 + 2*T (mod 4*T)
  [(); verify_modulus::<F, T>() as usize - 1]:, {
  coefficients: [F; D],
  _basis:       core::marker::PhantomData<Basis>,
}

// Implementation for standard basis
impl<F: PrimeField, const D: usize, const T: usize> Ring<F, D, T, StandardBasis>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  pub fn ntt(self) -> Ring<F, D, T, NTTBasis> {
    let omega = F::MULTIPLICATIVE_GENERATOR;
    let mut result = self.coefficients;

    // Length is D which is power of 2
    let n = D;
    let mut m = 1; // Subarray size starts at 1

    // Cooley-Tukey NTT
    while m < n {
      let w_m = omega.pow([(n / (2 * m)) as u64]);

      for k in (0..n).step_by(2 * m) {
        let mut w = F::ONE;

        for j in 0..m {
          let t = w * result[k + j + m];
          result[k + j + m] = result[k + j] - t;
          result[k + j] += t;
          w *= w_m;
        }
      }

      m *= 2;
    }

    Ring { coefficients: result, _basis: core::marker::PhantomData }
  }
}

impl<F: PrimeField, const D: usize, const T: usize, B> Add for Ring<F, D, T, B>
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

#[cfg(test)]
mod tests {
  use core::marker::PhantomData;

  use super::*;

  #[test]
  fn test_valid_dimensions() {
    let _: Ring<MockField, 8, 8> = create_ring([1; 8]);
    let _: Ring<MockField, 16, 8> = create_ring([1; 16]);
  }

  // Helper function to create Ring instances more easily
  fn create_ring<const D: usize, const T: usize>(values: [u64; D]) -> Ring<MockField, D, T>
  where
    [(); D.is_power_of_two() as usize - 1]:,
    [(); (D % T == 0) as usize - 1]:,
    [(); verify_modulus::<MockField, T>() as usize - 1]:, {
    Ring { coefficients: values.map(MockField::from), _basis: PhantomData::<StandardBasis> }
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
    let ring_a: Ring<MockField, 16, 8> = create_ring(a);
    let ring_b: Ring<MockField, 16, 8> = create_ring(b);
    let expected_ring: Ring<MockField, 16, 8> = create_ring(expected);

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
    let ring_a: Ring<MockField, 8, 8> = create_ring(a);
    let ring_b: Ring<MockField, 8, 8> = create_ring(b);
    let expected_ring: Ring<MockField, 8, 8> = create_ring(expected);

    assert_eq!(
      (ring_a + ring_b).coefficients.map(|f| f.inner()[0]),
      expected_ring.coefficients.map(|f| f.inner()[0])
    );
  }
}
