#[cfg(test)]
use super::*;

use core::ops::{Add, Mul};
use ff::PrimeField;

#[derive(Debug, Clone, Copy)]
pub struct Ring<F: PrimeField, const D: usize>([F; D]);

impl<F: PrimeField, const D: usize> Add for Ring<F, D> {
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    Self(core::array::from_fn(|i| self.0[i] + rhs.0[i]))
  }
}

impl<F: PrimeField, const D: usize> Mul for Ring<F, D> {
  type Output = Self;

  fn mul(self, rhs: Self) -> Self::Output {
    todo!()
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[derive(PrimeField)]
  #[PrimeFieldModulus = "19"]
  #[PrimeFieldGenerator = "2"]
  #[PrimeFieldReprEndianness = "little"]
  struct MockField([u64; 1]);

  // Helper function to create Ring instances more easily
  fn create_ring<const D: usize>(values: [u64; D]) -> Ring<MockField, D> {
    Ring(values.map(MockField::from))
  }

  #[rstest]
  #[case::basic_addition(
        [1, 2],
        [3, 4],
        [4, 6]
    )]
  #[case::with_zero(
        [0, 5],
        [7, 0],
        [7, 5]
    )]
  #[case::wrapping_around_modulus(
        [10, 15],
        [12, 8],
        [3, 4]  // (10+12)%19=3, (15+8)%19=4
    )]
  #[case::large_numbers(
        [18, 18],
        [18, 18],
        [17, 17]  // (18+18)%19=17
    )]
  fn test_add<const D: usize>(
    #[case] a: [u64; D],
    #[case] b: [u64; D],
    #[case] expected: [u64; D],
  ) {
    let ring_a = create_ring(a);
    let ring_b = create_ring(b);
    let expected_ring = create_ring(expected);

    assert_eq!((ring_a + ring_b).0.map(|f| f.0[0]), expected_ring.0.map(|f| f.0[0]));
  }

  // Test specific cases for dimension 1
  #[rstest]
  fn test_add_dimension_one() {
    let ring_a: Ring<MockField, 1> = create_ring([5]);
    let ring_b: Ring<MockField, 1> = create_ring([7]);
    let expected: Ring<MockField, 1> = create_ring([12]);

    assert_eq!((ring_a + ring_b).0.map(|f| f.0[0]), expected.0.map(|f| f.0[0]));
  }

  // Test specific cases for dimension 3
  #[rstest]
  fn test_add_dimension_three() {
    let ring_a: Ring<MockField, 3> = create_ring([1, 2, 3]);
    let ring_b: Ring<MockField, 3> = create_ring([4, 5, 6]);
    let expected: Ring<MockField, 3> = create_ring([5, 7, 9]);

    assert_eq!((ring_a + ring_b).0.map(|f| f.0[0]), expected.0.map(|f| f.0[0]));
  }
}
