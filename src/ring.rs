use core::ops::{Add, Mul};

use super::*;

#[derive(Debug, Clone, Copy)]
pub struct Ring<F: PrimeField, const D: usize, const T: usize>
where
  // D is a power of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:, {
  inner: [F; D],
}

impl<F: PrimeField, const D: usize, const T: usize> Add for Ring<F, D, T>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
{
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    Self { inner: core::array::from_fn(|i| self.inner[i] + rhs.inner[i]) }
  }
}

// impl<F: PrimeField, const D: usize> Mul for Ring<F, D>
// where [(); D.is_power_of_two() as usize - 1]:
// {
//   type Output = Self;

//   fn mul(self, rhs: Self) -> Self::Output { todo!() }
// }

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_valid_dimensions() {
    let _: Ring<MockField, 1, 1> = create_ring([1]); // Works
    let _: Ring<MockField, 2, 2> = create_ring([1, 2]); // Works
    let _: Ring<MockField, 4, 2> = create_ring([1, 2, 3, 4]); // Works
  }

  // Helper function to create Ring instances more easily
  fn create_ring<const D: usize, const T: usize>(values: [u64; D]) -> Ring<MockField, D, T>
  where
    [(); D.is_power_of_two() as usize - 1]:,
    [(); (D % T == 0) as usize - 1]:, {
    Ring { inner: values.map(MockField::from) }
  }

  #[rstest]
  #[case::basic_addition(
    [5],
    [7],
    [12]
)]
  #[case::with_zero(
    [0],
    [8],
    [8]
)]
  #[case::wrapping_around_modulus(
    [15],
    [10],
    [6]  // (15+10)%19=6
)]
  #[case::large_numbers(
    [18],
    [18],
    [17]  // (18+18)%19=17
)]
  fn test_add_dimension_one(#[case] a: [u64; 1], #[case] b: [u64; 1], #[case] expected: [u64; 1]) {
    let ring_a: Ring<MockField, 1, 1> = create_ring(a);
    let ring_b: Ring<MockField, 1, 1> = create_ring(b);
    let expected_ring: Ring<MockField, 1, 1> = create_ring(expected);

    assert_eq!(
      (ring_a + ring_b).inner.map(|f| f.inner()[0]),
      expected_ring.inner.map(|f| f.inner()[0])
    );
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
  fn test_add_dimension_two(#[case] a: [u64; 2], #[case] b: [u64; 2], #[case] expected: [u64; 2]) {
    let ring_a: Ring<MockField, 2, 1> = create_ring(a);
    let ring_b: Ring<MockField, 2, 1> = create_ring(b);
    let expected_ring: Ring<MockField, 2, 1> = create_ring(expected);

    assert_eq!(
      (ring_a + ring_b).inner.map(|f| f.inner()[0]),
      expected_ring.inner.map(|f| f.inner()[0])
    );
  }

  #[rstest]
  #[case::basic_addition(
      [1, 2, 3, 4],
      [5, 6, 7, 8],
      [6, 8, 10, 12]
  )]
  #[case::with_zero(
      [0, 2, 0, 4],
      [5, 0, 7, 0],
      [5, 2, 7, 4]
  )]
  #[case::wrapping_around_modulus(
      [15, 16, 17, 18],
      [10, 11, 12, 13],
      [6, 8, 10, 12]  // All values mod 19
  )]
  #[case::large_numbers(
      [18, 18, 18, 18],
      [18, 18, 18, 18],
      [17, 17, 17, 17]  // (18+18)%19=17
  )]
  fn test_add_dimension_four(#[case] a: [u64; 4], #[case] b: [u64; 4], #[case] expected: [u64; 4]) {
    let ring_a: Ring<MockField, 4, 2> = create_ring(a);
    let ring_b: Ring<MockField, 4, 2> = create_ring(b);
    let expected_ring: Ring<MockField, 4, 2> = create_ring(expected);

    assert_eq!(
      (ring_a + ring_b).inner.map(|f| f.inner()[0]),
      expected_ring.inner.map(|f| f.inner()[0])
    );
  }
}
