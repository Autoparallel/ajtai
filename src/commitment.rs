use comptime::verify_modulus;
use rand_core::{OsRng, RngCore};
use ring::{CyclotomicRing, StandardBasis};

use super::*;

pub struct Commitment<
  F: PrimeField,
  const D: usize,
  const T: usize,
  const K: usize,
  const M: usize,
  const B: usize,
> {
  matrix: [[CyclotomicRing<F, D, T, StandardBasis>; K]; M],
}

impl<
    F: PrimeField,
    const D: usize,
    const T: usize,
    const K: usize,
    const M: usize,
    const B: usize,
  > Commitment<F, D, T, K, M, B>
where
  // D is a power of two
  [(); D.is_power_of_two() as usize - 1]:,
  // T is a divisor of D
  [(); (D % T == 0) as usize - 1]:,
  // q = 1 + 2*T (mod 4*T)
  [(); verify_modulus::<F, T>() as usize - 1]:,
{
  pub fn setup() -> Self {
    let rng = OsRng;
    let mut matrix = [[CyclotomicRing::<F, D, T, StandardBasis>::DEFAULT; K]; M];

    for i in 0..M {
      for j in 0..K {
        // Generate random coefficients for each ring element
        let coefficients = core::array::from_fn(|_| F::random(rng));

        matrix[i][j] = CyclotomicRing::from(coefficients);
      }
    }

    Self { matrix }
  }
}

#[cfg(test)]
mod tests {
  use ff::Field;

  use super::*;

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_setup() {
    // Create a commitment instance
    let commitment = Commitment::<MockField, 16, 8, 4, 3, 1>::setup();

    // Check dimensions
    assert_eq!(commitment.matrix.len(), 3); // M rows
    assert_eq!(commitment.matrix[0].len(), 4); // K columns

    // Verify elements are not all zero (with high probability)
    let mut all_zero = true;
    'outer: for i in 0..3 {
      for j in 0..4 {
        if commitment.matrix[i][j].coefficients.iter().any(|&x| x != MockField::ZERO) {
          all_zero = false;
          break 'outer;
        }
      }
    }
    assert!(
      !all_zero,
      "All elements were zero, which is extremely unlikely with random generation"
    );

    // Test that multiple calls generate different matrices
    let commitment2 = Commitment::<MockField, 16, 8, 4, 3, 1>::setup();
    assert!(
      commitment.matrix != commitment2.matrix,
      "Two random setups produced identical matrices"
    );
  }
}
