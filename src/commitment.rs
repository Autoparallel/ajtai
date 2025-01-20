use comptime::verify_modulus;
use rand_core::OsRng;
use ring::{CyclotomicRing, StandardBasis};

use super::*;

mod sealed {
  pub trait Sealed {}
}

pub trait IsHiding: sealed::Sealed {
  type HidingMatrix;
  type HidingOpening;
}

pub struct NonHiding;

pub struct Hiding<R, const M: usize, const V: usize>(PhantomData<R>);

impl<R, const M: usize, const V: usize> sealed::Sealed for Hiding<R, M, V> {}
impl sealed::Sealed for NonHiding {}

impl<R, const M: usize, const V: usize> IsHiding for Hiding<R, M, V>
where [(); M + V]:
{
  type HidingMatrix = [[R; M]; M + V];
  type HidingOpening = [R; M + V];
}
impl IsHiding for NonHiding {
  type HidingMatrix = ();
  type HidingOpening = ();
}

pub struct CommitmentKey<
  F: PrimeField,
  const D: usize,
  const T: usize,
  const K: usize,
  const M: usize,
  H: IsHiding,
> {
  matrix:        [[CyclotomicRing<F, D, T, StandardBasis>; M]; K],
  hiding_matrix: H::HidingMatrix,
}

impl<F: PrimeField, const D: usize, const T: usize, const K: usize, const M: usize>
  CommitmentKey<F, D, T, K, M, NonHiding>
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
    let mut matrix = [[CyclotomicRing::<F, D, T, StandardBasis>::DEFAULT; M]; K];

    for i in 0..K {
      for j in 0..M {
        // Generate random coefficients for each ring element
        let coefficients = core::array::from_fn(|_| F::random(rng));

        matrix[i][j] = CyclotomicRing::from(coefficients);
      }
    }

    Self { matrix, hiding_matrix: () }
  }

  // TODO: We can implement a more optimal mat mul?
  pub fn commit(
    &self,
    message: [CyclotomicRing<F, D, T, StandardBasis>; M],
  ) -> Commitment<F, D, T, K, M, NonHiding> {
    let mut commitment = [CyclotomicRing::<F, D, T, StandardBasis>::DEFAULT; K];
    for i in 0..K {
      for j in 0..M {
        commitment[i] += self.matrix[i][j] * message[j];
      }
    }
    Commitment { commitment, opening_message: message, opening_randomness: () }
  }
}

pub struct Commitment<
  F: PrimeField,
  const D: usize,
  const T: usize,
  const K: usize,
  const M: usize,
  H: IsHiding,
> {
  commitment:         [CyclotomicRing<F, D, T, StandardBasis>; K],
  opening_message:    [CyclotomicRing<F, D, T, StandardBasis>; M],
  opening_randomness: H::HidingOpening,
}

impl<F: PrimeField, const D: usize, const T: usize, const K: usize, const M: usize>
  Commitment<F, D, T, K, M, NonHiding>
where
  [(); D.is_power_of_two() as usize - 1]:,
  [(); (D % T == 0) as usize - 1]:,
  [(); verify_modulus::<F, T>() as usize - 1]:,
  [(); ((F::NUM_BITS as usize + 7) / 8) * 8]:,
{
  pub fn sup_norm(&self) -> u64 {
    self.commitment.iter().map(CyclotomicRing::sup_norm).max().unwrap_or(0)
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
    let commitment = CommitmentKey::<MockField, 16, 8, 10, 20, NonHiding>::setup();

    // Check dimensions
    assert_eq!(commitment.matrix[0].len(), 20); // K rows
    assert_eq!(commitment.matrix.len(), 10); // M columns

    // Verify elements are not all zero (with high probability)
    let mut all_zero = true;
    'outer: for i in 0..10 {
      for j in 0..20 {
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
    let commitment2 = CommitmentKey::<MockField, 16, 8, 10, 20, NonHiding>::setup();
    assert!(
      commitment.matrix != commitment2.matrix,
      "Two random setups produced identical matrices"
    );
  }

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_commitment_sup_norm() {
    // Create a commitment with known values
    let commitment = [
      CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([1, 2, 3, 4, 5, 6, 7, 8]),
      CyclotomicRing::<MockField, 8, 8, StandardBasis>::new([0, 16, 15, 14, 13, 12, 11, 10]),
    ];
    let commitment = Commitment::<MockField, 8, 8, 2, 1, NonHiding> {
      opening_message: [CyclotomicRing::DEFAULT; 1],
      commitment,
      opening_randomness: (),
    };

    // The first ring has max value 8
    // The second ring's values in minimal representation are:
    // 0, 1 (17-16), 2 (17-15), 3 (17-14), 4 (17-13), 5 (17-12), 6 (17-11), 7 (17-10)
    // So the maximum across both rings should be 8
    assert_eq!(commitment.sup_norm(), 8);
  }
}
