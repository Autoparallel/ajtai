//! Homomorphic commitment scheme implementation based on Ring-LWE.
//!
//! This module provides an implementation of a lattice-based commitment scheme that is
//! homomorphic with respect to addition. The scheme can be instantiated in two modes:
//! - A non-hiding variant that provides only binding security
//! - A hiding variant that provides both hiding and binding security
//!
//! The scheme operates over cyclotomic rings and uses a random matrix for commitments.
//! Security is based on the Ring-LWE (Ring Learning With Errors) problem.
//!
//! # Examples
//!
//! ```ignore
//! use ajtai::{
//!   commitment::{CommitmentKey, NonHiding},
//!   ring::{CyclotomicRing, StandardBasis},
//! };
//!
//! // Create a commitment key with:
//! // - Field: MockField (modulus 17)
//! // - Ring degree: D = 8
//! // - Parameter T = 8
//! // - Matrix dimensions: 2 x 3
//! let key = CommitmentKey::<MockField, 8, 8, 2, 3, NonHiding>::setup();
//!
//! // Create a message to commit to
//! let message = [
//!   CyclotomicRing::new([1, 0, 0, 0, 0, 0, 0, 0]),
//!   CyclotomicRing::new([0, 1, 0, 0, 0, 0, 0, 0]),
//!   CyclotomicRing::new([0, 0, 1, 0, 0, 0, 0, 0]),
//! ];
//!
//! // Commit to the message
//! let commitment = key.commit(message);
//! ```
//!
//! # Security Parameters
//!
//! The security of the scheme depends on several parameters:
//! - `D`: The degree of the cyclotomic ring (must be a power of 2)
//! - `T`: A parameter dividing D that determines properties of the modulus
//! - `K`: Number of rows in the commitment matrix
//! - `M`: Number of columns in the commitment matrix
//!
//! # Type Parameters
//!
//! * `F`: The prime field used for ring coefficients
//! * `D`: Degree of the cyclotomic ring (must be a power of 2)
//! * `T`: Parameter dividing D that determines modulus properties
//! * `K`: Number of rows in the commitment matrix
//! * `M`: Number of columns in the commitment matrix
//! * `H`: Type indicating whether the scheme is hiding or non-hiding

use comptime::verify_modulus;
use rand_core::OsRng;
use ring::{CyclotomicRing, StandardBasis};

use super::*;

mod sealed {
  pub trait Sealed {}
}

/// Marker trait for distinguishing between hiding and non-hiding commitment variants.
///
/// This trait is sealed and cannot be implemented outside this crate.
pub trait IsHiding: sealed::Sealed {
  /// The type to use for the additional matrix when doing hiding
  type HidingMatrix;
  /// The type to use for the additional vector when doing hiding
  type HidingOpening;
}

/// Type indicating a non-hiding commitment scheme variant.
///
/// In this variant, commitments provide binding security but not hiding security.
pub struct NonHiding;

/// Type indicating a hiding commitment scheme variant with parameters for the hiding mechanism.
///
/// * `R`: The ring element type used for the hiding matrix
/// * `M`: Number of columns in the commitment matrix
/// * `V`: Additional parameter for the hiding mechanism
pub struct Hiding<R, const M: usize, const V: usize>(PhantomData<R>);

impl<R, const M: usize, const V: usize> sealed::Sealed for Hiding<R, M, V> {}
impl<R, const M: usize, const V: usize> IsHiding for Hiding<R, M, V>
where [(); M + V]:
{
  type HidingMatrix = [[R; M]; M + V];
  type HidingOpening = [R; M + V];
}

impl sealed::Sealed for NonHiding {}
impl IsHiding for NonHiding {
  type HidingMatrix = ();
  type HidingOpening = ();
}

/// A commitment key containing the public parameters for the commitment scheme.
///
/// The key consists of:
/// - A random matrix of ring elements
/// - Additional parameters for hiding in the hiding variant
pub struct CommitmentKey<
  F: PrimeField,
  const D: usize,
  const T: usize,
  const K: usize,
  const M: usize,
  H: IsHiding,
> {
  /// Random matrix used whether hiding or not
  matrix:        [[CyclotomicRing<F, D, T, StandardBasis>; M]; K],
  /// Random matrix also used when hiding
  #[allow(dead_code)]
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
  /// Generates a new commitment key with random parameters.
  ///
  /// Creates a random K×M matrix where each entry is a ring element
  /// with random coefficients. This matrix will be used to compute
  /// commitments.
  ///
  /// # Returns
  ///
  /// A new `CommitmentKey` instance with randomly generated parameters.
  pub fn setup() -> Self {
    let rng = OsRng;
    let mut matrix = [[CyclotomicRing::<F, D, T, StandardBasis>::DEFAULT; M]; K];

    #[allow(clippy::needless_range_loop)]
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
  /// Commits to a message using the commitment key.
  ///
  /// Computes the commitment as a matrix-vector product between the commitment
  /// key matrix and the message vector. The result is a binding commitment
  /// to the message.
  ///
  /// # Arguments
  ///
  /// * `message` - Array of M ring elements representing the message to commit to
  ///
  /// # Returns
  ///
  /// A `Commitment` containing the commitment value and opening information.
  pub fn commit(
    &self,
    message: [CyclotomicRing<F, D, T, StandardBasis>; M],
  ) -> Commitment<F, D, T, K, M, NonHiding> {
    let mut commitment = [CyclotomicRing::<F, D, T, StandardBasis>::DEFAULT; K];
    #[allow(clippy::needless_range_loop)]
    for i in 0..K {
      #[allow(clippy::needless_range_loop)]
      for j in 0..M {
        commitment[i] += self.matrix[i][j] * message[j];
      }
    }
    Commitment { commitment, opening_message: message, opening_randomness: () }
  }
}

/// A commitment along with data needed for opening.
///
/// Contains:
/// - The commitment value itself (a vector of ring elements)
/// - The original message
/// - Randomness used for hiding in the hiding variant
pub struct Commitment<
  F: PrimeField,
  const D: usize,
  const T: usize,
  const K: usize,
  const M: usize,
  H: IsHiding,
> {
  /// The underlying commitment
  commitment:         [CyclotomicRing<F, D, T, StandardBasis>; K],
  /// The message used to generate the commitment
  #[allow(dead_code)]
  opening_message:    [CyclotomicRing<F, D, T, StandardBasis>; M],
  /// The vector used to add randomness when hiding
  #[allow(dead_code)]
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
  /// Computes the supremum norm of the commitment.
  ///
  /// The supremum norm is the maximum of the sup norms of all ring elements
  /// in the commitment. For each ring element, its sup norm is the maximum
  /// of the minimal absolute values of its coefficients modulo q.
  ///
  /// # Returns
  ///
  /// The supremum norm of the commitment as a u64.
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
