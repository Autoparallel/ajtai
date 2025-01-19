#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![no_std]
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

#[cfg(test)]
#[macro_use]
extern crate std;

use ff::PrimeField;
#[cfg(test)] use {mock::MockField, rstest::rstest};

pub mod comptime;
pub mod ring;

#[cfg(test)]
mod mock {
  use super::*;

  #[derive(PrimeField)]
  #[PrimeFieldModulus = "17"]
  #[PrimeFieldGenerator = "2"]
  #[PrimeFieldReprEndianness = "little"]
  pub struct MockField([u64; 1]);

  impl MockField {
    pub const fn inner(&self) -> &[u64] { &self.0 }
  }
}
