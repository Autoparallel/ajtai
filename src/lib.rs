#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![no_std]
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use core::marker::PhantomData;

use ff::PrimeField;
#[cfg(all(target_arch = "wasm32", test))]
use wasm_bindgen_test::wasm_bindgen_test;
#[cfg(test)] use {mock::MockField, rstest::rstest};

#[cfg(test)]
#[macro_use]
extern crate std;

pub mod commitment;
pub mod comptime;
pub mod ring;

#[cfg(test)]
mod mock {
  use super::*;

  #[derive(PrimeField)]
  #[PrimeFieldModulus = "17"]
  // NOTE: This does NOT verify the number is actually a generator! FML
  #[PrimeFieldGenerator = "3"]
  #[PrimeFieldReprEndianness = "little"]
  pub struct MockField([u64; 1]);
}
