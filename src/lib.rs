#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![no_std]
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

#[cfg(test)]
#[macro_use]
extern crate std;

#[cfg(test)] use rstest::rstest;

pub mod ring;
