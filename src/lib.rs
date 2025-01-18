#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![no_std]

#[cfg(test)]
#[macro_use]
extern crate std;

#[cfg(test)]
use rstest::rstest;

pub mod ring;
