//! Compile-time utilities for validating and computing field parameters.
//!
//! This module provides const functions for parsing and validating moduli used in
//! cyclotomic ring implementations, particularly for verifying properties needed
//! for the Number Theoretic Transform (NTT).

use super::*;

/// Converts a single ASCII hex character to its numeric value.
///
/// # Arguments
///
/// * `c` - An ASCII character representing a hexadecimal digit
///
/// # Returns
///
/// * `Some(value)` if the character is a valid hex digit (0-9, a-f, A-F)
/// * `None` otherwise
const fn hex_to_digit(c: u8) -> Option<u64> {
  match c {
    b'0'..=b'9' => Some((c - b'0') as u64),
    b'a'..=b'f' => Some((c - b'a' + 10) as u64),
    b'A'..=b'F' => Some((c - b'A' + 10) as u64),
    _ => None,
  }
}

/// Parses a hexadecimal string into a u64 at compile time.
///
/// Supports optional "0x" or "0X" prefix. Returns None if the string is empty,
/// contains invalid characters, or represents a value larger than [`u64::MAX`].
///
/// # Arguments
///
/// * `s` - The string to parse, optionally beginning with "0x" or "0X"
///
/// # Returns
///
/// * `Some(value)` if the string is a valid hexadecimal number
/// * `None` if the string is invalid or empty
const fn parse_hex(s: &str) -> Option<u64> {
  let bytes = s.as_bytes();
  if bytes.is_empty() {
    return None;
  }

  let mut i = if bytes.len() >= 2 && bytes[0] == b'0' && (bytes[1] == b'x' || bytes[1] == b'X') {
    if bytes.len() == 2 {
      return None;
    }
    2 // Skip "0x" prefix
  } else {
    0
  };

  let mut value = 0u64;
  while i < bytes.len() {
    if let Some(digit) = hex_to_digit(bytes[i]) {
      value = value * 16 + digit;
    } else {
      return None;
    }
    i += 1;
  }

  Some(value)
}

pub const fn modulus_to_le_bytes<F: PrimeField>() -> [u8; ((F::NUM_BITS as usize + 7) / 8) * 8] {
  let bytes = F::MODULUS.as_bytes();
  let mut out = [0u8; ((F::NUM_BITS as usize + 7) / 8) * 8];

  let mut i = if bytes.len() >= 2 && bytes[0] == b'0' && (bytes[1] == b'x' || bytes[1] == b'X') {
    if bytes.len() == 2 {
      panic!();
    }
    2 // Skip "0x" prefix
  } else {
    0
  };

  let mut out_idx = 0;
  while i < bytes.len() {
    let hi = hex_to_digit(bytes[i]).unwrap();
    let lo = if i + 1 < bytes.len() { hex_to_digit(bytes[i + 1]).unwrap() } else { 0 };
    out[out_idx] = ((hi << 4) | lo) as u8;
    out_idx += 1;
    i += 2;
  }
  out
}

/// Computes value mod 4t at compile time.
///
/// # Arguments
///
/// * `value` - The value to reduce
/// * `t` - The parameter t where we compute modulo 4t
///
/// # Returns
///
/// The remainder when value is divided by 4t
const fn get_remainder_4t(value: u64, t: usize) -> u64 {
  let four_t = 4 * (t as u64);
  value % four_t
}

/// Computes (value - 1)/d at compile time.
///
/// Used to compute powers for roots of unity in the NTT implementation.
///
/// # Arguments
///
/// * `value` - The modulus value (typically prime)
/// * `d` - The degree parameter of the cyclotomic ring
///
/// # Returns
///
/// The quotient (value - 1)/d
const fn get_unity_power(value: u64, d: usize) -> u64 { (value - 1) / (d as u64) }

/// Verifies at compile time that a prime field's modulus satisfies
/// the congruence condition needed for the NTT implementation.
///
/// Specifically, checks if q ≡ 1 + 2T (mod 4T) where q is the field's
/// modulus and T is the supplied parameter.
///
/// # Type Parameters
///
/// * `F` - The prime field whose modulus we're checking
/// * `T` - The parameter T that determines the congruence condition
///
/// # Returns
///
/// * `true` if the modulus satisfies the congruence condition
/// * `false` otherwise
pub const fn verify_modulus<F: PrimeField, const T: usize>() -> bool {
  if let Some(value) = parse_hex(F::MODULUS) {
    let remainder = get_remainder_4t(value, T);
    let target = 1 + 2 * (T as u64);
    remainder == target % (4 * (T as u64))
  } else {
    false
  }
}

/// Computes at compile time the power needed for the NTT's root of unity.
///
/// For a field of order q and degree parameter D, computes (q-1)/D.
/// This value is used to determine the appropriate root of unity for
/// the NTT implementation.
///
/// # Type Parameters
///
/// * `F` - The prime field we're working in
/// * `D` - The degree parameter of the cyclotomic ring
///
/// # Returns
///
/// The power (q-1)/D where q is the field's modulus
pub const fn unity_power<F: PrimeField, const D: usize>() -> u64 {
  match parse_hex(F::MODULUS) {
    Some(value) => get_unity_power(value, D),
    None => 0,
  }
}

#[cfg(test)]
mod tests {
  use ff::Field;
  use rstest::rstest;

  use super::*;

  #[rstest]
  #[case(b'0', 0)]
  #[case(b'1', 1)]
  #[case(b'9', 9)]
  #[case(b'a', 10)]
  #[case(b'f', 15)]
  #[case(b'A', 10)]
  #[case(b'F', 15)]
  fn test_hex_to_digit_valid(#[case] input: u8, #[case] expected: u64) {
    assert_eq!(hex_to_digit(input).unwrap(), expected);
  }

  #[rstest]
  #[case(b'g')]
  #[case(b'G')]
  #[case(b'x')]
  #[case(b' ')]
  #[case(b'$')]
  fn test_hex_to_digit_invalid(#[case] input: u8) {
    assert!(hex_to_digit(input).is_none());
  }

  #[rstest]
  #[case("0x11", 17)]
  #[case("0xff", 255)]
  #[case("17", 23)]
  #[case("0x0", 0)]
  fn test_parse_hex_valid(#[case] input: &str, #[case] expected: u64) {
    assert_eq!(parse_hex(input).unwrap(), expected);
  }

  #[rstest]
  #[case("")] // Empty string
  #[case("0x")] // Only prefix
  #[case("0xG")] // Invalid hex digit
  #[case("0x 1")] // Space in number
  #[case("$17")] // Invalid prefix character
  fn test_parse_hex_invalid(#[case] input: &str) {
    assert!(parse_hex(input).is_none());
  }

  #[rstest]
  #[case(17, 8, 2)] // (17-1)/8 = 2
  #[case(41, 8, 5)] // (41-1)/8 = 5
  #[case(16, 4, 3)] // (16-1)/4 = 3
  fn test_get_unity_power(#[case] value: u64, #[case] d: usize, #[case] expected: u64) {
    assert_eq!(get_unity_power(value, d), expected);
  }

  #[rstest]
  #[case(17, 4, 1)] // 17 ≡ 1 (mod 4*4)
  #[case(41, 10, 1)] // 41 ≡ 1 (mod 4*10)
  fn test_get_remainder_4t(#[case] value: u64, #[case] t: usize, #[case] expected: u64) {
    assert_eq!(get_remainder_4t(value, t), expected);
  }

  #[test]
  fn test_verify_modulus() {
    assert!(verify_modulus::<MockField, 8>());
  }

  #[test]
  fn test_unity_power_with_field() {
    assert_eq!(unity_power::<MockField, 8>(), 2);
  }

  #[test]
  #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test)]
  fn test_modulus_to_le_bytes() {
    // For modulus "17" (0x11):
    // In binary: 0001 0001
    // In little endian: 11 00 00 00 ... (remaining bytes are 0)
    let bytes = modulus_to_le_bytes::<MockField>();
    assert_eq!(bytes, [0x11, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);
  }
}
