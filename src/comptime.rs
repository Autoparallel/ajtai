use super::*;

const fn hex_to_digit(c: u8) -> Option<u64> {
  match c {
    b'0'..=b'9' => Some((c - b'0') as u64),
    b'a'..=b'f' => Some((c - b'a' + 10) as u64),
    b'A'..=b'F' => Some((c - b'A' + 10) as u64),
    _ => None,
  }
}

const fn get_remainder_4T<const T: usize>(s: &str) -> Option<u64> {
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

  let mut remainder = 0u64;
  let four_t = 4 * (T as u64);

  while i < bytes.len() {
    if let Some(digit) = hex_to_digit(bytes[i]) {
      remainder = ((remainder * 16) % four_t + digit % four_t) % four_t;
    } else {
      return None;
    }
    i += 1;
  }

  Some(remainder)
}

pub const fn verify_modulus<F: PrimeField, const T: usize>() -> bool {
  if let Some(remainder) = get_remainder_4T::<T>(F::MODULUS) {
    let target = 1 + 2 * (T as u64);
    remainder == target % (4 * (T as u64))
  } else {
    false
  }
}

#[cfg(test)]
mod tests {
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
  #[case("0x11", 17)] // 0x11 = 17 (mod 32)
  #[case("0xFF", 31)] // 0xFF = 255 ≡ 31 (mod 32)
  #[case("17", 23)] // 17 (mod 32) = 17
  #[case("0x17", 23)] // 0x17 = 23 (mod 32)
  #[case("0x10", 16)] // 0x10 = 16 (mod 32)
  fn test_get_remainder_valid(#[case] input: &str, #[case] expected: u64) {
    assert_eq!(get_remainder_4T::<8>(input).unwrap(), expected);
  }

  #[rstest]
  #[case("")] // Empty after prefix
  #[case("0xG")] // Invalid hex digit
  #[case("0x 1")] // Space in number
  #[case("$17")] // Invalid prefix character
  fn test_get_remainder_invalid(#[case] input: &str) {
    assert!(get_remainder_4T::<8>(input).is_none());
  }

  #[test]
  fn test_verify_modulus() {
    // Test with T = 1 (smallest valid T)
    assert!(!verify_modulus::<MockField, 1>());

    // Test with larger T values
    assert!(!verify_modulus::<MockField, 16>());
  }
}
