use super::*;

const fn hex_to_digit(c: u8) -> Option<u64> {
  match c {
    b'0'..=b'9' => Some((c - b'0') as u64),
    b'a'..=b'f' => Some((c - b'a' + 10) as u64),
    b'A'..=b'F' => Some((c - b'A' + 10) as u64),
    _ => None,
  }
}

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

const fn get_remainder_4t(value: u64, t: usize) -> u64 {
  let four_t = 4 * (t as u64);
  value % four_t
}

const fn get_unity_power(value: u64, d: usize) -> u64 { (value - 1) / (d as u64) }

pub const fn verify_modulus<F: PrimeField, const T: usize>() -> bool {
  if let Some(value) = parse_hex(F::MODULUS) {
    let remainder = get_remainder_4t(value, T);
    let target = 1 + 2 * (T as u64);
    remainder == target % (4 * (T as u64))
  } else {
    false
  }
}

pub const fn unity_power<F: PrimeField, const D: usize>() -> u64 {
  match parse_hex(F::MODULUS) {
    Some(value) => get_unity_power(value, D),
    None => 0,
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
}
