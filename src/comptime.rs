use super::*;

const fn hex_to_digit(c: u8) -> Option<u64> {
  match c {
    b'0'..=b'9' => Some((c - b'0') as u64),
    b'a'..=b'f' => Some((c - b'a' + 10) as u64),
    b'A'..=b'F' => Some((c - b'A' + 10) as u64),
    _ => None,
  }
}

// Parse hex string into array of u64 limbs (little-endian)
const fn parse_hex_to_limbs<const N: usize>(s: &str) -> Option<[u64; N]> {
  let mut result = [0u64; N];
  let bytes = s.as_bytes();
  let mut i = if bytes.len() >= 2 && bytes[0] == b'0' && (bytes[1] == b'x' || bytes[1] == b'X') {
    2 // Skip "0x" prefix
  } else {
    0
  };

  let mut value = 0u64;
  while i < bytes.len() {
    if let Some(digit) = hex_to_digit(bytes[i]) {
      value = value * 16 + digit;
    } else {
      return None; // Invalid hex digit
    }
    i += 1;
  }

  result[0] = value;
  Some(result)
}

const fn get_remainder_4T<const T: usize>(s: &str) -> Option<u64> {
  let bytes = s.as_bytes();
  let mut i = if bytes.len() >= 2 && bytes[0] == b'0' && (bytes[1] == b'x' || bytes[1] == b'X') {
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

const fn verify_modulus<F: PrimeField, const T: usize>() -> bool {
  if let Some(remainder) = get_remainder_4T::<T>(F::MODULUS) {
    let target = 1 + 2 * (T as u64);
    remainder == target % (4 * (T as u64))
  } else {
    false
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_hex_to_digit() {
    let digit = hex_to_digit(b"c"[0]).unwrap();
    assert_eq!(digit, 12)
  }

  #[test]
  fn test_get_remainder() {
    let remainder = get_remainder_4T::<8>(MockField::MODULUS).unwrap();
    assert_eq!(remainder, 17);
  }

  #[test]
  fn test_verify_modulus() {
    let verified = verify_modulus::<MockField, 8>();
    assert!(verified);
  }

  //   #[test]
  //   fn test_parse_hex_to_limbs() {
  //     let limbs = parse_hex_to_limbs(MockField::MODULUS).unwrap();
  //     assert_eq!(limbs, [19; 1]);
  //   }
}
