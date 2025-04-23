//! Unchecked free functions, these are called by implementations and the [`checked`](crate::checked) functions

#[inline(always)]
pub const fn nside(depth: u8) -> u32 {
    1u32 << depth
}
#[inline(always)]
pub const fn nside_square(depth: u8) -> u32 {
    1u32 << (depth << 1)
}

/// Same as [checked::n_hash](crate::checked::n_hash) except that this version does not panic if the given `depth` is
/// out of range.
#[inline]
pub const fn n_hash(depth: u8) -> u64 {
    12u64 << (depth << 1u8)
}

/// Transforms the given NESTED hash value into its uniq representation, i.e. the depth
/// is encoded together with the hash value such that each possible (deph, hash) pair
/// gives a unique number.
/// In practice, the unique representation uses a sentinel bit (set to one) to code the depth.
/// The sentinel bit (set to 1) is the least significant bit (LSB) among the unused bits,
/// i.e. the most significant bit (MSB) located just after the hash MSB.
/// Said differently, the sentinel bit is the `(1 + 4 + 2*depth)^th` MSB
/// The encoding, in the case of the nested scheme, is thus `0...0sbbbb112233...`, with:
/// * `0...0`: unused bits
/// * `s` : sentinel bit
/// * `bbbb`: the 4 bits coding the base cell
/// * `11`: the 2 bits coding depth 1
/// * `22`: the 2 bits coding depth 2
/// * `33`: the 2 bits coding depth 3
/// * ...
///
/// # Example
///
/// ```rust
/// use cdshealpix::nested::{get, Layer};
/// let l0 = get(0);
/// assert_eq!(l0.to_uniq(0), 16);
/// ```
#[inline]
pub const fn to_uniq(depth: u8, hash: u64) -> u64 {
    (16_u64 << (depth << 1)) | hash
}

/// Same as [to_uniq](fn.to_uniq.html), but
/// following the [IVOA](http://ivoa.net/documents/MOC/) convention.
/// It does not rely on a sentinel bit and use one less bit.
#[inline]
pub const fn to_uniq_ivoa(depth: u8, hash: u64) -> u64 {
    (4_u64 << (depth << 1)) + hash
}
