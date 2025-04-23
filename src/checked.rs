use super::*;

#[inline(always)]
pub const fn nside(depth: u8) -> u32 {
    assert_valid_depth(depth);
    unchecked::nside(depth)
}
#[inline(always)]
pub const fn nside_square(depth: u8) -> u32 {
    assert_valid_depth(depth);
    unchecked::nside_square(depth)
}

/// Returns the number of distinct hash value (the number of cells or pixel the unit sphere is
/// devided in) at the given `depth`.
///
/// # Input
/// - `depth` must be in `[0, 29]`
///
/// # Output
/// - `n_hash` = `12 * nside^2`
///
/// # Panics
/// If `depth` is not valid (see [is_depth](fn.is_depth.html)), this method panics.
///
/// # Examples
///
/// ```rust
/// use cdshealpix::{n_hash};
///
/// assert_eq!(12u64, n_hash(0u8));
/// assert_eq!(48u64, n_hash(1u8));
/// assert_eq!(192u64, n_hash(2u8));
/// assert_eq!(768u64, n_hash(3u8));
/// assert_eq!(3072u64, n_hash(4u8));
/// assert_eq!(12288u64, n_hash(5u8));
/// assert_eq!(49152u64, n_hash(6u8));
/// assert_eq!(196608u64, n_hash(7u8));
/// assert_eq!(786432u64, n_hash(8u8));
/// assert_eq!(3145728u64, n_hash(9u8));
/// assert_eq!(12582912u64, n_hash(10u8));
/// assert_eq!(50331648u64, n_hash(11u8));
/// assert_eq!(201326592u64, n_hash(12u8));
/// assert_eq!(805306368u64, n_hash(13u8));
/// assert_eq!(3221225472u64, n_hash(14u8));
/// assert_eq!(12884901888u64, n_hash(15u8));
/// assert_eq!(51539607552u64, n_hash(16u8));
/// assert_eq!(206158430208u64, n_hash(17u8));
/// assert_eq!(824633720832u64, n_hash(18u8));
/// assert_eq!(3298534883328u64, n_hash(19u8));
/// assert_eq!(13194139533312u64, n_hash(20u8));
/// assert_eq!(52776558133248u64, n_hash(21u8));
/// assert_eq!(211106232532992u64, n_hash(22u8));
/// assert_eq!(844424930131968u64, n_hash(23u8));
/// assert_eq!(3377699720527872u64, n_hash(24u8));
/// assert_eq!(13510798882111488u64, n_hash(25u8));
/// assert_eq!(54043195528445952u64, n_hash(26u8));
/// assert_eq!(216172782113783808u64, n_hash(27u8));
/// assert_eq!(864691128455135232u64, n_hash(28u8));
/// assert_eq!(3458764513820540928u64, n_hash(29u8));
/// ```
///
#[inline]
pub const fn n_hash(depth: u8) -> u64 {
    assert_valid_depth(depth);
    unchecked::n_hash(depth)
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
pub const fn to_uniq(depth: u8, hash: u64) -> u64 {
    assert_valid_depth(depth);
    unchecked::to_uniq(depth, hash)
}

/// Same as [to_uniq](fn.to_uniq.html), but
/// following the [IVOA](http://ivoa.net/documents/MOC/) convention.
/// It does not rely on a sentinel bit and use one less bit.
pub const fn to_uniq_ivoa(depth: u8, hash: u64) -> u64 {
    assert_valid_depth(depth);
    unchecked::to_uniq_ivoa(depth, hash)
}
