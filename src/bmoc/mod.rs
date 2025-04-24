use triomphe::Arc;

pub mod iter;
mod ops;

/// A [BMOC](Bmoc) that has mutable access to its entries and therefore, is able to do in-place operations.
///
/// This type has a marker boolean to check if it's valid. While this isn't used for safety, it's used to enforce validity.
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
pub struct MutableBmoc<const VALID: bool = true> {
    max_depth: u8,
    entries: Vec<u64>,
}
impl Bmoc for MutableBmoc<true> {
    #[inline(always)]
    fn max_depth(&self) -> u8 {
        self.max_depth
    }
    #[inline(always)]
    fn entries(&self) -> &[u64] {
        &self.entries
    }
}
impl<const VALID: bool> MutableBmoc<VALID> {
    #[inline(always)]
    pub const fn new_empty(depth: u8) -> Self {
        Self {
            max_depth: depth,
            entries: Vec::new(),
        }
    }
    pub fn new_allsky(depth: u8) -> Self {
        let mut out = MutableBmoc::<false>::new_empty(depth);
        out.push_all_unchecked(0, 0, 12, true);
        out.into_whatever_unchecked()
    }
    pub fn with_capacity(max_depth: u8, capacity: usize) -> Self {
        Self {
            max_depth,
            entries: Vec::with_capacity(capacity),
        }
    }
    pub const fn max_depth(&self) -> u8 {
        self.max_depth
    }
    pub fn entries(&self) -> &[u64] {
        &self.entries
    }
    /// Postfix form of [`std::mem::take`]
    pub fn take(&mut self) -> Self {
        std::mem::take(self)
    }
    pub const fn as_valid_unchecked(&self) -> &MutableBmoc<true> {
        unsafe { std::mem::transmute(self) }
    }
    pub const fn as_invalid(&self) -> &MutableBmoc<false> {
        unsafe { std::mem::transmute(self) }
    }
    pub const fn as_valid_unchecked_mut(&mut self) -> &mut MutableBmoc<true> {
        unsafe { std::mem::transmute(self) }
    }
    pub const fn as_invalid_unchecked_mut(&mut self) -> &mut MutableBmoc<false> {
        unsafe { std::mem::transmute(self) }
    }
    pub const fn into_invalid(self) -> MutableBmoc<false> {
        // transmute used rather than destructuring so this can be const
        unsafe { std::mem::transmute(self) }
    }
    pub const fn into_valid_unchecked(self) -> MutableBmoc<true> {
        // transmute used rather than destructuring so this can be const
        unsafe { std::mem::transmute(self) }
    }
    /// Like `into_invalid` or `into_valid_checked`, but generic
    pub const fn into_whatever_unchecked<const V2: bool>(self) -> MutableBmoc<V2> {
        unsafe { std::mem::transmute(self) }
    }
    /// Pack the cells in an invalid BMOC, making it valid.
    pub fn pack(&mut self) {
        if !VALID {
            self.entries.sort();
            // On-place pack
            let mut prev_to_index = 0_usize;
            let mut curr_to_index = self.entries.len();
            while prev_to_index != curr_to_index {
                // changes occurs
                prev_to_index = curr_to_index;
                let mut i_prev_moc = 0_usize;
                let mut i_curr_moc = 0_usize;
                while i_prev_moc < prev_to_index {
                    let mut curr_cell = self.entries[i_prev_moc];
                    i_prev_moc += 1;
                    let mut curr_cell_depth = get_depth(curr_cell, self.max_depth);
                    let mut curr_cell_hash =
                        get_hash_from_delta_depth(curr_cell, self.max_depth - curr_cell_depth);
                    // Look for the first cell of the larger cell (depth - 1)  (=> 2 last bits = 00), the cell must be FULL
                    while i_prev_moc < prev_to_index
                        && (curr_cell_depth == 0
                            || is_partial(curr_cell)
                            || is_not_first_cell_of_larger_cell(curr_cell_hash))
                    {
                        if i_curr_moc != i_prev_moc {
                            self.entries[i_curr_moc] = curr_cell;
                            i_curr_moc += 1;
                        }
                        curr_cell = self.entries[i_prev_moc];
                        i_prev_moc += 1;
                        curr_cell_depth = get_depth(curr_cell, self.max_depth);
                        curr_cell_hash =
                            get_hash_from_delta_depth(curr_cell, self.max_depth - curr_cell_depth);
                    }
                    // Look at the 3 siblings
                    if i_prev_moc + 2 < prev_to_index
                        && self.entries[i_prev_moc]
                            == encode_raw_value(
                                curr_cell_depth,
                                curr_cell_hash | 1,
                                true,
                                self.max_depth,
                            )
                        && self.entries[i_prev_moc + 1]
                            == encode_raw_value(
                                curr_cell_depth,
                                curr_cell_hash | 2,
                                true,
                                self.max_depth,
                            )
                        && self.entries[i_prev_moc + 2]
                            == encode_raw_value(
                                curr_cell_depth,
                                curr_cell_hash | 3,
                                true,
                                self.max_depth,
                            )
                    {
                        self.entries[i_curr_moc] = encode_raw_value(
                            curr_cell_depth - 1,
                            curr_cell_hash >> 2,
                            true,
                            self.max_depth,
                        );
                        i_curr_moc += 1;
                        i_prev_moc += 3;
                    } else if i_curr_moc != i_prev_moc {
                        self.entries[i_curr_moc] = curr_cell;
                        i_curr_moc += 1;
                    }
                }
                curr_to_index = i_curr_moc;
            }
        }
    }
    pub fn into_packed(mut self) -> MutableBmoc<true> {
        self.pack();
        self.into_valid_unchecked()
    }
}
impl MutableBmoc<true> {
    /// Create a [`SharedBmoc`] from this mutable one.
    pub fn into_shared(self) -> SharedBmoc {
        SharedBmoc {
            max_depth: self.max_depth,
            entries: self.entries.into(),
        }
    }
    fn low_depth_raw_val_at_lower_depth(&self, raw_value: u64, new_depth: u8) -> u64 {
        debug_assert!(get_depth(raw_value, self.max_depth) <= new_depth);
        debug_assert!(new_depth <= self.max_depth);
        let twice_delta_depth = (self.max_depth - new_depth) << 1;
        (raw_value >> twice_delta_depth) | (raw_value & 1_u64)
    }
    /// Lower the maximum depth of the current BMOC
    pub fn lower_depth(&mut self, new_depth: u8) -> &mut Self {
        assert!(
            new_depth < self.max_depth,
            "The given depth must be lower than the depth max of the BMOC"
        );
        let mut i_new = 0_usize;
        let mut prev_hash_at_new_depth = loop {
            if i_new == self.entries.len() {
                // All cells have a depth <= new_depth
                break None;
            }
            let raw_value = self.entries[i_new];
            let depth = get_depth(raw_value, self.max_depth);
            if depth <= new_depth {
                self.entries[i_new] = self.low_depth_raw_val_at_lower_depth(raw_value, new_depth);
                i_new += 1;
            } else {
                break Some(get_hash_from_delta_depth(
                    raw_value,
                    self.max_depth - new_depth,
                ));
            }
        };
        for i in (i_new + 1)..self.entries.len() {
            let raw_value = self.entries[i];
            let depth = get_depth(raw_value, self.max_depth);
            if depth <= new_depth {
                if prev_hash_at_new_depth.is_some() {
                    self.entries[i_new] = (prev_hash_at_new_depth.take().unwrap() << 2) | 2_u64;
                    i_new += 1;
                }
                self.entries[i_new] = self.low_depth_raw_val_at_lower_depth(raw_value, new_depth);
                i_new += 1;
            } else {
                let curr_hash_at_new_depth =
                    get_hash_from_delta_depth(raw_value, self.max_depth - new_depth);
                if let Some(prev_val_at_new_depth) = prev_hash_at_new_depth {
                    if prev_val_at_new_depth != curr_hash_at_new_depth {
                        self.entries[i_new] = (prev_val_at_new_depth << 2) | 2_u64; // sentinel bit + flag = 0
                        i_new += 1;
                        prev_hash_at_new_depth.replace(curr_hash_at_new_depth);
                    }
                } else {
                    prev_hash_at_new_depth.replace(curr_hash_at_new_depth);
                }
            }
        }
        if prev_hash_at_new_depth.is_some() {
            self.entries[i_new] = (prev_hash_at_new_depth.take().unwrap() << 2) | 2_u64;
            i_new += 1;
        }
        self.entries.truncate(i_new);
        self.max_depth = new_depth;
        self
    }
    #[inline(always)]
    pub fn into_cells(self) -> iter::CellIter<iter::VecIter> {
        iter::CellIter::new(self.max_depth, self.entries.into_iter())
    }
    #[inline(always)]
    pub fn into_flat_cells(self) -> iter::FlatCellIter<iter::VecIter> {
        iter::FlatCellIter::new(self.max_depth, self.deep_size(), self.entries.into_iter())
    }
    #[inline(always)]
    pub fn into_flat_iter(self) -> iter::FlatIter<iter::VecIter> {
        iter::FlatIter::new(self.max_depth, self.deep_size(), self.entries.into_iter())
    }
    #[inline(always)]
    pub fn into_ranges(self) -> iter::RangeIter<iter::VecIter> {
        iter::RangeIter::new(self.max_depth, self.entries.into_iter())
    }
    #[inline(always)]
    pub fn into_flagged_ranges(self) -> iter::FlaggedRangeIter<iter::VecIter> {
        iter::FlaggedRangeIter::new(self.max_depth, self.entries.into_iter())
    }
}
impl MutableBmoc<false> {
    pub const fn new(max_depth: u8, entries: Vec<u64>) -> Self {
        Self { max_depth, entries }
    }
    pub const fn set_max_depth(&mut self, depth: u8) -> &mut Self {
        self.max_depth = depth;
        self
    }
    pub const fn entries_mut(&mut self) -> &mut Vec<u64> {
        &mut self.entries
    }
    pub fn push_unchecked(&mut self, depth: u8, hash: u64, is_full: bool) -> &mut Self {
        self.entries
            .push(encode_raw_value(depth, hash, is_full, self.max_depth));
        self
    }
    pub fn push_all_unchecked(
        &mut self,
        depth: u8,
        hash_min: u64,
        hash_max: u64,
        is_full: bool,
    ) -> &mut Self {
        debug_assert!(hash_max >= hash_min, "Invalid hash range");
        self.entries.reserve((hash_max - hash_min) as _);
        self.entries.extend(
            (hash_min..hash_max).map(|hash| encode_raw_value(depth, hash, is_full, self.max_depth)),
        );
        self
    }
    pub fn push_raw_unchecked(&mut self, raw_value: u64) -> &mut Self {
        self.entries.push(raw_value);
        self
    }
}

/// A shared [BMOC](Bmoc) backed by an [`Arc`]. It's more memory-efficient, but can't be modified in-place.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SharedBmoc {
    pub max_depth: u8,
    pub entries: Arc<[u64]>,
}
impl SharedBmoc {
    #[inline(always)]
    pub fn new_empty(depth: u8) -> Self {
        MutableBmoc::new_empty(depth).into_shared()
    }
    #[inline(always)]
    pub fn new_allsky(depth: u8) -> Self {
        MutableBmoc::new_allsky(depth).into_shared()
    }
    #[inline(always)]
    pub fn into_cells(self) -> iter::CellIter<iter::ArcIter> {
        iter::CellIter::new(self.max_depth, iter::ArcIter::new(self.entries))
    }
    #[inline(always)]
    pub fn into_flat_cells(self) -> iter::FlatCellIter<iter::ArcIter> {
        iter::FlatCellIter::new(
            self.max_depth,
            self.deep_size(),
            iter::ArcIter::new(self.entries),
        )
    }
    #[inline(always)]
    pub fn into_flat_iter(self) -> iter::FlatIter<iter::ArcIter> {
        iter::FlatIter::new(
            self.max_depth,
            self.deep_size(),
            iter::ArcIter::new(self.entries),
        )
    }
    #[inline(always)]
    pub fn into_ranges(self) -> iter::RangeIter<iter::ArcIter> {
        iter::RangeIter::new(self.max_depth, iter::ArcIter::new(self.entries))
    }
    #[inline(always)]
    pub fn into_flagged_ranges(self) -> iter::FlaggedRangeIter<iter::ArcIter> {
        iter::FlaggedRangeIter::new(self.max_depth, iter::ArcIter::new(self.entries))
    }
}
impl Bmoc for SharedBmoc {
    #[inline(always)]
    fn max_depth(&self) -> u8 {
        self.max_depth
    }
    #[inline(always)]
    fn entries(&self) -> &[u64] {
        &self.entries
    }
}
impl From<MutableBmoc> for SharedBmoc {
    fn from(value: MutableBmoc) -> Self {
        value.into_shared()
    }
}

/// A borrowed [BMOC](Bmoc) backed by a slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BorrowedBmoc<'a> {
    pub max_depth: u8,
    pub entries: &'a [u64],
}
impl Bmoc for BorrowedBmoc<'_> {
    #[inline(always)]
    fn max_depth(&self) -> u8 {
        self.max_depth
    }
    #[inline(always)]
    fn entries(&self) -> &[u64] {
        self.entries
    }
}
impl<'a> BorrowedBmoc<'a> {
    #[inline(always)]
    pub fn new<B: Bmoc + ?Sized>(from: &'a B) -> Self {
        Self {
            max_depth: from.max_depth(),
            entries: from.entries(),
        }
    }
    #[inline(always)]
    pub fn into_cells(&self) -> iter::CellIter<iter::SliceIter<'a>> {
        iter::CellIter::new(self.max_depth, self.entries.iter().copied())
    }
    #[inline(always)]
    pub fn into_flat_cells(&self) -> iter::FlatCellIter<iter::SliceIter<'a>> {
        iter::FlatCellIter::new(
            self.max_depth,
            self.deep_size(),
            self.entries.iter().copied(),
        )
    }
    #[inline(always)]
    pub fn into_flat_iter(&self) -> iter::FlatIter<iter::SliceIter<'a>> {
        iter::FlatIter::new(
            self.max_depth,
            self.deep_size(),
            self.entries.iter().copied(),
        )
    }
    #[inline(always)]
    pub fn into_ranges(&self) -> iter::RangeIter<iter::SliceIter<'a>> {
        iter::RangeIter::new(self.max_depth, self.entries.iter().copied())
    }
    #[inline(always)]
    pub fn into_flagged_ranges(&self) -> iter::FlaggedRangeIter<iter::SliceIter<'a>> {
        iter::FlaggedRangeIter::new(self.max_depth, self.entries.iter().copied())
    }
}

/// Whether or not a point is inside the MOC
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Status {
    Out,
    Unknown,
    In,
}

pub trait Bmoc {
    fn max_depth(&self) -> u8;
    fn entries(&self) -> &[u64];

    #[inline(always)]
    fn as_borrowed(&self) -> BorrowedBmoc<'_> {
        BorrowedBmoc::new(self)
    }
    #[inline(always)]
    fn cells(&self) -> iter::CellIter<iter::SliceIter<'_>> {
        iter::CellIter::new(self.max_depth(), self.entries().iter().copied())
    }
    #[inline(always)]
    fn flat_cells(&self) -> iter::FlatCellIter<iter::SliceIter<'_>> {
        iter::FlatCellIter::new(
            self.max_depth(),
            self.deep_size(),
            self.entries().iter().copied(),
        )
    }
    #[inline(always)]
    fn flat_iter(&self) -> iter::FlatIter<iter::SliceIter<'_>> {
        iter::FlatIter::new(
            self.max_depth(),
            self.deep_size(),
            self.entries().iter().copied(),
        )
    }
    #[inline(always)]
    fn ranges(&self) -> iter::RangeIter<iter::SliceIter<'_>> {
        iter::RangeIter::new(self.max_depth(), self.entries().iter().copied())
    }
    #[inline(always)]
    fn flagged_ranges(&self) -> iter::FlaggedRangeIter<iter::SliceIter<'_>> {
        iter::FlaggedRangeIter::new(self.max_depth(), self.entries().iter().copied())
    }
    /// Get the number of entries in this BMOC.
    #[inline(always)]
    fn size(&self) -> usize {
        self.entries().len()
    }
    /// Get the number of max-depth cells in this BMOC.
    #[inline(always)]
    fn deep_size(&self) -> usize {
        let max_depth = self.max_depth();
        self.entries()
            .iter()
            .map(|&raw| {
                crate::unchecked::nside_square(max_depth - get_depth(raw, max_depth)) as usize
            })
            .sum()
    }
    /// Compare this BMOC against another one generically.
    #[inline(always)]
    fn gen_eq<B2: Bmoc + ?Sized>(&self, other: &B2) -> bool
    where
        Self: Sized,
    {
        self.max_depth() == other.max_depth() && self.entries() == other.entries()
    }
    /// Test the given point and return its "Status": in, out of the MOC or maybe.
    #[inline(always)]
    fn test_coo(&self, lon: f64, lat: f64) -> Status {
        test_coo(self.max_depth(), self.entries(), lon, lat)
    }
    /// Test the given cell and return whether it's inside or outside of the MOC.
    ///
    /// TODO: test
    #[inline(always)]
    fn test_cell(&self, depth: u8, hash: u64) -> Status {
        test_cell(self.max_depth(), self.entries(), depth, hash)
    }

    /// Returns the BMOC complement:
    /// - cells with flag set to 1 (fully covered) are removed
    /// - cells with flag set to 0 (partially covered) are kept
    /// - empty cells are added with flag set to 1
    ///
    /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
    #[inline(always)]
    fn not(&self) -> MutableBmoc {
        ops::not(self.max_depth(), self.entries())
    }

    /// See [`and`](Bmoc::and).
    #[inline(always)]
    fn and_dyn(&self, rhs: BorrowedBmoc) -> MutableBmoc {
        ops::and(self.as_borrowed(), rhs)
    }

    /// Returns the intersection of this BMOC with the given BMOC:
    /// - all non overlapping cells are removed
    /// - when two cells are overlapping, the overlapping part is kept
    ///   - the value of the flag is the result of a logical AND between the flags of the merged cells.
    ///
    /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
    #[inline(always)]
    fn and<B2: Bmoc + ?Sized>(&self, rhs: &B2) -> MutableBmoc
    where
        Self: Sized,
    {
        ops::and(self.as_borrowed(), rhs.as_borrowed())
    }

    /// See [`or`](Bmoc::or).
    #[inline(always)]
    fn or_dyn(&self, rhs: BorrowedBmoc) -> MutableBmoc {
        ops::or(self.as_borrowed(), rhs)
    }

    /// Returns the union of this BMOC with the given BMOC:
    /// - all non overlapping cells in both BMOCs are kept;
    /// - overlapping cells are merged, the value of the flag is the result of a logical OR between.
    ///
    /// the flags of the merged cells.
    /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
    fn or<B2: Bmoc + ?Sized>(&self, rhs: &B2) -> MutableBmoc
    where
        Self: Sized,
    {
        ops::or(self.as_borrowed(), rhs.as_borrowed())
    }

    /// See [`xor`](Bmoc::xor).
    #[inline(always)]
    fn xor_dyn(&self, rhs: BorrowedBmoc) -> MutableBmoc {
        ops::xor(self.as_borrowed(), rhs)
    }

    /// Returns the symmetric difference of this BMOC with the given BMOC:
    /// - all non overlapping cells in both BMOCs are kept
    /// - when two cells are overlapping, the overlapping part is:
    ///   - removed if both flags = 1
    ///   - kept if one of the flags = 0 (since 0 meas partially covered but O don't know which part)
    ///
    /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
    fn xor<B2: Bmoc + ?Sized>(&self, rhs: &B2) -> MutableBmoc
    where
        Self: Sized,
    {
        ops::xor(self.as_borrowed(), rhs.as_borrowed())
    }

    /// See [`minus`](Bmoc::minus).
    #[inline(always)]
    fn minus_dyn(&self, rhs: BorrowedBmoc) -> MutableBmoc {
        ops::minus(self.as_borrowed(), rhs)
    }

    /// Returns the difference of this BMOC (left) with the given BMOC (right):
    /// - all non overlapping cells of this (left) BMOC are kept
    /// - non overlapping cells of the other (right) BMOC are removed
    ///   if full, and kept if partially covered (since A MINUS B = A AND (NOT(B))
    /// - when two cells are overlapping, the overlapping part is:
    ///   - removed if both flags = 1
    ///   - kept if one of the flags = 0 (since 0 meas partially covered but O don't know which part)
    ///
    /// Poor's man implementation: A MINUS B = A AND NOT(B).
    fn minus<B2: Bmoc + ?Sized>(&self, rhs: &B2) -> MutableBmoc
    where
        Self: Sized,
    {
        ops::minus(self.as_borrowed(), rhs.as_borrowed())
    }
}
impl<T: Bmoc> Bmoc for &T {
    fn max_depth(&self) -> u8 {
        T::max_depth(self)
    }
    fn entries(&self) -> &[u64] {
        T::entries(self)
    }
}
impl PartialEq for dyn Bmoc {
    fn eq(&self, other: &Self) -> bool {
        self.max_depth() == other.max_depth() && self.entries() == other.entries()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Cell {
    pub raw_value: u64,
    pub depth: u8,
    pub hash: u64,
    pub is_full: bool,
}

impl Cell {
    pub const fn decode(raw_value: u64, max_depth: u8) -> Cell {
        // Extract the flag
        let is_full = (raw_value & 1_u64) == 1_u64;
        // Remove the flag bit, then divide by 2 (2 bits per level)
        let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
        // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
        let hash = raw_value >> (2 + (delta_depth << 1));
        let depth = max_depth - delta_depth;
        Cell {
            raw_value,
            depth,
            hash,
            is_full,
        }
    }
}

/// Create a BMOC raw value coding the depth, the hash and a flag in a way such that
/// the natural ordering follow a z-order curve.
///
/// # Inputs
/// - `depth`: depth of the hash value
/// - `hash`: hash value
/// - `is_full`: must be `false` (not full) or `true` (full)
/// - `max_depth`: the depth of the BMOC (we can use 29 for a unique raw value, but it will work
///   only with languages supporting unsigned 64 bit integers)
///
/// # Outputs
/// - the value coded like this:
///   - BBBBxx...xxS00...00F if depth < max_depth
///   - BBBBxx...xxxx...xxSF if depth = depht_max
///   - with in bith cases:
///     -  B: the 4 bits coding the base hash [0- 11]
///     - xx: the 2 bits of level x
///     -  S: the sentinel bit coding the depth
///     - 00: if (depth != depht_max) those bits are unused bits
///     -  F: the flag bit (0: partial, 1: full)
#[inline]
pub const fn encode_raw_value(depth: u8, hash: u64, is_full: bool, max_depth: u8) -> u64 {
    // Set the sentinel bit
    let mut hash = (hash << 1) | 1_u64;
    // Shift according to the depth and add space for the flag bit
    hash <<= 1 + ((max_depth - depth) << 1);
    // Set the flag bit if needed
    hash | (is_full as u64)
}

#[inline]
const fn rm_flag(raw_value: u64) -> u64 {
    raw_value >> 1
}

#[inline]
const fn is_partial(raw_value: u64) -> bool {
    (raw_value & 1_u64) == 0_u64
}

#[inline]
const fn is_not_first_cell_of_larger_cell(hash: u64) -> bool {
    (hash & 3_u64) != 0_u64
}

#[inline]
const fn get_depth(raw_value: u64, max_depth: u8) -> u8 {
    get_depth_no_flag(rm_flag(raw_value), max_depth)
}

#[inline]
const fn get_depth_no_flag(raw_value_no_flag: u64, max_depth: u8) -> u8 {
    max_depth - (raw_value_no_flag.trailing_zeros() >> 1) as u8
}

#[inline]
const fn get_hash_from_delta_depth(raw_value: u64, delta_depth: u8) -> u64 {
    raw_value >> (2 + (delta_depth << 1))
}

/// Returns `true` if the given high resolution cell is in the low resolution cell
#[inline]
fn is_in(low_resolution: &Cell, high_resolution: &Cell) -> bool {
    low_resolution.depth <= high_resolution.depth
        && low_resolution.hash
            == (high_resolution.hash >> ((high_resolution.depth - low_resolution.depth) << 1))
}

// This is implemented as a free function to avoid monomorphization bloat
#[inline(never)]
fn test_coo(max_depth: u8, entries: &[u64], lon: f64, lat: f64) -> Status {
    let h_raw = encode_raw_value(
        max_depth,
        crate::Layer::get(max_depth).hash(lon, lat),
        true,
        max_depth,
    );
    match entries.binary_search(&h_raw) {
        Ok(i) => {
            if is_partial(entries[i]) {
                Status::Unknown
            } else {
                Status::In
            }
        }
        Err(i) => {
            let cell = Cell::decode(h_raw, max_depth);
            // look in next or previous cels
            if i > 0 && is_in(&Cell::decode(entries[i - 1], max_depth), &cell) {
                if is_partial(entries[i - 1]) {
                    Status::Unknown
                } else {
                    Status::In
                }
            } else if i < entries.len() && is_in(&Cell::decode(entries[i], max_depth), &cell) {
                if is_partial(entries[i]) {
                    Status::Unknown
                } else {
                    Status::In
                }
            } else {
                Status::Out
            }
        }
    }
}

#[inline(never)]
fn test_cell(max_depth: u8, entries: &[u64], depth: u8, hash: u64) -> Status {
    let h_raw = encode_raw_value(depth, hash, true, max_depth);
    match entries.binary_search(&h_raw) {
        Ok(i) => {
            if is_partial(entries[i]) {
                Status::Unknown
            } else {
                Status::In
            }
        }
        Err(i) => {
            let cell = Cell::decode(h_raw, max_depth);
            // look in next or previous cels
            if i > 0 && is_in(&Cell::decode(entries[i - 1], max_depth), &cell) {
                if is_partial(entries[i - 1]) {
                    Status::Unknown
                } else {
                    Status::In
                }
            } else if i < entries.len() && is_in(&Cell::decode(entries[i], max_depth), &cell) {
                if is_partial(entries[i]) {
                    Status::Unknown
                } else {
                    Status::In
                }
            } else {
                Status::Out
            }
        }
    }
}
