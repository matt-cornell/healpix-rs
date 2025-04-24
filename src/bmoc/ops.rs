use super::*;
use std::cmp::{Ordering, max};

type CellIter<'a> = iter::CellIter<iter::SliceIter<'a>>;

#[inline(never)]
pub fn not(max_depth: u8, entries: &[u64]) -> MutableBmoc {
    // Worst case: only 1 sub-cell by cell in the MOC (+11 for depth 0)
    let mut builder = MutableBmoc::<false>::with_capacity(max_depth, 3 * entries.len() + 12);
    // Empty MOC, easy
    if entries.is_empty() {
        return builder
            .push_all_unchecked(0_u8, 0, 12, true)
            .take()
            .into_valid_unchecked();
    }
    // Real case
    let mut d = 0_u8;
    let mut h = 0_u64;
    // Go down to first cell
    let mut cell = Cell::decode(entries[0], max_depth);
    go_down(&mut d, &mut h, cell.depth, cell.hash, true, &mut builder);
    if !cell.is_full {
        builder.push_raw_unchecked(cell.raw_value);
    }
    // Between first and last
    for &entry in &entries[1..] {
        cell = Cell::decode(entry, max_depth);
        let dd = dd_4_go_up(d, h, cell.depth, cell.hash);
        go_up(&mut d, &mut h, dd, true, &mut builder);
        go_down(&mut d, &mut h, cell.depth, cell.hash, true, &mut builder);
        if !cell.is_full {
            builder.push_raw_unchecked(cell.raw_value);
        }
    }
    // After last
    let delta_depth = d;
    go_up(&mut d, &mut h, delta_depth, true, &mut builder); // go up to depth 0
    builder.push_all_unchecked(0, 0, 12, true);
    builder.into_valid_unchecked()
}

#[inline(never)]
pub fn and(lhs: BorrowedBmoc, rhs: BorrowedBmoc) -> MutableBmoc {
    let mut builder = MutableBmoc::with_capacity(
        max(lhs.max_depth, rhs.max_depth),
        max(lhs.entries.len(), rhs.entries.len()),
    );
    let mut it_left = lhs.cells();
    let mut it_right = rhs.cells();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
        match l.depth.cmp(&r.depth) {
            Ordering::Less => {
                let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
                match l.hash.cmp(&hr_at_dl) {
                    Ordering::Less => left = it_left.next(),
                    Ordering::Greater => right = it_right.next(),
                    Ordering::Equal => {
                        debug_assert_eq!(l.hash, hr_at_dl);
                        builder.push_unchecked(r.depth, r.hash, r.is_full && l.is_full);
                        right = it_right.next()
                    }
                }
            }
            Ordering::Greater => {
                let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
                match hl_at_dr.cmp(&r.hash) {
                    Ordering::Less => left = it_left.next(),
                    Ordering::Greater => right = it_right.next(),
                    Ordering::Equal => {
                        debug_assert_eq!(hl_at_dr, r.hash);
                        builder.push_unchecked(l.depth, l.hash, r.is_full && l.is_full);
                        left = it_left.next()
                    }
                }
            }
            Ordering::Equal => {
                debug_assert_eq!(l.depth, r.depth);
                match l.hash.cmp(&r.hash) {
                    Ordering::Less => left = it_left.next(),
                    Ordering::Greater => right = it_right.next(),
                    Ordering::Equal => {
                        debug_assert_eq!(l.hash, r.hash);
                        builder.push_unchecked(l.depth, l.hash, r.is_full && l.is_full);
                        left = it_left.next();
                        right = it_right.next()
                    }
                }
            }
        }
    }
    builder.into_valid_unchecked()
}

#[inline(never)]
pub fn or(lhs: BorrowedBmoc, rhs: BorrowedBmoc) -> MutableBmoc {
    let mut builder = MutableBmoc::with_capacity(
        max(lhs.max_depth, rhs.max_depth),
        max(lhs.entries.len(), rhs.entries.len()),
    );
    let mut it_left = lhs.cells();
    let mut it_right = rhs.cells();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
        match l.depth.cmp(&r.depth) {
            Ordering::Less => {
                let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
                if l.hash < hr_at_dl {
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    left = it_left.next();
                } else if l.hash > hr_at_dl {
                    builder.push_unchecked(r.depth, r.hash, r.is_full);
                    right = it_right.next();
                } else if l.is_full {
                    debug_assert_eq!(l.hash, hr_at_dl);
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    right = consume_while_overlapped(l, &mut it_right);
                    left = it_left.next();
                } else {
                    debug_assert_eq!(l.hash, hr_at_dl);
                    debug_assert!(!l.is_full);
                    let mut is_overlapped = false;
                    right =
                        consume_while_overlapped_and_partial(l, &mut it_right, &mut is_overlapped);
                    if is_overlapped {
                        right = not_in_cell_4_or(l, right.unwrap(), &mut it_right, &mut builder);
                    } else {
                        // all flags set to 0 => put large cell with flag  = 0
                        builder.push_unchecked(l.depth, l.hash, false);
                    }
                    left = it_left.next();
                }
            }
            Ordering::Greater => {
                let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
                if hl_at_dr < r.hash {
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    left = it_left.next();
                } else if hl_at_dr > r.hash {
                    builder.push_unchecked(r.depth, r.hash, r.is_full);
                    right = it_right.next();
                } else if r.is_full {
                    debug_assert_eq!(hl_at_dr, r.hash);
                    builder.push_unchecked(r.depth, r.hash, r.is_full);
                    left = consume_while_overlapped(r, &mut it_left);
                    right = it_right.next();
                } else {
                    debug_assert_eq!(hl_at_dr, r.hash);
                    debug_assert!(!r.is_full);
                    let mut is_overlapped = false;
                    left =
                        consume_while_overlapped_and_partial(r, &mut it_left, &mut is_overlapped);
                    if is_overlapped {
                        left = not_in_cell_4_or(r, left.unwrap(), &mut it_left, &mut builder);
                    } else {
                        // all flags set to 0 => put large cell with flag  = 0
                        builder.push_unchecked(r.depth, r.hash, false);
                    }
                    right = it_right.next();
                }
            }
            Ordering::Equal => {
                debug_assert_eq!(l.depth, r.depth);
                match l.hash.cmp(&r.hash) {
                    Ordering::Less => {
                        builder.push_unchecked(l.depth, l.hash, l.is_full);
                        left = it_left.next();
                    }
                    Ordering::Greater => {
                        builder.push_unchecked(r.depth, r.hash, r.is_full);
                        right = it_right.next();
                    }
                    Ordering::Equal => {
                        debug_assert_eq!(l.hash, r.hash);
                        builder.push_unchecked(l.depth, l.hash, r.is_full || l.is_full);
                        left = it_left.next();
                        right = it_right.next();
                    }
                }
            }
        }
    }
    while let Some(l) = &left {
        debug_assert!(right.is_none());
        builder.push_unchecked(l.depth, l.hash, l.is_full);
        left = it_left.next();
    }
    while let Some(r) = &right {
        debug_assert!(left.is_none());
        builder.push_unchecked(r.depth, r.hash, r.is_full);
        right = it_right.next();
    }
    builder.pack();
    builder.into_valid_unchecked()
}

fn not_in_cell_4_or(
    low_resolution: &Cell,
    mut c: Cell,
    iter: &mut CellIter,
    builder: &mut MutableBmoc<false>,
) -> Option<Cell> {
    let mut d = low_resolution.depth;
    let mut h = low_resolution.hash;
    debug_assert!(c.is_full);
    go_down(&mut d, &mut h, c.depth, c.hash, false, builder);
    builder.push_unchecked(c.depth, c.hash, true);
    let mut is_overlapped = false;
    let mut cell;
    while {
        cell = consume_while_overlapped_and_partial(low_resolution, iter, &mut is_overlapped);
        is_overlapped
    } {
        c = cell.unwrap(); // if flag => right is not None
        let dd = dd_4_go_up(d, h, c.depth, c.hash);
        go_up(&mut d, &mut h, dd, false, builder);
        go_down(&mut d, &mut h, c.depth, c.hash, false, builder);
        builder.push_unchecked(c.depth, c.hash, true);
    }
    let dd = d - low_resolution.depth;
    go_up(&mut d, &mut h, dd, false, builder);
    go_down(
        &mut d,
        &mut h,
        low_resolution.depth,
        low_resolution.hash + 1,
        false,
        builder,
    );
    cell
}

#[inline(never)]
pub fn xor(lhs: BorrowedBmoc, rhs: BorrowedBmoc) -> MutableBmoc {
    let mut builder = MutableBmoc::with_capacity(
        max(lhs.max_depth, rhs.max_depth),
        max(lhs.entries.len(), rhs.entries.len()),
    );
    let mut it_left = lhs.cells();
    let mut it_right = rhs.cells();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
        match l.depth.cmp(&r.depth) {
            Ordering::Less => {
                let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
                if l.hash < hr_at_dl {
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    left = it_left.next();
                } else if l.hash > hr_at_dl {
                    builder.push_unchecked(r.depth, r.hash, r.is_full);
                    right = it_right.next();
                } else if l.is_full {
                    debug_assert_eq!(l.hash, hr_at_dl);
                    right = not_in_cell_4_xor(l, r, &mut it_right, &mut builder);
                    left = it_left.next();
                } else {
                    debug_assert_eq!(l.hash, hr_at_dl);
                    debug_assert!(!l.is_full);
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    right = consume_while_overlapped(l, &mut it_right);
                    left = it_left.next();
                }
            }
            Ordering::Greater => {
                let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
                if hl_at_dr < r.hash {
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    left = it_left.next();
                } else if hl_at_dr > r.hash {
                    builder.push_unchecked(r.depth, r.hash, r.is_full);
                    right = it_right.next();
                } else if r.is_full {
                    debug_assert_eq!(hl_at_dr, r.hash);
                    left = not_in_cell_4_xor(r, l, &mut it_left, &mut builder);
                    right = it_right.next();
                } else {
                    debug_assert_eq!(hl_at_dr, r.hash);
                    debug_assert!(!r.is_full);
                    builder.push_unchecked(r.depth, r.hash, r.is_full);
                    left = consume_while_overlapped(r, &mut it_left);
                    right = it_right.next();
                }
            }
            Ordering::Equal => {
                debug_assert_eq!(l.depth, r.depth);
                match l.hash.cmp(&r.hash) {
                    Ordering::Less => {
                        builder.push_unchecked(l.depth, l.hash, l.is_full);
                        left = it_left.next();
                    }
                    Ordering::Greater => {
                        builder.push_unchecked(r.depth, r.hash, r.is_full);
                        right = it_right.next();
                    }
                    Ordering::Equal => {
                        debug_assert_eq!(l.hash, r.hash);
                        let both_fully_covered = r.is_full && l.is_full;
                        if !both_fully_covered {
                            builder.push_unchecked(l.depth, l.hash, both_fully_covered);
                        }
                        left = it_left.next();
                        right = it_right.next();
                    }
                }
            }
        }
    }
    while let Some(l) = &left {
        debug_assert!(right.is_none());
        builder.push_unchecked(l.depth, l.hash, l.is_full);
        left = it_left.next();
    }
    while let Some(r) = &right {
        debug_assert!(left.is_none());
        builder.push_unchecked(r.depth, r.hash, r.is_full);
        right = it_right.next();
    }
    builder.pack();
    builder.into_valid_unchecked()
}

// add elements of the low resolution cell which are not in the c cell
fn not_in_cell_4_xor(
    low_resolution: &Cell,
    c: &Cell,
    iter: &mut CellIter,
    builder: &mut MutableBmoc<false>,
) -> Option<Cell> {
    let mut d = low_resolution.depth;
    let mut h = low_resolution.hash;
    go_down(&mut d, &mut h, c.depth, c.hash, true, builder);
    if !c.is_full {
        builder.push_unchecked(c.depth, c.hash, false);
    }
    let mut cell = iter.next();
    while let Some(c) = &cell {
        if !is_in(low_resolution, c) {
            break;
        }
        let dd = dd_4_go_up(d, h, c.depth, c.hash);
        go_up(&mut d, &mut h, dd, true, builder);
        go_down(&mut d, &mut h, c.depth, c.hash, true, builder);
        if !c.is_full {
            builder.push_unchecked(c.depth, c.hash, false);
        }
        cell = iter.next()
    }
    let dd = d - low_resolution.depth;
    go_up(&mut d, &mut h, dd, true, builder);
    go_down(
        &mut d,
        &mut h,
        low_resolution.depth,
        low_resolution.hash + 1,
        true,
        builder,
    );
    cell
}

#[inline(never)]
pub fn minus(lhs: BorrowedBmoc, rhs: BorrowedBmoc) -> MutableBmoc {
    let mut builder = MutableBmoc::with_capacity(
        max(lhs.max_depth, rhs.max_depth),
        max(lhs.entries.len(), rhs.entries.len()),
    );
    let mut it_left = lhs.cells();
    let mut it_right = rhs.cells();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
        match l.depth.cmp(&r.depth) {
            Ordering::Less => {
                // The l cell is larger than the r cell
                // - degrade r cell at the l cell depth
                let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
                if l.hash < hr_at_dl {
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    left = it_left.next();
                } else if l.hash > hr_at_dl {
                    right = it_right.next();
                } else if l.is_full {
                    debug_assert_eq!(l.hash, hr_at_dl);
                    // add elements of the l cell which are not in common with the r cell
                    right = not_in_cell_4_xor(l, r, &mut it_right, &mut builder);
                    left = it_left.next();
                } else {
                    debug_assert_eq!(l.hash, hr_at_dl);
                    debug_assert!(!l.is_full);
                    builder.push_unchecked(l.depth, l.hash, false);
                    right = consume_while_overlapped(l, &mut it_right);
                    left = it_left.next();
                }
            }
            Ordering::Greater => {
                // The r cell is larger than the l cell
                // - degrade l cell at the r cell depth
                let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
                if hl_at_dr < r.hash {
                    builder.push_unchecked(l.depth, l.hash, l.is_full);
                    left = it_left.next();
                } else if hl_at_dr > r.hash {
                    // remove the r cell
                    right = it_right.next();
                } else if !r.is_full || !l.is_full {
                    debug_assert_eq!(hl_at_dr, r.hash);
                    builder.push_unchecked(l.depth, l.hash, false);
                    left = it_left.next();
                }
            }
            Ordering::Equal => {
                debug_assert_eq!(l.depth, r.depth);
                match l.hash.cmp(&r.hash) {
                    Ordering::Less => {
                        builder.push_unchecked(l.depth, l.hash, l.is_full);
                        left = it_left.next();
                    }
                    Ordering::Greater => {
                        right = it_right.next();
                    }
                    Ordering::Equal => {
                        debug_assert_eq!(l.hash, r.hash);
                        let both_fully_covered = r.is_full && l.is_full;
                        if !both_fully_covered {
                            builder.push_unchecked(l.depth, l.hash, both_fully_covered);
                        }
                        left = it_left.next();
                        right = it_right.next();
                    }
                }
            }
        }
    }
    while let Some(l) = &left {
        debug_assert!(right.is_none());
        builder.push_unchecked(l.depth, l.hash, l.is_full);
        left = it_left.next();
    }
    builder.pack();
    builder.into_valid_unchecked()
}

#[inline]
fn consume_while_overlapped(low_resolution: &Cell, iter: &mut CellIter) -> Option<Cell> {
    let mut cell = iter.next();
    while {
        match &cell {
            Some(c) => is_in(low_resolution, c),
            None => false,
        }
    } {
        cell = iter.next();
    }
    cell
}

/// Returns boolean:
/// - false = returned cell do not overlap any more
/// - true =  returned cell overlap and its flag is 'full'
#[inline]
fn consume_while_overlapped_and_partial(
    low_resolution: &Cell,
    iter: &mut CellIter,
    res_is_overlapped: &mut bool,
) -> Option<Cell> {
    let mut cell = iter.next();
    while {
        match &cell {
            Some(c) => {
                if is_in(low_resolution, c) {
                    if c.is_full {
                        *res_is_overlapped = true;
                        false
                    } else {
                        true
                    }
                } else {
                    false
                }
            }
            None => false,
        }
    } {
        cell = iter.next();
    }
    cell
}

/// Fill with all cells from `start_hash` at `start_depth` to `start_hash_at_target_depth + 1`.
/// with `target_depth` = `start_depth - delta_depth`.
/// - `flag`: value of the is_full flag to be set in cells while going up
///
/// The output depth is the input depth minus delta_depth
/// The output hash value is the input hash at the output depth, plus one
fn go_up(
    start_depth: &mut u8,
    start_hash: &mut u64,
    delta_depth: u8,
    flag: bool,
    builder: &mut MutableBmoc<false>,
) {
    // let output_depth = *start_depth - delta_depth;       // For debug only
    // let output_hash = (*start_hash >> (delta_depth << 1)) + 1; // For debug only
    for _ in 0_u8..delta_depth {
        let target_hash = *start_hash | 3_u64;
        builder.push_all_unchecked(*start_depth, *start_hash + 1, target_hash + 1, flag);
        *start_hash >>= 2;
        *start_depth -= 1;
    }
    *start_hash += 1;
    // debug_assert_eq!(*start_depth, output_depth);
    // debug_assert_eq!(*start_hash, output_hash);
}

fn go_down(
    start_depth: &mut u8,
    start_hash: &mut u64,
    target_depth: u8,
    target_hash: u64,
    flag: bool,
    builder: &mut MutableBmoc<false>,
) {
    debug_assert!(target_depth >= *start_depth);
    let mut twice_dd = (target_depth - *start_depth) << 1;
    for d in *start_depth..=target_depth {
        //range(0, target_depth - start_depth).rev() {
        let target_h_at_d = target_hash >> twice_dd;
        builder.push_all_unchecked(d, *start_hash, target_h_at_d, flag);
        if d != target_depth {
            *start_hash = target_h_at_d << 2;
            twice_dd -= 2;
        }
    }
    *start_depth = target_depth;
    *start_hash = target_hash;
}

#[inline]
fn dd_4_go_up(d: u8, h: u64, next_d: u8, next_h: u64) -> u8 {
    // debug_assert!(d != next_d || h != next_h);
    let target_h_at_d = if next_d < d {
        // previous hash deeper than current hash => need to go up
        next_h << ((d - next_d) << 1)
    } else {
        // current hash deeper then (or equal to) previous hash => need to go up only if current hash
        next_h >> ((next_d - d) << 1)
    };
    // - look at the difference to see if we have to go up to add lower level cells
    // We look at the depth of the deeper common cell (i.e. all most significant bits are the same)
    // With XOR (^), we only set to 1 the bits which are set to 1 in a value and 0 in the other.
    // If number of leading = 64 => the two cell are identical, WRONG :/
    // If number of leading zero = 63 or 62 => are in the same cell => dd = 0
    // If number of leading zero = 60 or 61 => dd = 1
    // We just have to add .min(d) since base cells are coded on 4 bits (not 2)
    let xor = h ^ target_h_at_d;
    if xor != 0 {
        ((63_u8 - (xor.leading_zeros() as u8)) >> 1).min(d)
    } else {
        0
    }
}
