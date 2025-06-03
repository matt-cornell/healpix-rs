use super::Cell;
use std::iter::Copied;
use std::ops::Range;
use triomphe::Arc;

pub type SliceIter<'a> = Copied<std::slice::Iter<'a, u64>>;
pub type VecIter = std::vec::IntoIter<u64>;
/// An iterator over an `Arc<[u64]>`.
#[derive(Debug, Clone)]
pub struct ArcIter {
    inner: Arc<[u64]>,
    offset: usize,
}
impl ArcIter {
    pub const fn new(inner: Arc<[u64]>) -> Self {
        Self { inner, offset: 0 }
    }
}
impl Iterator for ArcIter {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .get(self.offset)
            .copied()
            .inspect(|_| self.offset += 1)
    }
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if let Some(&v) = self.inner.get(self.offset + n) {
            self.offset += n;
            Some(v)
        } else {
            self.offset = self.inner.len();
            None
        }
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.inner.len() - self.offset;
        (n, Some(n))
    }
}
impl ExactSizeIterator for ArcIter {
    fn len(&self) -> usize {
        self.inner.len() - self.offset
    }
}

/// An iterator adaptor that decodes cells from a raw values.
#[derive(Debug, Clone)]
pub struct CellIter<I> {
    raw_val_iter: I,
    max_depth: u8,
}
impl<I> CellIter<I> {
    pub(super) const fn new(max_depth: u8, raw_val_iter: I) -> Self {
        Self {
            max_depth,
            raw_val_iter,
        }
    }
}
impl<I: Iterator<Item = u64>> Iterator for CellIter<I> {
    type Item = Cell;

    fn next(&mut self) -> Option<Self::Item> {
        self.raw_val_iter
            .next()
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.raw_val_iter
            .nth(n)
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.raw_val_iter.size_hint()
    }
}
impl<I: ExactSizeIterator<Item = u64>> ExactSizeIterator for CellIter<I> {
    fn len(&self) -> usize {
        self.raw_val_iter.len()
    }
}
impl<I: DoubleEndedIterator<Item = u64>> DoubleEndedIterator for CellIter<I> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.raw_val_iter
            .next_back()
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.raw_val_iter
            .nth_back(n)
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
}

/// An iterator adaptor that wraps a stream of raw values iterates over max-level hashes.
#[derive(Debug, Clone)]
pub struct FlatIter<I> {
    max_depth: u8,
    deep_size: usize,
    raw_val_iter: I,
    curr_val: Option<u64>,
    curr_val_max: u64,
    n_returned: usize,
}

impl<I: Iterator<Item = u64>> FlatIter<I> {
    pub(super) fn new(max_depth: u8, deep_size: usize, raw_val_iter: I) -> Self {
        let mut flat_iter = Self {
            max_depth,
            deep_size,
            raw_val_iter,
            curr_val: None,
            curr_val_max: 0_u64,
            n_returned: 0_usize,
        };
        flat_iter.next_cell();
        flat_iter
    }

    pub fn deep_size(&self) -> usize {
        self.deep_size
    }

    pub fn depth(&self) -> u8 {
        self.max_depth
    }

    fn next_cell(&mut self) -> Option<u64> {
        match self.raw_val_iter.next() {
            None => self.curr_val.take(),
            Some(raw_value) => {
                // Remove the flag bit, then divide by 2 (2 bits per level)
                let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
                let twice_delta_depth = delta_depth << 1;
                // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
                let hash = raw_value >> (2 + twice_delta_depth);
                let val = hash << twice_delta_depth;
                self.curr_val_max = val | ((1_u64 << twice_delta_depth) - 1_u64);
                self.curr_val.replace(val)
                /*// Remove the flag bit, then divide by 2 (2 bits per level)
                let twice_delta_depth = (raw_value >> 1).trailing_zeros() as u8;
                // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
                let mask = 0xFFFFFFFFFFFFFFFC_u64 << twice_delta_depth;
                let min = raw_value & mask;
                self.curr_val_max = min | ((!mask) >> 1);
                self.curr_val.replace(min)*/
            }
        }
    }
}

impl<I: Iterator<Item = u64>> Iterator for FlatIter<I> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        if let Some(val) = self.curr_val {
            self.n_returned += 1;
            if val < self.curr_val_max {
                self.curr_val.replace(val + 1)
            } else {
                self.next_cell()
            }
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.deep_size - self.n_returned;
        (n, Some(n))
    }
}
impl<I: Iterator<Item = u64>> ExactSizeIterator for FlatIter<I> {
    fn len(&self) -> usize {
        self.deep_size - self.n_returned
    }
}

/// An iterator adaptor that wraps a stream of raw values iterates over max-level cells.
#[derive(Debug, Clone)]
pub struct FlatCellIter<I> {
    max_depth: u8,
    deep_size: usize,
    raw_val_iter: I,
    curr_val: Option<Cell>,
    curr_val_max: u64,
    n_returned: usize,
}

impl<I: Iterator<Item = u64>> FlatCellIter<I> {
    pub(super) fn new(max_depth: u8, deep_size: usize, raw_val_iter: I) -> Self {
        let mut flat_iter = Self {
            max_depth,
            deep_size,
            raw_val_iter,
            curr_val: None,
            curr_val_max: 0_u64,
            n_returned: 0_usize,
        };
        flat_iter.next_cell();
        flat_iter
    }

    pub fn deep_size(&self) -> usize {
        self.deep_size
    }

    pub fn depth(&self) -> u8 {
        self.max_depth
    }

    fn next_cell(&mut self) -> Option<Cell> {
        match self.raw_val_iter.next() {
            None => self.curr_val.take(),
            Some(raw_value) => {
                // Remove the flag bit, then divide by 2 (2 bits per level)
                let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
                let twice_delta_depth = delta_depth << 1;
                // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
                let hash = raw_value >> (2 + twice_delta_depth);
                let val = hash << twice_delta_depth;
                self.curr_val_max = val | ((1_u64 << twice_delta_depth) - 1_u64);
                self.curr_val.replace(Cell {
                    raw_value,
                    depth: self.max_depth,
                    hash: val,
                    is_full: (raw_value & 1_u64) == 1_u64,
                })
            }
        }
    }
}

impl<I: Iterator<Item = u64>> Iterator for FlatCellIter<I> {
    type Item = Cell;

    fn next(&mut self) -> Option<Cell> {
        if let Some(cell) = &self.curr_val {
            self.n_returned += 1;
            if cell.hash < self.curr_val_max {
                let new_cell = Cell {
                    raw_value: cell.raw_value,
                    depth: self.max_depth,
                    hash: cell.hash + 1,
                    is_full: cell.is_full,
                };
                self.curr_val.replace(new_cell)
            } else {
                self.next_cell()
            }
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.deep_size - self.n_returned;
        (n, Some(n))
    }
}
impl<I: Iterator<Item = u64>> ExactSizeIterator for FlatCellIter<I> {
    fn len(&self) -> usize {
        self.deep_size - self.n_returned
    }
}

/// An iterator over the ranges of cells in a BMOC.
#[derive(Debug, Clone)]
pub struct RangeIter<I> {
    max_depth: u8,
    raw_val_iter: I,
    prev_min: u64,
    prev_max: u64,
}
impl<I> RangeIter<I> {
    pub const fn new(max_depth: u8, raw_val_iter: I) -> Self {
        Self {
            max_depth,
            raw_val_iter,
            prev_min: 0,
            prev_max: 0,
        }
    }
}
impl<I: Iterator<Item = u64>> Iterator for RangeIter<I> {
    type Item = Range<u64>;

    fn next(&mut self) -> Option<Self::Item> {
        self.raw_val_iter
            .find_map(|raw| {
                let cell = Cell::decode(raw, self.max_depth);
                let mut res = None;
                if cell.depth < self.max_depth {
                    let range = crate::to_range(cell.hash, self.max_depth - cell.depth);
                    if range.start != self.prev_max {
                        if self.prev_min != self.prev_max {
                            res = Some(self.prev_min..self.prev_max);
                        }
                        self.prev_min = range.start;
                    }
                    self.prev_max = range.end;
                } else if cell.hash == self.prev_max {
                    self.prev_max += 1;
                } else {
                    if self.prev_min != self.prev_max {
                        res = Some(self.prev_min..self.prev_max);
                    }
                    self.prev_min = cell.hash;
                    self.prev_max = cell.hash + 1;
                }
                res
            })
            .or_else(|| {
                let res = (self.prev_min != self.prev_max).then_some(self.prev_min..self.prev_max);
                self.prev_min = self.prev_max;
                res
            })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (mut min, max) = self.raw_val_iter.size_hint();
        if min == 0 && self.prev_min != self.prev_max {
            min = 1;
        }
        (min, max.map(|v| v + 1))
    }
}

/// An iterator over the ranges of cells in a BMOC, along with their flag.
#[derive(Debug, Clone)]
pub struct FlaggedRangeIter<I> {
    max_depth: u8,
    raw_val_iter: I,
    prev_min: u64,
    prev_max: u64,
    prev_flag: Option<bool>,
}
impl<I> FlaggedRangeIter<I> {
    pub const fn new(max_depth: u8, raw_val_iter: I) -> Self {
        Self {
            max_depth,
            raw_val_iter,
            prev_min: 0,
            prev_max: 0,
            prev_flag: None,
        }
    }
}
impl<I: Iterator<Item = u64>> Iterator for FlaggedRangeIter<I> {
    type Item = (bool, Range<u64>);

    fn next(&mut self) -> Option<Self::Item> {
        self.raw_val_iter
            .find_map(|raw| {
                let cell = Cell::decode(raw, self.max_depth);
                let prev_flag = self.prev_flag.get_or_insert(cell.is_full);
                let mut res = None;
                if cell.depth < self.max_depth {
                    let range = crate::to_range(cell.hash, self.max_depth - cell.depth);
                    if range.start == self.prev_max
                        && (self.prev_max == 0 || cell.is_full == *prev_flag)
                    {
                        self.prev_max = range.end;
                    } else {
                        if self.prev_min != self.prev_max {
                            res = Some((*prev_flag, self.prev_min..self.prev_max));
                        }
                        self.prev_min = range.start;
                        self.prev_max = range.end;
                        *prev_flag = cell.is_full;
                    }
                } else if cell.hash == self.prev_max && cell.is_full == *prev_flag {
                    self.prev_max += 1;
                } else {
                    if self.prev_min != self.prev_max {
                        res = Some((*prev_flag, self.prev_min..self.prev_max));
                    }
                    self.prev_min = cell.hash;
                    self.prev_max = cell.hash + 1;
                    *prev_flag = cell.is_full;
                }
                res
            })
            .or_else(|| {
                let res = (self.prev_min != self.prev_max)
                    .then_some((self.prev_flag.unwrap(), self.prev_min..self.prev_max));
                self.prev_min = self.prev_max;
                res
            })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (mut min, max) = self.raw_val_iter.size_hint();
        if min == 0 && self.prev_min != self.prev_max {
            min = 1;
        }
        (min, max.map(|v| v + 1))
    }
}
