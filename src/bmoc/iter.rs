use super::Cell;
use std::iter::Copied;
use triomphe::Arc;

pub type SliceIter<'a> = Copied<std::slice::Iter<'a, u64>>;
pub type VecIter = std::vec::IntoIter<u64>;
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

#[derive(Debug, Clone)]
pub struct CellIter<I> {
    pub(super) inner: I,
    pub(super) max_depth: u8,
}
impl<I: Iterator<Item = u64>> Iterator for CellIter<I> {
    type Item = Cell;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.inner
            .nth(n)
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}
impl<I: ExactSizeIterator<Item = u64>> ExactSizeIterator for CellIter<I> {
    fn len(&self) -> usize {
        self.inner.len()
    }
}
impl<I: DoubleEndedIterator<Item = u64>> DoubleEndedIterator for CellIter<I> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner
            .next_back()
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.inner
            .nth_back(n)
            .map(|raw| Cell::decode(raw, self.max_depth))
    }
}

#[derive(Debug, Clone)]
pub struct FlatIter<I> {
    depth_max: u8,
    deep_size: usize,
    raw_val_iter: I,
    curr_val: Option<u64>,
    curr_val_max: u64,
    n_returned: usize,
}

impl<I: Iterator<Item = u64>> FlatIter<I> {
    pub(super) fn new(depth_max: u8, deep_size: usize, raw_val_iter: I) -> Self {
        let mut flat_iter = Self {
            depth_max,
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
        self.depth_max
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

#[derive(Debug, Clone)]
pub struct FlatCellIter<I> {
    depth_max: u8,
    deep_size: usize,
    raw_val_iter: I,
    curr_val: Option<Cell>,
    curr_val_max: u64,
    n_returned: usize,
}

impl<I: Iterator<Item = u64>> FlatCellIter<I> {
    pub(super) fn new(depth_max: u8, deep_size: usize, raw_val_iter: I) -> Self {
        let mut flat_iter = Self {
            depth_max,
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
        self.depth_max
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
                    depth: self.depth_max,
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
                    depth: self.depth_max,
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
