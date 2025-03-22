use super::*;
use std::array::IntoIter;
use std::fmt::{self, Debug, Formatter};
use std::iter::{Enumerate, Filter, Map};
use std::marker::PhantomData;

/// A set of enum values, with an API similar to the standard [`HashSet`](std::collections::HashSet).
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EnumSet<K, const N: usize> {
    values: [bool; N],
    _marker: PhantomData<[K; N]>,
}
impl<K, const N: usize> EnumSet<K, N> {
    pub const fn new() -> Self {
        Self {
            values: [false; N],
            _marker: PhantomData,
        }
    }
}
impl<K: Enum, const N: usize> EnumSet<K, N> {
    pub fn iter(&self) -> SetIter<K, N> {
        self.into_iter()
    }
    pub fn iter_mut(&mut self) -> SetIter<K, N> {
        self.into_iter()
    }
    pub fn insert(&mut self, key: K) -> bool {
        std::mem::replace(&mut self.values[key.into_enum_index()], true)
    }
    pub fn remove(&mut self, key: K) -> bool {
        std::mem::replace(&mut self.values[key.into_enum_index()], false)
    }
    pub fn contains(&mut self, key: K) -> bool {
        self.values[key.into_enum_index()]
    }
}
impl<K, const N: usize> Default for EnumSet<K, N> {
    fn default() -> Self {
        Self::new()
    }
}
impl<K, const N: usize> Clone for EnumSet<K, N> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<K, const N: usize> Copy for EnumSet<K, N> {}
impl<K: Debug + Enum, const N: usize> Debug for EnumSet<K, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_set().entries(self.iter()).finish()
    }
}
impl<K: Enum, const N: usize> IntoIterator for EnumSet<K, N> {
    type IntoIter = SetIter<K, N>;
    type Item = K;

    fn into_iter(self) -> Self::IntoIter {
        fn filter(p: &(usize, bool)) -> bool {
            p.1
        }
        fn map<K: Enum>(p: (usize, bool)) -> K {
            K::from_enum_index(p.0)
        }
        self.values
            .into_iter()
            .enumerate()
            .filter(filter as _)
            .map(map as _)
    }
}
impl<K: Enum, const N: usize> IntoIterator for &'_ EnumSet<K, N> {
    type IntoIter = SetIter<K, N>;
    type Item = K;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
impl<K: Enum, const N: usize> IntoIterator for &'_ mut EnumSet<K, N> {
    type IntoIter = SetIter<K, N>;
    type Item = K;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

pub type SetIter<K, const N: usize> =
    Map<Filter<Enumerate<IntoIter<bool, N>>, fn(&(usize, bool)) -> bool>, fn((usize, bool)) -> K>;

pub type EnumSet4<K> = EnumSet<K, 4>;
pub type EnumSet8<K> = EnumSet<K, 8>;
pub type CardinalSet = EnumSet4<Cardinal>;
pub type OrdinalSet = EnumSet4<Ordinal>;
pub type DirectionSet = EnumSet8<Direction>;
