use super::*;
use std::array::IntoIter;
use std::fmt::{self, Debug, Formatter};
use std::iter::{Enumerate, FilterMap};
use std::marker::PhantomData;
use std::slice::{Iter, IterMut};

/// A map of enum keys to values, with an API similar to the standard [`HashMap`](std::collections::HashMap).
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EnumMap<K, V, const N: usize> {
    values: [Option<V>; N],
    _marker: PhantomData<[K; N]>,
}
impl<K, V, const N: usize> EnumMap<K, V, N> {
    pub const fn new() -> Self {
        Self {
            values: [const { None }; N],
            _marker: PhantomData,
        }
    }
}
impl<K: Enum, V, const N: usize> EnumMap<K, V, N> {
    pub fn iter(&self) -> MapIter<'_, K, V> {
        self.values
            .iter()
            .enumerate()
            .filter_map(|(n, v)| v.as_ref().map(|v| (K::from_enum_index(n), v)))
    }
    pub fn iter_mut(&mut self) -> MapIterMut<'_, K, V> {
        self.values
            .iter_mut()
            .enumerate()
            .filter_map(|(n, v)| v.as_mut().map(|v| (K::from_enum_index(n), v)))
    }
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        self.values[key.into_enum_index()].replace(value)
    }
    pub fn remove(&mut self, key: K) -> Option<V> {
        self.values[key.into_enum_index()].take()
    }
    pub fn get(&mut self, key: K) -> Option<&V> {
        self.values[key.into_enum_index()].as_ref()
    }
    pub fn get_mut(&mut self, key: K) -> Option<&mut V> {
        self.values[key.into_enum_index()].as_mut()
    }
}
impl<K, V, const N: usize> Default for EnumMap<K, V, N> {
    fn default() -> Self {
        Self::new()
    }
}
impl<K, V: Clone, const N: usize> Clone for EnumMap<K, V, N> {
    fn clone(&self) -> Self {
        Self {
            values: self.values.clone(),
            _marker: PhantomData,
        }
    }
    fn clone_from(&mut self, source: &Self) {
        self.values.clone_from(&source.values);
    }
}
impl<K, V: Copy, const N: usize> Copy for EnumMap<K, V, N> {}
impl<K: Debug + Enum, V: Debug, const N: usize> Debug for EnumMap<K, V, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_map().entries(self.iter()).finish()
    }
}
impl<K: Enum, V, const N: usize> IntoIterator for EnumMap<K, V, N> {
    type IntoIter = MapIntoIter<K, V, N>;
    type Item = (K, V);

    fn into_iter(self) -> Self::IntoIter {
        self.values
            .into_iter()
            .enumerate()
            .filter_map(|(n, v)| v.map(|v| (K::from_enum_index(n), v)))
    }
}
impl<'a, K: Enum, V, const N: usize> IntoIterator for &'a EnumMap<K, V, N> {
    type IntoIter = MapIter<'a, K, V>;
    type Item = (K, &'a V);

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
impl<'a, K: Enum, V, const N: usize> IntoIterator for &'a mut EnumMap<K, V, N> {
    type IntoIter = MapIterMut<'a, K, V>;
    type Item = (K, &'a mut V);

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

pub type MapIntoIter<K, V, const N: usize> =
    FilterMap<Enumerate<IntoIter<Option<V>, N>>, fn((usize, Option<V>)) -> Option<(K, V)>>;
pub type MapIter<'a, K, V> =
    FilterMap<Enumerate<Iter<'a, Option<V>>>, fn((usize, &'a Option<V>)) -> Option<(K, &'a V)>>;
pub type MapIterMut<'a, K, V> = FilterMap<
    Enumerate<IterMut<'a, Option<V>>>,
    fn((usize, &'a mut Option<V>)) -> Option<(K, &'a mut V)>,
>;

pub type EnumMap4<K, V> = EnumMap<K, V, 4>;
pub type EnumMap8<K, V> = EnumMap<K, V, 8>;
pub type CardinalMap<V> = EnumMap4<Cardinal, V>;
pub type OrdinalMap<V> = EnumMap4<Ordinal, V>;
pub type DirectionMap<V> = EnumMap8<Direction, V>;
