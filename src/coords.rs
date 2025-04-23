//! Longitude and latitude coordinate utilities

use std::f64::consts::{FRAC_PI_2, TAU};

/// A longitude-latitude pair of coordinates. Normalizes and handles conversions between `f32` and `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LonLat {
    pub lon: f64,
    pub lat: f64,
}
impl LonLat {
    pub const fn from_f64s(lon: f64, lat: f64) -> Self {
        Self { lon, lat }
    }
    pub const fn from_f32s(lon: f32, lat: f32) -> Self {
        Self {
            lon: lon as _,
            lat: lat as _,
        }
    }
    pub fn normalized(self) -> Self {
        Self {
            lon: self.lon % TAU,
            lat: self.lat.clamp(-FRAC_PI_2, FRAC_PI_2),
        }
    }
    pub fn normalize(&mut self) {
        self.lon %= TAU;
        self.lat = self.lat.clamp(-FRAC_PI_2, FRAC_PI_2);
    }
    pub fn as_f64s(self) -> [f64; 2] {
        [self.lon % TAU, self.lat.clamp(-FRAC_PI_2, FRAC_PI_2)]
    }
    pub fn as_f32s(self) -> [f32; 2] {
        use std::f32::consts::{FRAC_PI_2, TAU};
        [
            self.lon as f32 % TAU,
            (self.lat as f32).clamp(-FRAC_PI_2, FRAC_PI_2),
        ]
    }
}
impl From<(f32, f32)> for LonLat {
    fn from(value: (f32, f32)) -> Self {
        Self {
            lon: value.0 as _,
            lat: value.1 as _,
        }
    }
}
impl From<(f64, f64)> for LonLat {
    fn from(value: (f64, f64)) -> Self {
        Self {
            lon: value.0,
            lat: value.1,
        }
    }
}
impl From<LonLat> for (f32, f32) {
    fn from(value: LonLat) -> Self {
        let [a, b] = value.as_f32s();
        (a, b)
    }
}
impl From<LonLat> for (f64, f64) {
    fn from(value: LonLat) -> Self {
        let [a, b] = value.as_f64s();
        (a, b)
    }
}
impl From<[f32; 2]> for LonLat {
    fn from(value: [f32; 2]) -> Self {
        Self {
            lon: value[0] as _,
            lat: value[1] as _,
        }
    }
}
impl From<[f64; 2]> for LonLat {
    fn from(value: [f64; 2]) -> Self {
        Self {
            lon: value[0],
            lat: value[1],
        }
    }
}
impl From<LonLat> for [f32; 2] {
    fn from(value: LonLat) -> Self {
        value.as_f32s()
    }
}
impl From<LonLat> for [f64; 2] {
    fn from(value: LonLat) -> Self {
        value.as_f64s()
    }
}

/// A type that acts like a longitude-latitude pair
pub trait LonLatT {
    fn lon(&self) -> f64;
    fn lat(&self) -> f64;
}
impl LonLatT for LonLat {
    fn lon(&self) -> f64 {
        self.lon
    }
    fn lat(&self) -> f64 {
        self.lat
    }
}
impl LonLatT for (f64, f64) {
    fn lon(&self) -> f64 {
        self.0
    }
    fn lat(&self) -> f64 {
        self.1
    }
}
impl LonLatT for [f64; 2] {
    fn lon(&self) -> f64 {
        self[0]
    }
    fn lat(&self) -> f64 {
        self[1]
    }
}
impl LonLatT for (f32, f32) {
    fn lon(&self) -> f64 {
        self.0 as _
    }
    fn lat(&self) -> f64 {
        self.1 as _
    }
}
impl LonLatT for [f32; 2] {
    fn lon(&self) -> f64 {
        self[0] as _
    }
    fn lat(&self) -> f64 {
        self[1] as _
    }
}
