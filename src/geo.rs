//! Geodesic operations for converting between longitude/latitude and azimuth/distance.
use crate::coords::*;

/// Get the offset of a point from a reference starting point.
pub fn relative(start: impl LonLatT, end: impl LonLatT) -> [f64; 2] {
    relative_impl(start.as_lonlat(), end.as_lonlat())
}
fn relative_impl(start: LonLat, end: LonLat) -> [f64; 2] {
    let LonLat { lon: lo1, lat: la1 } = start;
    let LonLat { lon: lo2, lat: la2 } = end;
    let dlo = lo2 - lo1;
    let dla = la2 - la1;
    let (dlosin, dlocos) = dlo.sin_cos();
    let (la1sin, la1cos) = la1.sin_cos();
    let (la2sin, la2cos) = la2.sin_cos();
    let azimuth = (dlosin * la2cos).atan2(la1cos * la2sin - la1sin * la2cos * dlocos);
    let a = (dla * 0.5).sin().powi(2) + la1cos * la2cos * (dlo * 0.5).sin().powi(2);
    let dist = a.sqrt().atan2((1.0 - a).sqrt());
    let (azsin, azcos) = (-azimuth).sin_cos();
    [azsin * dist * 2.0, azcos * dist * 2.0]
}

/// Get the ending point on a sphere given an offset and starting point.
pub fn absolute(start: impl LonLatT, offset: [f64; 2]) -> LonLat {
    absolute_impl(start.as_lonlat(), offset)
}
fn absolute_impl(start: LonLat, offset: [f64; 2]) -> LonLat {
    let LonLat { lon: lo1, lat: la1 } = start;
    let azimuth = -offset[0].atan2(offset[1]);
    let dist = (offset[0] * offset[0] + offset[1] * offset[1]).sqrt();
    let (lasin, lacos) = la1.sin_cos();
    let (azsin, azcos) = azimuth.sin_cos();
    let (dsin, dcos) = dist.sin_cos();
    let la2 = (lasin * dcos + lacos * dsin * azcos).asin();
    let lo2 = lo1 + (azsin * dsin * lacos).atan2(dcos - lasin * la2.sin());
    LonLat::from_f64s(lo2, la2)
}

/// Get the distance between two points, in radians.
/// ```rust
///
/// use healpix::coords::Degrees;
/// use healpix::geo::{distance, relative};
///
/// let a = Degrees(0.0, 0.0);
/// let b = Degrees(15.0, 45.0);
///
/// let [x, y] = relative(a, b);
/// let vec_dist = (x * x + y * y).sqrt();
/// let dist = distance(a, b);
/// assert!((vec_dist - dist).abs() < 1.0e-8);
/// ```
pub fn distance(start: impl LonLatT, end: impl LonLatT) -> f64 {
    distance_impl(start.as_lonlat(), end.as_lonlat())
}
fn distance_impl(start: LonLat, end: LonLat) -> f64 {
    let LonLat { lon: lo1, lat: la1 } = start;
    let LonLat { lon: lo2, lat: la2 } = end;
    let dlo = lo2 - lo1;
    let dla = la2 - la1;
    let la1cos = la1.cos();
    let la2cos = la2.cos();
    let a = (dla * 0.5).sin().powi(2) + la1cos * la2cos * (dlo * 0.5).sin().powi(2);
    let dist = a.sqrt().atan2((1.0 - a).sqrt());
    dist * 2.0
}

// TODO: figure out why this doesn't work
// /// Get the azimuth between two points, in radians.
// /// ```rust
// ///
// /// use healpix::coords::Degrees;
// /// use healpix::geo::{azimuth, relative};
// ///
// /// let a = Degrees(0.0, 0.0);
// /// let b = Degrees(45.0, 45.0);
// ///
// /// let [x, y] = relative(a, b);
// /// let vec_angle = dbg!(y.atan2(x));
// /// let angle = dbg!(azimuth(a, b));
// /// assert!(dbg!(vec_angle - angle).abs() < 1.0e-8);
// /// ```
// pub fn azimuth(start: impl LonLatT, end: impl LonLatT) -> f64 {
//     azimuth_impl(start.as_lonlat(), end.as_lonlat())
// }
// fn azimuth_impl(start: LonLat, end: LonLat) -> f64 {
//     let LonLat { lon: lo1, lat: la1 } = start;
//     let LonLat { lon: lo2, lat: la2 } = end;
//     let dlo = lo2 - lo1;
//     let (dlosin, dlocos) = dlo.sin_cos();
//     let (la1sin, la1cos) = la1.sin_cos();
//     let (la2sin, la2cos) = la2.sin_cos();
//     PI - (dlosin * la2cos).atan2(la1cos * la2sin - la1sin * la2cos * dlocos)
// }
