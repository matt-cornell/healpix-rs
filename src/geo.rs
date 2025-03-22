//! Geodesic operations for converting between longitude/latitude and azimuth/distance.
use super::LonLat;

/// Get the offset of a point from a reference starting point.
pub fn relative(start: LonLat, end: LonLat) -> [f64; 2] {
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
pub fn absolute(start: LonLat, offset: [f64; 2]) -> LonLat {
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

/// Get the distance between two points. Same as `relative(start, end).length()` but faster.
pub fn distance(start: LonLat, end: LonLat) -> f64 {
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

/// Azimuth from one point to another. Same as `relative(start, end).to_angle()` but faster.
pub fn azimuth(start: LonLat, end: LonLat) -> f64 {
    let LonLat { lon: lo1, lat: la1 } = start;
    let LonLat { lon: lo2, lat: la2 } = end;
    let dlo = lo2 - lo1;
    let (dlosin, dlocos) = dlo.sin_cos();
    let (la1sin, la1cos) = la1.sin_cos();
    let (la2sin, la2cos) = la2.sin_cos();
    (dlosin * la2cos).atan2(la1cos * la2sin - la1sin * la2cos * dlocos)
}
