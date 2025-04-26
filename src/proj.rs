use super::{
    F64_BUT_SIGN_BIT_MASK, F64_SIGN_BIT_MASK, FRAC_PI_2, FRAC_PI_4, SQRT_6, TRANSITION_LATITUDE,
    TRANSITION_Z,
};

/// For each HEALPix depth, stores the smallest distance from an edge of a cell to the opposite
/// edge of the same cell. If the radius of a cone is smaller than this distance, we know that
/// it will overlap maximum 9 pixels (the pixel containing the center of the cone plus
/// the 8 neighbours).  
/// In practice, this distance if the distance between the point of coordinate
/// (0, TRANSITION_LATITUDE) and it nearest point on the Northeast edge of the
/// cell of base hash 0 and coordinates in the base hash (x=0, y=nside-1).  
/// IMPORTANT REMARK:  
/// - this value is larger than the smallest center to vertex distance
/// - this value x2 is larger than the smallest diagonal (NS or EW)
/// - this value x2 is larger than the smallest edge
/// - BUT there is no case in which the value is larger than the four center-to-vertex distance
/// - BUT there is no case in which the value x2 is larger than both diagonals
///
/// So this radius is smaller than the smaller circumcircle radius (=> no cone having the smaller
/// edge-to-opposite-edge-radius radius can contains the 4 vertices of a cell (but 3 is ok))
pub(crate) static SMALLER_EDGE2OPEDGE_DIST: [f64; 30] = [
    8.410686705679302e-1,  // depth = 0
    3.7723631722170065e-1, // depth = 1
    1.8203364957037313e-1, // depth = 2
    8.91145416330163e-2,   // depth = 3
    4.3989734509169175e-2, // depth = 4
    2.1817362566054977e-2, // depth = 5
    1.0854009694242892e-2, // depth = 6
    5.409888140793663e-3,  // depth = 7
    2.6995833266547898e-3, // depth = 8
    1.3481074874673246e-3, // depth = 9
    6.735240905806414e-4,  // depth = 10
    3.365953703015157e-4,  // depth = 11
    1.682452196838741e-4,  // depth = 12
    8.410609042173736e-5,  // depth = 13
    4.204784317861652e-5,  // depth = 14
    2.1022283297961136e-5, // depth = 15
    1.0510625670060442e-5, // depth = 16
    5.255150320257332e-6,  // depth = 17
    2.6275239729465538e-6, // depth = 18
    1.3137458638808036e-6, // depth = 19
    6.568678535571394e-7,  // depth = 20
    3.284323270983175e-7,  // depth = 21
    1.642156595517884e-7,  // depth = 22
    8.21076709163242e-8,   // depth = 23
    4.105378528139296e-8,  // depth = 24
    2.0526876713226626e-8, // depth = 25
    1.0263433216329513e-8, // depth = 26
    5.131714858175969e-9,  // depth = 27
    2.5658567623093986e-9, // depth = 28
    1.2829280665188905e-9, // depth = 29
];

/// Latitude, in the equatorial region, for which the distance from the cell center to its four
/// vertices is almost equal on the sky (i.e. the shape of the cell on the sky is close to a square).
/// The larger the depth, the better the approximation (based on differential calculus).
/// > dX = dY = 1 / nside (center to vertex distance)
/// > X = 4/pi * lon     => dX = 4/pi dlon
/// > Y = 3/2 * sin(lat) => dY = 3/2 * cos(lat) dlat
/// > dlon * cos(lat) = dlat (same distance on the sky)
/// > => cos^2(lat) = 2/3 * 4/pi
/// > => lat = arccos(sqrt(2/3 * 4/pi)) ~= 22.88 deg ~= 0.39934 rad
///
/// ```rust
/// use cdshealpix::{TRANSITION_Z, FOUR_OVER_PI, LAT_OF_SQUARE_CELL};
/// assert!(f64::abs(f64::acos(f64::sqrt(TRANSITION_Z * FOUR_OVER_PI)) - LAT_OF_SQUARE_CELL) < 1e-15_f64);
/// ```
pub const LAT_OF_SQUARE_CELL: f64 = 0.399_340_199_478_977_75_f64;

/// Returns `true` if the function [best_starting_depth](fn.best_starting_depth.html) is valid
/// for the given argument `d_max_rad`. So if `d_max_rad < ~48 deg`. `d_max_rad` is given in radians.
///
/// ```rust
/// use cdshealpix::{has_best_starting_depth};
/// use std::f64::consts::PI;
///
/// assert!(!has_best_starting_depth(PI / 3f64));
/// assert!(has_best_starting_depth(PI / 4f64));
/// ```
#[inline]
pub fn has_best_starting_depth(d_max_rad: f64) -> bool {
    d_max_rad < SMALLER_EDGE2OPEDGE_DIST[0]
}

/// Returns the the smallest depth (in `[0, 29]`) at which a shape having the given largest distance
/// from its center to a border overlaps a maximum of 9 cells (the cell containing the center of
/// the shape plus the 8 neighbouring cells).  
/// Info: internally, unrolled binary search loop on 30 pre-computed values (one by depth).
///
/// Returns -1 if the given distance is very large (> ~48deg), else returns the smallest depth
/// (in [0, 29]) at which a shape having the given largest distance from its center to a border
/// overlaps a maximum of 9 cells (the cell containing the center of the shape plus the 8
/// neighbouring cells).
///
/// # Input
/// - `d_max_rad` largest possible distance, in radians, between the center and the border of a shape
///
/// # Output
/// - `depth` = the smallest depth (in `[0, 29]`) at which a shape having the given largest distance
///   from its center to a border overlaps a maximum of 9 cells (the cell containing the center of the
///   shape plus the 8 neighbouring cells).
///
/// # Panics
/// If the given distance is very large (> ~48deg), this function is not valid since the 12 base
/// cells could be overlaped by the shape
/// (see [has_best_starting_depth](fn.has_best_starting_depth.html)). Thus it panics.
///
/// # Examples
///
/// ```rust
/// use cdshealpix::{best_starting_depth};
/// use std::f64::consts::PI;
///
/// assert_eq!(0, best_starting_depth(PI / 4f64)); // 45 deg
/// assert_eq!(5, best_starting_depth(0.0174533)); //  1 deg
/// assert_eq!(7, best_starting_depth(0.0043632)); // 15 arcmin
/// assert_eq!(9, best_starting_depth(0.0013));    // 4.469 arcmin
/// assert_eq!(15, best_starting_depth(1.454E-5)); // 3 arcsec
/// assert_eq!(20, best_starting_depth(6.5E-7));   // 0.134 arcsec
/// assert_eq!(22, best_starting_depth(9.537E-8)); // 20 mas
/// ```
#[inline]
pub fn best_starting_depth(d_max_rad: f64) -> u8 {
    // Could have used an Option
    assert!(
        d_max_rad < SMALLER_EDGE2OPEDGE_DIST[0],
        "Too large value, use first function has_best_starting_depth"
    );
    // Unrolled binary search loop
    if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[29] {
        29
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[15] {
        if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[22] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[25] {
                if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[27] {
                    if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[28] {
                        28
                    } else {
                        27
                    }
                } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[26] {
                    26
                } else {
                    25
                }
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[24] {
                24
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[23] {
                23
            } else {
                22
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[18] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[20] {
                if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[21] {
                    21
                } else {
                    20
                }
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[19] {
                19
            } else {
                18
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[17] {
            17
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[16] {
            16
        } else {
            15
        }
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[7] {
        if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[11] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[13] {
                if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[14] {
                    14
                } else {
                    13
                }
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[12] {
                12
            } else {
                11
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[9] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[10] {
                10
            } else {
                9
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[8] {
            8
        } else {
            7
        }
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[3] {
        if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[5] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[6] {
                6
            } else {
                5
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[4] {
            4
        } else {
            3
        }
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[2] {
        2
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[1] {
        1
    } else {
        0
    }
}

/// Performs the HEALPix projection: `(x, y) = proj(lon, lat)`.  
/// The chosen scale is such that: base cell vertices and center coordinates are integers;
/// the distance from a cell center to its vertices equals one.  
/// This projection is multi-purpose in the sense that if `lon` is in `[-pi, pi]`, then
/// `x` is in `[-4, 4]` and if `lon` is in `[0, 2pi]`, then `x` is in `[0, 8]`.  
/// It means that a same position on the sphere can lead to different positions in the projected
/// Euclidean plane.
///
/// Simplified projection formulae are:
///  - Equatorial region
/// ```math
/// \boxed{
///   \left\{
///     \begin{array}{lcl}
///       X & = & \alpha \times \frac{4}{\pi} \\
///       Y & = & \sin(\delta) \times \frac{3}{2}
///     \end{array}
///   \right.
/// }
/// \Rightarrow
/// \left\{
///   \begin{array}{lcl}
///     \alpha \in [0, 2\pi] & \leadsto &  X \in [0, 8] \\
///     \sin\delta \in [-\frac{2}{3}, \frac{2}{3}] & \leadsto & Y \in [-1, 1]
///   \end{array}
/// \right.
/// ```
///  - Polar caps:
/// ```math
/// \boxed{
///   \left\{
///     \begin{array}{lcl}
///       t & = & \sqrt{3(1-\sin\delta)} \\
///       X & = & (\alpha\frac{4}{\pi} - 1)t+1 \\
///       Y & = & 2 - t
///     \end{array}
///   \right.
/// }
/// \Rightarrow
/// \left\{
///   \begin{array}{l}
///     \alpha \in [0, \frac{\pi}{2}] \\
///     \sin\delta \in ]\frac{2}{3}, 1]
///   \end{array}
/// \right.
/// \leadsto
/// \begin{array}{l}
///    t \in [0, 1[  \\
///    X \in ]0, 2[ \\
///    Y \in ]1, 2]
/// \end{array}
/// ```
///
/// It is the responsibility of the caller to homogenize the result according to its needs.
/// ![Proj](https://raw.githubusercontent.com/cds-astro/cds-healpix-rust/master/resources/4doc/hpx_proj.png)
///
/// # Inputs
/// - `lon` longitude in radians, support positive and negative reasonably large values
///   (naive approach, no Cody-Waite nor Payne Hanek range reduction).
/// - `lat` latitude in radians, must be in `[-pi/2, pi/2]`
///
/// # Output
/// - `(x, y)` the projected planar Euclidean coordinates of the point
///   of given coordinates `(lon, lat)` on the unit sphere
///     - `lon` &le; `0` => `x in [-8, 0]`
///     - `lon` &ge; `0` => `x in [0, 8]`
///     - `y in [-2, 2]`
///
/// # Panics
/// If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
///
/// # Examples
/// To obtain the WCS projection (see Calabretta2007), you can write:
/// ```rust
/// use cdshealpix::proj;
/// use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};
///
/// let lon = 25.1f64;
/// let lat = 46.7f64;
///
/// let (mut x, mut y) = proj(lon.to_radians(), lat.to_radians());
/// if x > 4f64 {
///   x -= 8f64;
/// }
/// x *= FRAC_PI_4;
/// y *= FRAC_PI_4;
///
/// assert!(-PI <= x && x <= PI);
/// assert!(-FRAC_PI_2 <= y && y <= FRAC_PI_2);
/// ```
///
/// Other test example:
/// ```rust
/// use cdshealpix::{TRANSITION_LATITUDE, proj};
/// use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
///
/// let (x, y) = proj(0.0, 0.0);
/// assert_eq!(0f64, x);
/// assert_eq!(0f64, y);
///
/// assert_eq!((0.0, 1.0), proj(0.0 * FRAC_PI_2, TRANSITION_LATITUDE));
///
/// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
///     f64::sqrt((p2.0 - p1.0) * (p2.0 - p1.0) + (p2.1 - p1.1) * (p2.1 - p1.1))
/// }
/// assert!(dist((0.0, 0.0), proj(0.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((0.0, 1.0), proj(0.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((1.0, 2.0), proj(0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((2.0, 1.0), proj(1.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((3.0, 2.0), proj(1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((4.0, 1.0), proj(2.0 * FRAC_PI_2 , TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((5.0, 2.0), proj(2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((6.0, 1.0), proj(3.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((7.0, 2.0), proj(3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((0.0, 1.0), proj(4.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((0.0, 0.0), proj(4.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((0.0, -1.0), proj(0.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((1.0, -2.0), proj(0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((2.0, -1.0), proj(1.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((3.0, -2.0), proj(1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((4.0, -1.0), proj(2.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((5.0, -2.0), proj(2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((6.0, -1.0), proj(3.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((7.0, -2.0), proj(3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((0.0, -1.0), proj(4.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
///
/// assert!(dist((-0.0, 0.0), proj(-0.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((-0.0, 1.0), proj(-0.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-1.0, 2.0), proj(-0.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-2.0, 1.0), proj(-1.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-3.0, 2.0), proj(-1.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-4.0, 1.0), proj(-2.0 * FRAC_PI_2 , TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-5.0, 2.0), proj(-2.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-6.0, 1.0), proj(-3.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-7.0, 2.0), proj(-3.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-0.0, 1.0), proj(-4.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-0.0, 0.0), proj(-4.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((-0.0, -1.0), proj(-0.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-1.0, -2.0), proj(-0.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-2.0, -1.0), proj(-1.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-3.0, -2.0), proj(-1.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-4.0, -1.0), proj(-2.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-5.0, -2.0), proj(-2.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-6.0, -1.0), proj(-3.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-7.0, -2.0), proj(-3.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-0.0, -1.0), proj(-4.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// ```
#[inline]
pub fn proj(lon: f64, lat: f64) -> (f64, f64) {
    super::check_lat(lat);
    let lon = abs_sign_decompose(lon);
    let lat = abs_sign_decompose(lat);
    let x = pm1_offset_decompose(lon.abs / FRAC_PI_4);
    let mut xy = (x.pm1, lat.abs);
    if is_in_equatorial_region(lat.abs) {
        proj_cea(&mut xy);
    } else {
        proj_collignon(&mut xy);
    }
    apply_offset_and_signs(&mut xy, x.offset, lon.sign, lat.sign);
    xy
}

/// Returns the hash of the base cell given the coordinates of a points in the Euclidean projection
/// plane.
/// The purpose so far is just to test and compare both speed and precision with the 45deg
/// rotation solution.
///
/// # Input
/// - `(x, y)` the coordinates in the Euclidean projection plane, i.e. $x \in [0, 8[$
///   and $y \in [-2, 2]$
///
/// # Ouput
/// - `d0h` the hash value of the base cell (i.e. the depth 0 / nside 1 cell)
///
/// # Example
/// Simple example based on the center of each base cell.
/// ```rust
/// use cdshealpix::base_cell_from_proj_coo;
///
/// assert_eq!(base_cell_from_proj_coo(1.0,  1.0),  0);
/// assert_eq!(base_cell_from_proj_coo(3.0,  1.0),  1);
/// assert_eq!(base_cell_from_proj_coo(5.0,  1.0),  2);
/// assert_eq!(base_cell_from_proj_coo(7.0,  1.0),  3);
/// assert_eq!(base_cell_from_proj_coo(0.0,  0.0),  4);
/// assert_eq!(base_cell_from_proj_coo(2.0,  0.0),  5);
/// assert_eq!(base_cell_from_proj_coo(4.0,  0.0),  6);
/// assert_eq!(base_cell_from_proj_coo(6.0,  0.0),  7);
/// assert_eq!(base_cell_from_proj_coo(1.0, -1.0),  8);
/// assert_eq!(base_cell_from_proj_coo(3.0, -1.0),  9);
/// assert_eq!(base_cell_from_proj_coo(5.0, -1.0), 10);
/// assert_eq!(base_cell_from_proj_coo(7.0, -1.0), 11);
/// ```
pub fn base_cell_from_proj_coo(x: f64, y: f64) -> u8 {
    let mut x = 0.5 * x.rem_euclid(8.0);
    let mut y = 0.5 * (y + 3.0);
    let mut i = x as u8;
    debug_assert!(i < 4); // if can be == 4, then (x as u8) & 3
    let mut j = (y as u8) << 1;
    debug_assert!(j == 0 || j == 2 || j == 4);
    x -= i as f64;
    debug_assert!((0.0..1.0).contains(&x));
    y -= (j >> 1) as f64;
    debug_assert!((0.0..1.0).contains(&y));
    let in_northwest = (x <= y) as u8; // 1/0
    let in_southeast = (x >= 1.0 - y) as u8; // 0\1
    i += in_southeast >> in_northwest; // <=> in_southeast & (1 - in_northwest) => 0 or 1
    j += in_northwest + in_southeast;
    if j == 6 {
        j = 4;
    } // Very rare case (North pole), so few risks of branch miss-prediction
    debug_assert!(j == 2 || j == 3 || j == 4);
    ((4 - j) << 2) + i
}

/// Unproject the given HEALPix projected points.  
/// This unprojection is multi-purpose in the sense that:
///  - if input `x` in `[-8, 0[`, then output `lon` in `[-2pi, 0]`
///  - if input `x` in `[ 0, 8]`, then output `lon` in `[0, 2pi]`
///  - output `lat` always in `[-pi/2, pi/2]`
///
/// # Inputs
/// - `x` the projected coordinate along the x-axis, supports positive and negative reasonably
///   large values with a naive approach (no Cody-Waite nor Payne Hanek range reduction).
/// - `y` the projected coordinate along te x-axis, must be in `[-2, 2]`
///
/// # Output
/// -  `(lon, lat)` in radians, the position on the unit sphere whose projected coordinates are
///    the input coordinates `(x, y)`.  
///   - if `x <= 0`, then `lon` in `[-2pi, 0]`;
///   - else if `x >= 0`, the  `lon` in `[0, 2pi]`
///   - `lat` always in `[-pi/2, pi/2]`.
///
/// # Panics
/// If `y` **not in** `[-2, 2]`, this method panics.
///
/// # Examples
/// To obtain the WCS un-projection (see Calabretta2007), you can write:
/// ```rust
/// use cdshealpix::{FOUR_OVER_PI, unproj};
/// use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};
///
/// let x = 2.1f64;
/// let y = 0.36f64;
///
/// let (mut lon, mut lat) = unproj(x * FOUR_OVER_PI, y * FOUR_OVER_PI);
/// if lon < 0f64 {
///     lon += 2f64 * PI;
/// }
///
/// assert!(0f64 <= lon && lon <= 2f64 * PI);
/// assert!(-FRAC_PI_2 <= lat && lat <= FRAC_PI_2);
/// ```
///
/// Other test example:
/// ```rust
/// use cdshealpix::{TRANSITION_LATITUDE, proj, unproj};
/// use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
///
/// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
///     let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
///     let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
///     2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
/// }
///
/// let points: [(f64, f64); 40] = [
///     (0.0 * FRAC_PI_2, 0.0),
///     (1.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (2.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (3.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (4.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (4.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (0.0 * FRAC_PI_2, 0.0),
///     (1.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (2.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (3.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (4.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (4.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-0.0 * FRAC_PI_2, 0.0),
///     (-1.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-2.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-3.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-4.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-4.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-0.0 * FRAC_PI_2, 0.0),
///     (-1.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-2.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-3.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-4.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-4.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)
/// ];
///
/// for (lon, lat) in points.iter() {
///     let (x, y): (f64, f64) = proj(*lon, *lat);
///     assert!(dist((*lon, *lat), unproj(x, y)) < 1e-15);
/// }
/// ```
#[inline]
pub fn unproj(x: f64, y: f64) -> (f64, f64) {
    check_y(y);
    let x = abs_sign_decompose(x);
    let y = abs_sign_decompose(y);
    let lon = pm1_offset_decompose(x.abs);
    let mut lonlat = (lon.pm1, y.abs);
    if is_in_projected_equatorial_region(y.abs) {
        deproj_cea(&mut lonlat);
    } else {
        deproj_collignon(&mut lonlat);
    }
    apply_offset_and_signs(&mut lonlat, lon.offset, x.sign, y.sign);
    lonlat.0 *= FRAC_PI_4;
    lonlat
}

/// Verify that the projected y coordinate is in [-2, 2], panics if not.
#[inline]
fn check_y(y: f64) {
    assert!((-2f64..=2f64).contains(&y));
}

/// Returns `true` if the point of given (absolute value of) latitude is in the equatorial region,
/// and `false` if it is located in one of the two polar caps
#[inline]
pub fn is_in_equatorial_region(abs_lat: f64) -> bool {
    abs_lat <= TRANSITION_LATITUDE
}

/// Returns `true` if the point of given (absolute value of) y coordinate in the projected plane
/// is in the equatorial region, and `false` if it is located in one of the two polar caps
#[inline]
pub fn is_in_projected_equatorial_region(abs_y: f64) -> bool {
    abs_y <= 1.0
}

struct AbsAndSign {
    abs: f64,
    sign: u64,
}
#[inline]
fn abs_sign_decompose(x: f64) -> AbsAndSign {
    let bits = f64::to_bits(x);
    AbsAndSign {
        abs: f64::from_bits(bits & F64_BUT_SIGN_BIT_MASK),
        sign: bits & F64_SIGN_BIT_MASK,
    }
}
// Decompose the given positive real value in
// --* an integer offset in [1, 3, 5, 7] (*PI/4) and
// --* a real value in [-1.0, 1.0] (*PI/4)
pub(crate) struct OffsetAndPM1 {
    offset: u8, // = 1, 3, 5 or 7
    pm1: f64,   // in [-1.0, 1.0]
}
#[inline]
pub(crate) fn pm1_offset_decompose(x: f64) -> OffsetAndPM1 {
    let floor: u8 = x as u8;
    let odd_floor: u8 = floor | 1u8;
    OffsetAndPM1 {
        offset: odd_floor & 7u8, // value modulo 8
        pm1: x - (odd_floor as f64),
    }
}

// Cylindrical Equal Area projection
#[inline]
pub(crate) fn proj_cea(xy: &mut (f64, f64)) {
    let (_, ref mut y) = *xy;
    *y = f64::sin(*y) / TRANSITION_Z;
}
#[inline]
fn deproj_cea(lonlat: &mut (f64, f64)) {
    let (_, ref mut lat) = *lonlat;
    // Using asin is OK here since |lat*TRANSITION_Z| < 2/3, so not near from 1.
    *lat = f64::asin((*lat) * TRANSITION_Z);
}

// Collignon projection
#[inline]
pub(crate) fn proj_collignon(xy: &mut (f64, f64)) {
    let (ref mut x, ref mut y) = *xy;
    *y = SQRT_6 * f64::cos(0.5 * *y + FRAC_PI_4);
    *x *= *y;
    *y = 2.0 - *y;
}
#[inline]
fn deproj_collignon(lonlat: &mut (f64, f64)) {
    let (ref mut lon, ref mut lat) = *lonlat;
    *lat = 2.0 - *lat;
    if *lat > 1e-13 {
        // Rare, so few risks of branch miss-prediction
        *lon /= *lat;
        *lon = lon.clamp(-1.0, 1.0);
    } // in case of pole, lon = lat = 0 (we avoid NaN due to division by lat=0)
    *lat /= SQRT_6;
    // Using acos is OK here since lat < 1/sqrt(6), so not near from 1.
    *lat = 2.0 * f64::acos(*lat) - FRAC_PI_2;
}

// Shift x by the given offset and apply lon and lat signs to x and y respectively
#[inline]
fn apply_offset_and_signs(ab: &mut (f64, f64), off: u8, a_sign: u64, b_sign: u64) {
    let (ref mut a, ref mut b) = *ab;
    *a += off as f64;
    *a = f64::from_bits(f64::to_bits(*a) | a_sign);
    *b = f64::from_bits(f64::to_bits(*b) | b_sign);
}
