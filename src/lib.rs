pub mod bmoc;
pub mod checked;
pub mod coords;
pub mod coverage;
pub mod dir;
pub mod geo;
mod hash;
pub mod neighbor;
pub mod proj;
mod special_points_finder;
pub mod sph_geom;
pub mod unchecked;
pub mod zoc;

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
const SQRT_6: f64 = 2.449489743;

pub use coords::LonLat;
use zoc::{ZOC, ZOrderCurve};

use dir::map::CardinalMap;
use dir::set::CardinalSet;
use dir::{Cardinal, Direction};

/// Constant = 29, i.e. the largest possible depth we can store on a signed positive long
/// (4 bits for base cells + 2 bits per depth + 2 remaining bits (1 use in the unique notation).
///
/// ```rust
/// use healpix::{DEPTH_MAX};
/// assert_eq!(29, DEPTH_MAX);
/// ```
pub const MAX_DEPTH: u8 = 29;
/// Constant = nside(29), i.e. the largest possible nside available when we store HEALPix hash
/// on a u64.
///
/// ```rust
/// use healpix::{DEPTH_MAX, NSIDE_MAX, nside};
/// assert_eq!(nside(DEPTH_MAX), NSIDE_MAX);
/// ```
pub const NSIDE_MAX: u32 = 536870912;

/// Limit on the latitude (in radians) between the equatorial region and the polar caps.
/// Equals asin(2/3) = 0.7297276562269663 radians ~= 41,81 degrees.
/// Written $\theta_X$ in Calabretta2007.
///
/// ```rust
/// use healpix::{TRANSITION_LATITUDE};
/// assert_eq!(f64::asin(2f64 / 3f64), TRANSITION_LATITUDE);
/// ```
pub const TRANSITION_LATITUDE: f64 = 0.729_727_656_226_966_3_f64; // asin(2/3)
/// Limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
/// Equals 2/3, see Eq. (1) in Gorsky2005.
pub const TRANSITION_Z: f64 = 2_f64 / 3_f64;
/// Inverse of the limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
/// Equals 1/(2/3) = 1.5, see Eq. (1) in Gorsky2005.
pub const ONE_OVER_TRANSITION_Z: f64 = 1.5_f64;

/// Mask to keep only the f64 sign
const F64_SIGN_BIT_MASK: u64 = 0x8000000000000000;
/// Equals !F64_SIGN_BIT_MASK (the inverse of the f64 sign mask)
const F64_BUT_SIGN_BIT_MASK: u64 = 0x7FFFFFFFFFFFFFFF;

#[inline(always)]
pub const fn is_valid_depth(depth: u8) -> bool {
    depth <= MAX_DEPTH
}
#[inline(always)]
pub const fn assert_valid_depth(depth: u8) {
    debug_assert!(is_valid_depth(depth), "Invalid HEALPix depth");
}

const ONE_OVER_NSIDE: [f64; 30] = [
    1.0,
    0.5,
    0.25,
    0.125,
    0.0625,
    0.03125,
    0.015625,
    0.0078125,
    0.00390625,
    0.001953125,
    0.0009765625,
    0.00048828125,
    0.000244140625,
    0.0001220703125,
    0.00006103515625,
    0.000030517578125,
    0.0000152587890625,
    0.00000762939453125,
    0.000003814697265625,
    0.0000019073486328125,
    0.00000095367431640625,
    0.000000476837158203125,
    0.0000002384185791015625,
    0.00000011920928955078125,
    0.00000005960464477539063,
    0.000000029802322387695313,
    0.000000014901161193847656,
    0.000000007450580596923828,
    0.000000003725290298461914,
    0.000000001862645149230957,
];

const LAYERS: [Layer; 30] = [
    Layer::new(0, ZOC::EMPTY),
    Layer::new(1, ZOC::SMALL),
    Layer::new(2, ZOC::SMALL),
    Layer::new(3, ZOC::SMALL),
    Layer::new(4, ZOC::SMALL),
    Layer::new(5, ZOC::SMALL),
    Layer::new(6, ZOC::SMALL),
    Layer::new(7, ZOC::SMALL),
    Layer::new(8, ZOC::SMALL),
    Layer::new(9, ZOC::MEDIUM),
    Layer::new(10, ZOC::MEDIUM),
    Layer::new(11, ZOC::MEDIUM),
    Layer::new(12, ZOC::MEDIUM),
    Layer::new(13, ZOC::MEDIUM),
    Layer::new(14, ZOC::MEDIUM),
    Layer::new(15, ZOC::MEDIUM),
    Layer::new(16, ZOC::MEDIUM),
    Layer::new(17, ZOC::LARGE),
    Layer::new(18, ZOC::LARGE),
    Layer::new(19, ZOC::LARGE),
    Layer::new(20, ZOC::LARGE),
    Layer::new(21, ZOC::LARGE),
    Layer::new(22, ZOC::LARGE),
    Layer::new(23, ZOC::LARGE),
    Layer::new(24, ZOC::LARGE),
    Layer::new(25, ZOC::LARGE),
    Layer::new(26, ZOC::LARGE),
    Layer::new(27, ZOC::LARGE),
    Layer::new(28, ZOC::LARGE),
    Layer::new(29, ZOC::LARGE),
];

/// Get the [`Layer`] corresponding to the given depth.
pub fn get(depth: u8) -> &'static Layer {
    &LAYERS[depth as usize]
}

pub const fn to_range(hash: u64, delta_depth: u8) -> std::ops::Range<u64> {
    let twice_delta_depth = delta_depth << 1;
    (hash << twice_delta_depth)..((hash + 1) << twice_delta_depth)
}

/// mask ...010101
/// ```rust
/// use healpix::nested::{x_mask};
/// assert_eq!(x_mask(3), 0b00010101);
/// ```
#[inline]
const fn x_mask(depth: u8) -> u64 {
    0x0555555555555555_u64 >> (60 - (depth << 1))
}

/// mask ...101010
/// ```rust
/// use healpix::nested::{y_mask, x_mask};
/// assert_eq!(y_mask(3), 0b00101010);
/// assert_eq!(y_mask(3), x_mask(3) << 1);
/// ```
#[inline]
const fn y_mask(depth: u8) -> u64 {
    0x0AAAAAAAAAAAAAAA_u64 >> (60 - (depth << 1))
}

/// mask ...111111
/// ```rust
/// use healpix::nested::{xy_mask};
/// assert_eq!(xy_mask(3), 0b00111111);
/// let depth = 3_u8;
/// assert_eq!(0xFFFFFFFFFFFFFFFF_u64 >> (64 - (depth << 1)), (1_u64 << (depth << 1)) - 1_u64);
/// ```
#[inline]
const fn xy_mask(depth: u8) -> u64 {
    // 0x0FFFFFFFFFFFFFFF_u64 >> (60 - (depth << 1))
    (1_u64 << (depth << 1)) - 1_u64
}

#[inline(always)]
const fn check_lat(lat: f64) {
    debug_assert!(
        lat >= -FRAC_PI_2 && lat <= FRAC_PI_2,
        "Latitude needs to be on [FRAC_PI_2, FRAC_PI_2]"
    );
}

/// Defines an HEALPix layer in the NESTED scheme.
/// A layer is simply an utility structure containing all constants and methods related
/// to a given depth.
pub struct Layer {
    depth: u8,
    nside: u32,
    nside_minus_1: u32,
    n_hash: u64,
    twice_depth: u8,
    d0h_mask: u64,
    x_mask: u64,
    y_mask: u64,
    xy_mask: u64,
    time_half_nside: i64,
    one_over_nside: f64,
    // z_order_curve: &'static dyn ZOrderCurve,
    z_order_curve: ZOC,
}

impl Layer {
    const fn new(depth: u8, z_order_curve: ZOC) -> Layer {
        // const onef64_to_bits: u64 = 4607182418800017408_u64; // 1_f64.to_bits()
        // let time_nside = (depth as u64) << 52;
        let twice_depth: u8 = depth << 1u8;
        let nside: u32 = unchecked::nside(depth);
        Layer {
            depth,
            nside,
            nside_minus_1: nside - 1,
            n_hash: unchecked::n_hash(depth),
            twice_depth,
            d0h_mask: 15_u64 << twice_depth,
            x_mask: x_mask(depth),
            y_mask: y_mask(depth), //x_mask << 1,
            xy_mask: xy_mask(depth),
            time_half_nside: (depth as i64 - 1) << 52,
            one_over_nside: ONE_OVER_NSIDE[depth as usize], // f64::from_bits(onef64_to_bits - time_nside),
            z_order_curve,
        }
    }

    /// Get the [`Layer`] corresponding to the given depth.
    pub const fn get(depth: u8) -> &'static Layer {
        &LAYERS[depth as usize]
    }

    /// Returns the depth of the Layer (i.e. the HEALPix *order*)
    #[inline]
    pub fn depth(&self) -> u8 {
        self.depth
    }

    /// Returns the number of hash value of the Layer, i.e. the number of cells.
    #[inline]
    pub fn n_hash(&self) -> u64 {
        self.n_hash
    }

    /// Transform the input longitude, in radians, in a value `x` in `[-1, 1[` plus a quarter in `[0, 3]`,
    /// such that `lon = (x + 1) * PI / 4 + q * PI / 2`.
    #[inline]
    fn xpm1_and_q(lon: f64) -> (f64, u8) {
        let lon_bits = lon.to_bits();
        let lon_abs = f64::from_bits(lon_bits & F64_BUT_SIGN_BIT_MASK);
        let lon_sign = lon_bits & F64_SIGN_BIT_MASK;
        let x = lon_abs / FRAC_PI_4;
        let q = x as u8 | 1_u8;
        // Remark: to avoid the branch, we could have copied lon_sign on x - q,
        //         but I so far lack of idea to deal with q efficiently.
        //         And we are not supposed to have negative longitudes in ICRS
        //         (the most used reference system in astronomy).
        if lon_sign == 0 {
            // => lon >= 0
            (x - (q as f64), (q & 7_u8) >> 1)
        } else {
            // case lon < 0 should be rare => few risks of branch miss-prediction
            // Since q in [0, 3]: 3 - (q >> 1)) <=> 3 & !(q >> 1)
            // WARNING: BE SURE TO HANDLE THIS CORRECTLY IN THE REMAINING OF THE CODE!
            //  - Case lon =  3/4 pi = 270 deg => x = -1, q=3
            //  - Case lon = -1/2 pi = -90 deg => x =  1, q=2
            (q as f64 - x, 3 - ((q & 7_u8) >> 1))
        }
    }

    /// Computes the position on the unit sphere of the cell vertex located at the given *direction*
    /// with respect to the center of the cell.
    ///   
    /// # Input
    /// - `hash`: the hash value of the cell we look for the position of a vertex
    /// - `vertex_direction`: the direction of the wanted vertex coordiantes
    ///
    /// # Output
    /// - `(lon, lat)` in radians, the position (on the unit sphere) of the vertex
    ///   - `lon`, longitude in `[0, 2pi]` radians;
    ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    ///
    /// # Example
    /// ```rust
    /// use std::f64::consts::{PI};
    /// use healpix::compass_point::{Cardinal};
    /// use healpix::nested::{get, Layer};
    ///
    /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
    ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
    ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
    /// }
    ///
    /// let depth = 0u8;
    /// let nested0 = get(depth);
    ///
    /// assert!(dist((PI / 4f64, 0.0) , nested0.vertex(0, Cardinal::S)) < 1e-15);
    /// ```
    #[inline]
    pub fn vertex(&self, hash: u64, vertex_direction: Cardinal) -> (f64, f64) {
        let (x, y) = self.center_of_projected_cell(hash);
        self.vertex_lonlat(x, y, vertex_direction)
    }

    /// Computes the positions on the unit sphere of the 4 vertices of the given cell.
    /// If you want to access the position for a given direction, use method
    /// [vertices_map](#method.vertices_map).
    ///   
    /// # Input
    /// - `hash`: the hash value of the cell we look for the positions of its vertices
    ///
    /// # Output
    /// - `[(lon_S, lat_S), (lon_E, lat_E), (lon_N, lat_N), (lon_W, lat_W)]` in radians,
    ///   the positions (on the unit sphere) of the vertices
    ///   - `lon`, longitude in `[0, 2pi]` radians;
    ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    ///
    /// # Example
    /// ```rust
    /// use std::f64::consts::{PI};
    /// use healpix::compass_point::{Cardinal};
    /// use healpix::nested::{get, Layer};
    ///
    /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
    ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
    ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
    /// }
    ///
    /// let depth = 0u8;
    /// let nested0 = get(depth);
    ///
    /// assert!(dist((PI / 4f64, 0.0) , nested0.vertices(0)[0]) < 1e-15);
    /// ```
    #[inline]
    pub fn vertices(&self, hash: u64) -> [(f64, f64); 4] {
        let (x, y) = self.center_of_projected_cell(hash);
        let t = self.one_over_nside;
        [
            proj::unproj(x, y - t),         // S
            proj::unproj(x + t, y),         // E
            proj::unproj(x, y + t),         // N
            proj::unproj((x - t) % 8.0, y), // W
        ]
    }

    /// Computes the projected coordinates of the 4 vertcies of the given cell.
    ///  
    /// # Input
    /// - `hash`: the hash value of the cell we look for the positions of its vertices
    ///
    /// # Output
    /// - `[(x_S, y_S), (x_E, y_E), (x_N, y_N), (x_W, y_W)]` are the projected
    ///   the positions of the vertices.
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    #[inline]
    pub fn projected_vertices(&self, hash: u64) -> [(f64, f64); 4] {
        let (x, y) = self.center_of_projected_cell(hash);
        let t = self.one_over_nside;
        [
            (x, y - t), // S
            (x + t, y), // E
            (x, y + t), // N
            (x - t, y), // W
        ]
    }

    /// Computes the positions on the unit sphere of the vertices of the given cell which direction
    /// are in the given set.
    /// If you don't care about the association between position and direction,
    /// you should use [vertices](#method.vertices).
    ///   
    /// # Input
    /// - `hash`: the hash value of the cell we look for the positions of its vertices
    ///
    /// # Output
    /// - the vertices position stored in a map associating each vertex direction with its position.
    ///   - `lon`, longitude in `[0, 2pi]` radians;
    ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    ///
    /// # Example
    /// ```rust
    /// use std::f64::consts::{PI};
    /// use healpix::compass_point::{Cardinal, CardinalSet};
    /// use healpix::nested::{get, Layer};
    ///
    /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
    ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
    ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
    /// }
    ///
    /// let depth = 0u8;
    /// let nested0 = get(depth);
    ///
    /// assert!(dist((PI / 4f64, 0.0) , *nested0.vertices_map(0, CardinalSet::all()).get(Cardinal::S).unwrap()) < 1e-15);
    /// ```
    #[inline]
    pub fn vertices_map(&self, hash: u64, directions: CardinalSet) -> CardinalMap<(f64, f64)> {
        let (x, y) = self.center_of_projected_cell(hash);
        let mut result_map = CardinalMap::new();
        for direction in directions {
            let vertex = self.vertex_lonlat(x, y, direction);
            result_map.insert(direction, vertex);
        }
        result_map
    }

    /// Computes a list of positions on a given side of a given HEALPix cell on the unit sphere.
    ///
    /// # Input
    /// - `hash`: the hash value of the cell we look for side path on the unit sphere.
    /// - `from_vertex`: direction (from the cell center) of the path starting vertex
    /// - `to_vertex`: direction (from the cell center) of the path ending vertex
    /// - `include_to_vertex`: if set to *false*, the result contains `n_segments` points and do
    ///   not include the ending vertex.
    ///   Else the result contains `n_segments + 1` points.
    /// - `n_segments`: number of segments in the path from the starting vertex to the ending vertex
    ///
    /// # Output
    /// - the list of positions on the given side of the given HEALPix cell on the unit sphere.
    ///
    pub fn path_along_cell_side(
        &self,
        hash: u64,
        from_vertex: Cardinal,
        to_vertex: Cardinal,
        include_to_vertex: bool,
        n_segments: u32,
    ) -> Box<[(f64, f64)]> {
        let n_points: usize = if include_to_vertex {
            n_segments + 1
        } else {
            n_segments
        } as usize;
        let mut path_points: Vec<(f64, f64)> = Vec::with_capacity(n_points);
        let proj_center = self.center_of_projected_cell(hash);
        self.path_along_cell_side_internal(
            proj_center,
            from_vertex,
            to_vertex,
            include_to_vertex,
            n_segments,
            &mut path_points,
        );
        path_points.into_boxed_slice()
    }

    fn path_along_cell_side_internal(
        &self,
        proj_center: (f64, f64),
        from_vertex: Cardinal,
        to_vertex: Cardinal,
        include_to_vertex: bool,
        n_segments: u32,
        path_points: &mut Vec<(f64, f64)>,
    ) {
        let n_points: usize = if include_to_vertex {
            n_segments + 1
        } else {
            n_segments
        } as usize;
        // Compute starting point offsets
        let from_offset_x = from_vertex.x() * self.one_over_nside;
        let from_offset_y = from_vertex.y() * self.one_over_nside;
        // Compute stepX and stepY
        let step_x = (to_vertex.x() * self.one_over_nside - from_offset_x) / (n_segments as f64);
        let step_y = (to_vertex.y() * self.one_over_nside - from_offset_y) / (n_segments as f64);
        // Compute intermediary vertices
        for i in 0..n_points {
            let k = i as f64;
            let x = proj_center.0 + from_offset_x + k * step_x;
            let y = proj_center.1 + from_offset_y + k * step_y;
            path_points.push(proj::unproj(x % 8.0, y));
        }
    }

    /// Computes a list of positions on the edge of a given HEALPix cell on the unit sphere.
    ///
    /// # Input
    /// - `hash`: the hash value of the cell we look for side path on the unit sphere.
    /// - `starting_vertex`: direction (from the cell center) of the path starting vertex
    /// - `clockwise_direction`: tells if the path is in the clockwise or anti-clockwise direction
    /// - `n_segments_by_side`: number of segments in each each side. Hence, the total number of
    ///   points in the path equals *4 x n_segments_by_side*.
    /// # Output
    /// - the list of positions on the given side of the given HEALPix cell on the unit sphere.
    ///
    pub fn path_along_cell_edge(
        &self,
        hash: u64,
        starting_vertex: Cardinal,
        clockwise_direction: bool,
        n_segments_by_side: u32,
    ) -> Box<[(f64, f64)]> {
        // Prepare space for the result
        let mut path_points: Vec<(f64, f64)> =
            Vec::with_capacity((n_segments_by_side << 2) as usize);
        // Compute center
        let proj_center = self.center_of_projected_cell(hash);
        // Unrolled loop over successive sides
        // - compute vertex sequence

        let [v1, v2, v3, v4] = {
            let mut buf = Cardinal::VALUES;
            buf.rotate_left(starting_vertex.index() as _);
            if !clockwise_direction {
                buf.reverse();
            }
            buf
        };
        // - make the four sides
        self.path_along_cell_side_internal(
            proj_center,
            v1,
            v2,
            false,
            n_segments_by_side,
            &mut path_points,
        );
        self.path_along_cell_side_internal(
            proj_center,
            v2,
            v3,
            false,
            n_segments_by_side,
            &mut path_points,
        );
        self.path_along_cell_side_internal(
            proj_center,
            v3,
            v4,
            false,
            n_segments_by_side,
            &mut path_points,
        );
        self.path_along_cell_side_internal(
            proj_center,
            v4,
            v1,
            false,
            n_segments_by_side,
            &mut path_points,
        );
        path_points.into_boxed_slice()
    }

    /// Computes the positions on the sky of each points located on a regular grid in the projection
    /// plane. The grid x-axis is the South-to-east axis and the y-axis is the south-to-west axis.
    /// The return array contains the square fo n_segments_by_side + 1 elements.
    ///
    /// # Input
    /// - `hash`: the hash value of the cell we look for the grid on the unit sphere.
    /// - `n_segments_by_side`: number of segments in each each side. Hence, the total number of
    ///   points in the path equals *(n_segments_by_side + 1)^2*.
    /// # Output
    /// - the list of positions on the given side of the given HEALPix cell on the unit sphere.
    ///
    /// # Motivation
    /// - to create a mesh in Unity
    pub fn grid(&self, hash: u64, n_segments_by_side: u16) -> Box<[(f64, f64)]> {
        let n_points_per_side = (n_segments_by_side as usize) + 1;
        // Prepare space for the result
        let mut grid: Vec<(f64, f64)> = Vec::with_capacity(n_points_per_side * n_points_per_side);
        // Compute center
        let proj_center = self.center_of_projected_cell(hash);
        // Compute grid
        for i in 0..n_points_per_side {
            let x = (i as f64) / (n_segments_by_side as f64); // in [0, 1]
            for j in 0..n_points_per_side {
                let y = (j as f64) / (n_segments_by_side as f64); // in [0, 1]
                let l = x - y;
                let h = x + y - 1.0;
                grid.push(proj::unproj(
                    proj_center.0 + l * self.one_over_nside,
                    proj_center.1 + h * self.one_over_nside,
                ));
            }
        }
        grid.into_boxed_slice()
    }

    /// Center of the given cell in the Euclidean projection space.
    /// # Output
    /// - `(x, y)` coordinates such that $x \in [0, 8[$ and $y \in [-2, 2]$.
    pub fn center_of_projected_cell(&self, hash: u64) -> (f64, f64) {
        self.check_hash(hash);
        let h_parts: HashParts = self.decode_hash(hash);
        let mut hl: (i32, i32) = (
            h_parts.i as i32 - h_parts.j as i32,
            (h_parts.i + h_parts.j) as i32,
        );
        self.shift_from_small_cell_center_to_base_cell_center(&mut hl);
        let (mut x, mut y) = self.scale_to_proj_dividing_by_nside(hl);
        let offset_y = 1 - (h_parts.d0h >> 2) as i8;
        let mut offset_x = (h_parts.d0h & 0b11) << 1u8;
        // +1 if the base cell is not equatorial
        offset_x |= (offset_y & 1_i8) as u8;
        x += offset_x as f64;
        y += offset_y as f64;
        // If x < 0, then x += 8; (happens only in case of base cell 4)
        x += ((x.to_bits() & F64_SIGN_BIT_MASK) >> 60) as f64;
        (x, y)
    }

    // Computes the position on the unit sphere of the vertex, located at the given direction,
    // of the cell of given center coordinate on the projection plane.
    #[inline]
    fn vertex_lonlat(
        &self,
        center_x: f64,
        center_y: f64,
        vertex_direction: Cardinal,
    ) -> (f64, f64) {
        let x = center_x + vertex_direction.x() * self.one_over_nside;
        let y = center_y + vertex_direction.y() * self.one_over_nside;
        proj::unproj(x % 8.0, y)
    }

    #[inline]
    fn is_hash(&self, hash: u64) -> bool {
        hash < self.n_hash
    }

    #[inline]
    fn check_hash(&self, hash: u64) {
        assert!(self.is_hash(hash), "Wrong hash value: too large.");
    }

    #[inline]
    fn h_2_d0h(&self, hash: u64) -> u8 {
        (hash >> self.twice_depth) as u8
    }

    #[inline]
    fn decode_hash(&self, hash: u64) -> HashParts {
        let ij: u64 = self.z_order_curve.h2ij(hash & self.xy_mask);
        HashParts {
            d0h: self.h_2_d0h(hash),
            i: self.z_order_curve.ij2i(ij),
            j: self.z_order_curve.ij2j(ij),
        }
    }

    #[inline]
    fn pull_bits_appart(&self, hash: u64) -> HashBits {
        HashBits {
            d0h: hash & self.d0h_mask,
            i: hash & self.x_mask,
            j: hash & self.y_mask,
        }
    }

    #[inline]
    fn shift_from_small_cell_center_to_base_cell_center(&self, ij: &mut (i32, i32)) {
        let (ref mut _i, ref mut j) = *ij;
        *j -= self.nside_minus_1 as i32;
    }

    #[inline]
    fn scale_to_proj_dividing_by_nside(&self, (x, y): (i32, i32)) -> (f64, f64) {
        (
            x as f64 * self.one_over_nside,
            y as f64 * self.one_over_nside,
        )
    }

    ////////////////////////////
    // Bilinear interpolation //
    ////////////////////////////

    /// See [wikipeida](https://en.wikipedia.org/wiki/Bilinear_interpolation) about bilinear interpolation.
    /// The main difficulty here are the corners of base cells for which the number of neighbors is not
    /// equals to 8.
    /// In the normal case we have:
    /// ```math
    /// f(x, y) = f(0, 0) (1 - x) (1 - y)
    ///         + f(1, 0) x (1 - y)         
    ///         + f(0, 1) (1 - x) y
    ///         + f(1, 1) x y
    /// ```
    /// If a neighbor is missing, we share equally its contribution between the 2 cells that do not
    /// contains the given coordinate, and we fill the array with the cell of the given coordinate
    /// with a weight of 0.  
    /// # Output
    /// - `[(cell, weigth), (cell, weigth), (cell, weigth), (cell, weigth)]` the cell number
    ///   together with their weight
    /// # Panics
    ///   If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
    pub fn bilinear_interpolation(&self, lon: f64, lat: f64) -> [(u64, f64); 4] {
        let (h, dx, dy) = self.hash_with_dxdy(lon, lat);
        // We can probably optimize here since we are interested in only 3 neighbors
        let neigbours_map = self.neighbors(h);
        // Look at the four pixels
        let xcoo = (dx > 0.5) as u8;
        let ycoo = (dy > 0.5) as u8;
        let quarter: u8 = (ycoo << 1) + xcoo;
        match quarter {
            0 => {
                // => S => (dx + 0.5, dy + 0.5, S, SE, SW, C)
                match neigbours_map.get(Direction::S) {
                    Some(nh) => [
                        (*nh, (0.5 - dx) * (0.5 - dy)),
                        (
                            *neigbours_map.get(Direction::SE).unwrap(),
                            (0.5 + dx) * (0.5 - dy),
                        ),
                        (
                            *neigbours_map.get(Direction::SW).unwrap(),
                            (0.5 - dx) * (0.5 + dy),
                        ),
                        (h, (0.5 + dx) * (0.5 + dy)),
                    ],
                    None => [
                        (h, 0.0),
                        (
                            *neigbours_map.get(Direction::SE).unwrap(),
                            (0.5 - dy) * (0.75 + 0.5 * dx),
                        ),
                        (
                            *neigbours_map.get(Direction::SW).unwrap(),
                            (0.5 - dx) * (0.75 + 0.5 * dy),
                        ),
                        (h, (0.5 + dx) * (0.5 + dy)),
                    ],
                }
            }
            1 =>
            // => E => (dx - 0.5, dy + 0.5, SE, E, C, NE)
            {
                match neigbours_map.get(Direction::E) {
                    Some(nh) => [
                        (
                            *neigbours_map.get(Direction::SE).unwrap(),
                            (1.5 - dx) * (0.5 - dy),
                        ),
                        (*nh, (dx - 0.5) * (0.5 - dy)),
                        (h, (1.5 - dx) * (0.5 + dy)),
                        (
                            *neigbours_map.get(Direction::NE).unwrap(),
                            (dx - 0.5) * (0.5 + dy),
                        ),
                    ],
                    None => [
                        (
                            *neigbours_map.get(Direction::SE).unwrap(),
                            (0.5 - dy) * (1.25 - 0.5 * dx),
                        ),
                        (h, 0.0),
                        (h, (1.5 - dx) * (0.5 + dy)),
                        (
                            *neigbours_map.get(Direction::NE).unwrap(),
                            (dx - 0.5) * (0.75 + 0.5 * dy),
                        ),
                    ],
                }
            }
            2 =>
            // => W => (dx + 0.5, dy - 0.5, SW, C, W, NW)
            {
                match neigbours_map.get(Direction::W) {
                    Some(nh) => [
                        (
                            *neigbours_map.get(Direction::SW).unwrap(),
                            (0.5 - dx) * (1.5 - dy),
                        ),
                        (h, (dx + 0.5) * (1.5 - dy)),
                        (*nh, (0.5 - dx) * (dy - 0.5)),
                        (
                            *neigbours_map.get(Direction::NW).unwrap(),
                            (0.5 + dx) * (dy - 0.5),
                        ),
                    ],
                    None => [
                        (
                            *neigbours_map.get(Direction::SW).unwrap(),
                            (0.5 - dx) * (1.25 - 0.5 * dy),
                        ),
                        (h, (dx + 0.5) * (1.5 - dy)),
                        (h, 0.0),
                        (
                            *neigbours_map.get(Direction::NW).unwrap(),
                            (dy - 0.5) * (0.5 * dx + 0.75),
                        ),
                    ],
                }
            }
            3 =>
            // => N => (dx - 0.5, dy - 0.5, C, NE, NW, N)
            {
                match neigbours_map.get(Direction::N) {
                    Some(nh) => [
                        (h, (1.5 - dx) * (1.5 - dy)),
                        (
                            *neigbours_map.get(Direction::NE).unwrap(),
                            (dx - 0.5) * (1.5 - dy),
                        ),
                        (
                            *neigbours_map.get(Direction::NW).unwrap(),
                            (1.5 - dx) * (dy - 0.5),
                        ),
                        (*nh, (dx - 0.5) * (dy - 0.5)),
                    ],
                    None => [
                        (h, (1.5 - dx) * (1.5 - dy)),
                        (
                            *neigbours_map.get(Direction::NE).unwrap(),
                            (dx - 0.5) * (1.25 - 0.5 * dy),
                        ),
                        (
                            *neigbours_map.get(Direction::NW).unwrap(),
                            (1.25 - 0.5 * dx) * (dy - 0.5),
                        ),
                        (h, 0.0),
                    ],
                }
            }
            _ => unreachable!(),
        }
    }
}

struct HashParts {
    d0h: u8, // base cell number (depth 0 hash value)
    i: u32,  // in the base cell, z-order curve coordinate along the x-axis
    j: u32,  // in the base cell, z-order curve coordinate along the x-axis
}

struct HashBits {
    d0h: u64, // base cell number (depth 0 hash value) bits
    i: u64,   // in the base cell, z-order curve coordinate along the x-axis bits
    j: u64,   // in the base cell, z-order curve coordinate along the y-axis bits
}
