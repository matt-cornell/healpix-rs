use crate::{FRAC_PI_2, FRAC_PI_4, SQRT_6, TRANSITION_LATITUDE, TRANSITION_Z, proj, unchecked};

impl super::Layer {
    /// Returns the cell number (hash value) associated with the given position on the unit sphere
    /// # Inputs
    /// - `lon`: longitude in radians, support reasonably large positive and negative values
    ///   producing accurate results with a naive range reduction like modulo 2*pi
    ///   (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
    /// - `lat`: latitude in radians, must be in `[-pi/2, pi/2]`
    /// # Output
    /// - the cell number (hash value) associated with the given position on the unit sphere,
    ///   in `[0, 12*nside^2[`
    /// # Panics
    ///   If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
    /// # Examples
    /// ```rust
    /// use healpix::{nside};
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 12_u8;
    /// let nside = nside(depth) as u64;
    /// let nested12: &Layer = get(depth);
    /// assert_eq!(nside * nside - 1, nested12.hash(12.5_f64.to_radians(), 89.99999_f64.to_radians()));
    /// ```
    pub fn hash(&self, lon: f64, lat: f64) -> u64 {
        crate::check_lat(lat);
        let (d0h, l_in_d0c, h_in_d0c) = Self::d0h_lh_in_d0c(lon, lat);
        // Coords inside the base cell
        //  - ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
        let i =
            f64::from_bits((self.time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64)
                as u32;
        let j =
            f64::from_bits((self.time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64)
                as u32;
        //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
        let i = if i == self.nside {
            self.nside_minus_1
        } else {
            i
        };
        let j = if j == self.nside {
            self.nside_minus_1
        } else {
            j
        };
        self.build_hash_from_parts(d0h, i, j)
    }

    pub fn try_hash(&self, lon: f64, lat: f64) -> Result<u64, f64> {
        if lat.abs() > FRAC_PI_2 {
            return Err(lat);
        }
        Ok(self.hash(lon, lat))
    }

    #[inline]
    fn d0h_lh_in_d0c(lon: f64, lat: f64) -> (u8, f64, f64) {
        let (x_pm1, q) = Self::xpm1_and_q(lon);
        if lat > TRANSITION_LATITUDE {
            // North polar cap, Collignon projection.
            // - set the origin to (PI/4, 0)
            let sqrt_3_one_min_z = SQRT_6 * (0.5 * lat + FRAC_PI_4).cos();
            let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, 2.0 - sqrt_3_one_min_z);
            let d0h = q;
            (d0h, x_proj, y_proj)
        } else if lat < -TRANSITION_LATITUDE {
            // South polar cap, Collignon projection
            // - set the origin to (PI/4, -PI/2)
            let sqrt_3_one_min_z = SQRT_6 * (0.5 * lat - FRAC_PI_4).cos(); // cos(-x) = cos(x)
            let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, sqrt_3_one_min_z);
            let d0h = q + 8;
            (d0h, x_proj, y_proj)
        } else {
            // Equatorial region, Cylindrical equal area projection
            // - set the origin to (PI/4, 0)               if q = 2
            // - set the origin to (PI/4, -PI/2)           if q = 0
            // - set the origin to (0, -TRANSITION_LAT)    if q = 3
            // - set the origin to (PI/2, -TRANSITION_LAT) if q = 1
            // let zero_or_one = (x_cea as u8) & 1;
            let y_pm1 = lat.sin() / TRANSITION_Z;
            // Inequalities have been carefully chosen so that S->E and S->W axis are part of the cell,
            // and not E->N and W->N

            // Version with branch
            // |\3/|
            // .2X1.
            // |/0\|
            /*let q13 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */  debug_assert!(q12 == 0 || q12 == 1);
            let q23 = (x_pm1 <=  y_pm1) as u8; /* 1/0 */  debug_assert!(q23 == 0 || q23 == 1);
            match q13 | (q23 << 1) {
              0 => ( q         , x_pm1      , y_pm1 + 2.0),
              1 => ((q + 5) & 7, x_pm1 - 1.0, y_pm1 + 1.0), // (q + 5) & 7 <=> (q + 1) | 4
              2 => ( q + 4     , x_pm1 + 1.0, y_pm1 + 1.0),
              3 => ( q + 8     , x_pm1      , y_pm1),
              _ => unreachable!(),
            }*/
            // Branch free version
            // |\2/|
            // .3X1.
            // |/0\|
            let q01 = (x_pm1 > y_pm1) as u8; /* 0/1 */
            debug_assert!(q01 == 0 || q01 == 1);
            let q12 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */
            debug_assert!(q12 == 0 || q12 == 1);
            let q1 = q01 & q12; /* = 1 if q1, 0 else */
            debug_assert!(q1 == 0 || q1 == 1);
            let q013 = q01 + (1 - q12); // = q01 + q03; /* 0/1 + 1\0 +  */
            // x: x_pm1 + 1 if q3 | x_pm1 - 1 if q1 | x_pm1 if q0 or q2
            let x_proj = x_pm1 - ((q01 + q12) as i8 - 1) as f64;
            // y: y_pm1 + 0 if q2 | y_pm1 + 1 if q1 or q3 | y_pm1 + 2 if q0
            let y_proj = y_pm1 + q013 as f64;
            // d0h: +8 if q0 | +4 if q3 | +5 if q1
            let d0h = (q013 << 2) + ((q + q1) & 3);
            (d0h, x_proj, y_proj)
        }
    }

    /// Returns the cell number (hash value) associated with the given position on the unit sphere,
    /// together with the offset `(dx, dy)` on the Euclidean plane of the projected position with
    /// respect to the origin of the cell (South vertex).
    /// # Inputs
    /// - `lon`: longitude in radians, support reasonably large positive and negative values
    ///   producing accurate results with a naive range reduction like modulo 2*pi
    ///   (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
    /// - `lat`: latitude in radians, must be in `[-pi/2, pi/2]`
    /// # Output
    /// - the cell number (hash value) associated with the given position on the unit sphere,
    ///   in `[0, 12*nside^2[`
    /// - `dx`: the positional offset $\in [0, 1[$ along the south-to-east axis
    /// - `dy`: the positional offset $\in [0, 1[$ along the south-to-west axis
    /// # Panics
    ///   If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
    /// # Examples
    /// ```rust
    /// use healpix::{nside};
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 12_u8;
    /// let nside = nside(depth) as u64;
    /// let nested12: &Layer = get(depth);
    /// let h_org = nside * nside - 1;
    /// let (h_ra, h_dec) = nested12.center(h_org);
    /// let (h, dx, dy) = nested12.hash_with_dxdy(h_ra, h_dec);
    /// assert_eq!(h_org, h);
    /// // A precision of 1e-12 in a cell of depth 12 (side of ~51.5 arcsec)
    /// // leads to an absolute precision of ~0.05 nanoarcsec
    /// assert!((dx - 0.5).abs() < 1e-12_f64);
    /// assert!((dy - 0.5).abs() < 1e-12_f64);
    /// ```
    pub fn hash_with_dxdy(&self, lon: f64, lat: f64) -> (u64, f64, f64) {
        crate::check_lat(lat);
        let (d0h, l_in_d0c, h_in_d0c) = Self::d0h_lh_in_d0c(lon, lat);
        let x =
            f64::from_bits((self.time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64);
        debug_assert!(
            -1e-14 < x || x < self.nside as f64 * (1.0_f64 + 1e-14),
            "x: {}, x_proj: {}; y_proj: {}",
            &x,
            &h_in_d0c,
            &l_in_d0c
        );
        let y =
            f64::from_bits((self.time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64);
        debug_assert!(
            -1e-14 < y || y < self.nside as f64 * (1.0_f64 + 1e-14),
            "y: {}, x_proj: {}; y_proj: {}",
            &y,
            &h_in_d0c,
            &l_in_d0c
        );
        // - ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
        let i = x as u32;
        let j = y as u32;
        //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
        let i = if i == self.nside {
            self.nside_minus_1
        } else {
            i
        };
        let j = if j == self.nside {
            self.nside_minus_1
        } else {
            j
        };
        (
            self.build_hash_from_parts(d0h, i, j),
            (x - (i as f64)),
            (y - (j as f64)),
        )
    }

    /// Compute the position on the unit sphere of the center (in the Euclidean projection plane)
    /// of the cell associated to the given hash value.
    ///
    /// # Input
    /// - `hash`: the hash value of the cell we look for the unprojected center
    ///
    /// # Output
    /// - `(lon, lat)` in radians, the unprojected position (on the unit sphere) of the center of
    ///   the cell in the Euclidean plane
    ///   - `lon`, longitude in `[0, 2pi]` radians;
    ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    ///
    /// # Example
    /// ```rust
    /// use std::f64::consts::{PI};
    /// use healpix::{TRANSITION_LATITUDE};
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
    /// assert!(dist((PI / 4f64, TRANSITION_LATITUDE) , nested0.center(0u64)) < 1e-15);
    /// ```
    ///
    #[inline]
    pub fn center(&self, hash: u64) -> (f64, f64) {
        let (x, y) = self.center_of_projected_cell(hash);
        proj::unproj(x, y)
    }

    /// Compute the position on the unit sphere of the position '(dx, dy)' from the south vertex of
    /// the HEALPix cell associated to the given hash value.
    /// The x-axis is the South-East axis while the y-axis is the south-west axis.
    ///
    /// # Input
    /// - `hash`: the hash value of the cell in which are defined `dx` and `dy`
    /// - `dx`: the positional offset $\in [0, 1[$ along the south-to-east axis
    /// - `dy`: the positional offset $\in [0, 1[$ along the south-to-west axis
    ///
    /// # Output
    /// - `(lon, lat)` in radians, the unprojected position (on the unit sphere) of the given position
    ///   inside the given cell in the Euclidean plane
    ///   - `lon`, longitude in `[0, 2pi]` radians;
    ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
    ///
    /// # Panics
    /// This method panics if either:
    /// - the given `hash` value is not in `[0, 12*nside^2[`,
    /// - `dx` or `dy` is not $\in [0, 1[$
    ///
    /// # Example
    /// ```rust
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
    /// assert!(dist(nested0.sph_coo(0, 0.5, 0.5) , nested0.center(0)) < 1e-15);
    /// ```
    ///
    pub fn sph_coo(&self, hash: u64, dx: f64, dy: f64) -> (f64, f64) {
        assert!((0.0..1.0).contains(&dx));
        assert!((0.0..1.0).contains(&dy));
        let (mut x, mut y) = self.center_of_projected_cell(hash);
        x += (dx - dy) * self.one_over_nside;
        y += (dx + dy - 1.0) * self.one_over_nside;
        proj::unproj(x.rem_euclid(8.0), y)
    }

    #[inline]
    pub(crate) fn build_hash_from_parts(&self, d0h: u8, i: u32, j: u32) -> u64 {
        self.build_hash((d0h as u64) << self.twice_depth, i, j)
    }

    #[inline]
    pub(crate) fn build_hash(&self, d0h_bits: u64, i: u32, j: u32) -> u64 {
        use crate::zoc::ZOrderCurve;
        debug_assert!(
            i < self.nside && j < self.nside,
            "nside: {}; i: {}, j: {}",
            self.nside,
            i,
            j
        );
        d0h_bits | self.z_order_curve.ij2h(i, j)
    }

    /// Conveniency function simply calling the [hash](fn.to_uniq.html) method with this layer *depth*.
    pub fn to_uniq(&self, hash: u64) -> u64 {
        // depth already tested, so we call the unsafe method
        unchecked::to_uniq(self.depth, hash)
    }

    /// Conveniency function simply calling the [hash](fn.to_uniq_ivoa.html) method with this layer *depth*.
    pub fn to_uniq_ivoa(&self, hash: u64) -> u64 {
        // depth already tested, so we call the unsafe method
        unchecked::to_uniq_ivoa(self.depth, hash)
    }
}
