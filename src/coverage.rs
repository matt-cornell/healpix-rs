impl super::Layer {
    /// Returns a hierarchical view of the list of cells overlapped by the given zone.
    ///
    /// # Input
    /// - `lon_min` the longitude of the bottom left corner
    /// - `lat_min` the latitude of the bottom left corner
    /// - `lon_max` the longitude of the upper left corner
    /// - `lat_max` the latitude of the upper left corner
    ///
    /// # Output
    /// - the list of cells overlapped by the given zone, in a BMOC (hierarchical view also telling
    ///   if a cell is fully or partially covered).
    ///
    /// # Remark
    /// - If `lon_min > lon_max` then we consider that the zone crosses the primary meridian.
    /// - The north pole is included only if `lon_min == 0 && lat_max == pi/2`
    ///
    /// # Panics
    /// * if `lon_min` or `lon_max` not in `[0, 2\pi[`
    /// * if `lat_min` or `lat_max` not in `[-\pi/2, \pi/2[`
    /// * `lat_min >= lat_max`.
    pub fn zone_coverage(&self, lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64) -> BMOC {
        let zone = Zone::new(lon_min, lat_min, lon_max, lat_max);
        let (depth_start, hashs_start) = zone
            .smallest_enclosing_cone()
            .map(|mec| {
                let center = mec.center().lonlat();
                let wider_than_180deg = zone.dlon() > PI;
                if !wider_than_180deg
                    && zone.contains(center.lon(), center.lat())
                    && has_best_starting_depth(mec.radius())
                {
                    let d = best_starting_depth(mec.radius()).min(self.depth);
                    let h = hash(d, center.lon(), center.lat());
                    (d, neighbors(d, h, true).sorted_values_vec())
                } else {
                    (0, (0..12).collect::<Vec<u64>>())
                }
            })
            .unwrap_or((0, (0..12).collect::<Vec<u64>>()));
        let zone_perimeter = 2.0 * (zone.dlon() + zone.dlat());
        let approx_cell_side_at_depth_max = 2.0 * (PI / n_hash(self.depth) as f64).sqrt();
        let mut bmoc_builder = BMOCBuilderUnsafe::new(
            self.depth,
            6 * (2 + (zone_perimeter / approx_cell_side_at_depth_max) as usize),
        );

        // What about considering the zone in the projection plane?
        // With 3 parts, in all 3 cases an HEALPix cell is a diamond
        // * north polar cap: zone = a trapezoid
        // * equatorial region: zone = a rectangle
        // * south polar cap: zone = a trapezoid

        let zone_vertices = zone.vertices();
        // The flag tells that the zone vertex is NOT on the South HEALPix vertex (and thus is not to be
        // ignore at the deepest level).
        let mut zone_vertices_hashs_flags = Vec::with_capacity(4);
        if zone_vertices[0].1 > -HALF_PI {
            let vertex_hash_sw = self.hash(zone_vertices[0].0, zone_vertices[0].1);
            zone_vertices_hashs_flags.push((vertex_hash_sw, true));
        }
        // Others are not included in the zone
        if zone_vertices[1].1 < HALF_PI {
            let (vertex_hash_nw, dx, dy) =
                self.hash_with_dxdy(zone_vertices[1].0, zone_vertices[1].1);
            zone_vertices_hashs_flags.push((
                vertex_hash_nw,
                dx > 1e-15 && dy > 1e-15 && dx < 1.0 && dy < 1.0,
            ));
        }
        if zone_vertices[2].1 < HALF_PI {
            let (vertex_hash_ne, dx, dy) =
                self.hash_with_dxdy(zone_vertices[2].0, zone_vertices[2].1);
            zone_vertices_hashs_flags.push((
                vertex_hash_ne,
                dx > 1e-15 && dy > 1e-15 && dx < 1.0 && dy < 1.0,
            ));
        }
        if zone_vertices[3].1 > -HALF_PI {
            let (vertex_hash_se, dx, dy) =
                self.hash_with_dxdy(zone_vertices[3].0, zone_vertices[3].1);
            zone_vertices_hashs_flags.push((
                vertex_hash_se,
                dx > 1e-15 && dy > 1e-15 && dx < 1.0 && dy < 1.0,
            ));
        }

        for h in hashs_start {
            self.zone_coverage_recur(
                depth_start,
                h,
                &zone,
                &zone_vertices_hashs_flags,
                &mut bmoc_builder,
            );
        }
        bmoc_builder.to_bmoc()
    }

    fn zone_coverage_recur(
        &self,
        depth: u8,
        hash: u64,
        zone: &Zone,
        zone_vertices_hashs_flags: &Vec<(u64, bool)>,
        bmoc_builder: &mut BMOCBuilderUnsafe,
    ) {
        let shift = (self.depth - depth) << 1; // <=> shift = 2 * (self.depth - depth)
        let matches_zone_vertex_hash = |hash: u64| {
            if shift == 0 {
                // At the deeper depth, we do not consider zone vertices on the edges of an HEALPix cell,
                // except for the SW zone vertex (which is included in the zone).
                for (h, flag) in zone_vertices_hashs_flags {
                    if *flag && *h == hash {
                        return true;
                    }
                }
            } else {
                for (h, _) in zone_vertices_hashs_flags {
                    if *h >> shift == hash {
                        return true;
                    }
                }
            }
            false
        };
        let [(l_s, b_s), (l_e, b_e), (l_n, b_n), (l_w, b_w)] = vertices(depth, hash);
        let n_in = zone.contains(l_s, b_s) as u8
            + zone.contains_exclusive(l_e, b_e) as u8
            + zone.contains_exclusive(l_n, b_n) as u8
            + zone.contains_exclusive(l_w, b_w) as u8;
        // eprintln!("depth: {}; hash: {}; n_in: {}", depth, hash, n_in);

        // A cell may intersect a zone without having a vertex inside:
        // * either a vertex of the zone is in the cell: we can easily test this
        // * or no vertex of the zone is int the cell: we must detect either
        //     + the cell NS  vertical line intersect the zone => (N_lat > lat_max && S_lat < lat_min) && NS_lon_NS in dlon
        //     + the cell EW horizontal line intersect the zone => (EW_lat in dlat) && E_lon in ...
        // Remark: in the projected plane:
        // * Equatorial region: the zone is a rectangle
        // * Polar caps: the zone is a trapezoid
        // * A cell is always a diamond

        if n_in == 4 && zone.is_lon_range_compatible(l_w, l_e) {
            bmoc_builder.push(depth, hash, true);
        } else if n_in > 0 || matches_zone_vertex_hash(hash)
        // we may have false positive in the polar caps around n*PI/2
        || zone.crossed_vertically(l_n, b_s, b_n)
        || zone.crossed_vertically(l_s, b_s, b_n)
        || zone.crossed_horizontally(l_w, l_e, b_w)
        || zone.crossed_horizontally(l_w, l_e, b_e)
        {
            if depth == self.depth {
                bmoc_builder.push(depth, hash, false);
            } else {
                let hash = hash << 2;
                let depth = depth + 1;
                self.zone_coverage_recur(
                    depth,
                    hash,
                    zone,
                    zone_vertices_hashs_flags,
                    bmoc_builder,
                );
                self.zone_coverage_recur(
                    depth,
                    hash | 1_u64,
                    zone,
                    zone_vertices_hashs_flags,
                    bmoc_builder,
                );
                self.zone_coverage_recur(
                    depth,
                    hash | 2_u64,
                    zone,
                    zone_vertices_hashs_flags,
                    bmoc_builder,
                );
                self.zone_coverage_recur(
                    depth,
                    hash | 3_u64,
                    zone,
                    zone_vertices_hashs_flags,
                    bmoc_builder,
                );
            }
        }
    }

    /// Returns a hierarchical view of the list of cells having their center in the given cone.
    /// All BMOC flags are set to 1.
    ///
    /// # Input
    /// - `cone_lon` the longitude of the center of the cone, in radians
    /// - `cone_lat` the latitude of the center of the cone, in radians
    /// - `cone_radius` the radius of the cone, in radians
    ///
    /// # Output
    /// - the list of cells having their center in the given cone, in a BMOC (hierarchical view with
    ///   all flags set to 1).
    ///
    /// # Example
    /// ```rust
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 4_u8;
    /// let nested = get(depth);
    ///
    /// let lon = 13.158329_f64.to_radians();
    /// let lat = -72.80028_f64.to_radians();
    /// let radius = 5.64323_f64.to_radians();
    ///
    /// let actual_res = nested.cone_coverage_centers(lon, lat, radius);
    /// let expected_res: [u64; 7] = [2058, 2059, 2080, 2081, 2082, 2083, 2088];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///      assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn cone_coverage_centers(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
        // Special case: the full sky is covered
        if cone_radius >= PI {
            return self.allsky_bmoc_builder().to_bmoc();
        }
        // Common variable
        let cos_cone_lat = cone_lat.cos();
        let shs_cone_radius = to_squared_half_segment(cone_radius);
        // Special case of very large radius: test the 12 base cells
        if !has_best_starting_depth(cone_radius) {
            let distances = largest_center_to_vertex_distances_with_radius(
                0,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius,
            );
            let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius),
            );
            for h in 0..12 {
                self.cone_coverage_centers_recur(
                    0,
                    h,
                    shs_cone_radius,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    &minmax_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            return bmoc_builder.to_bmoc_packing();
        }
        // Normal case
        let depth_start = best_starting_depth(cone_radius);
        if depth_start >= self.depth {
            let neigs: Vec<u64> = self
                .neighbors(self.hash(cone_lon, cone_lat), true)
                .values_vec()
                .into_iter()
                .filter(|neigh| {
                    let (lon, lat) = self.center(*neigh);
                    squared_half_segment(lon - cone_lon, lat - cone_lat, lon.cos(), cos_cone_lat)
                        <= shs_cone_radius
                })
                .collect();
            let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
            for neig in neigs {
                bmoc_builder.push(self.depth, neig, true);
            }
            bmoc_builder.to_bmoc()
        } else {
            let distances = largest_center_to_vertex_distances_with_radius(
                depth_start,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius,
            );
            let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
            let root_layer = get(depth_start);
            let root_center_hash = root_layer.hash(cone_lon, cone_lat);
            let neigs = root_layer.neighbors(root_center_hash, true);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius),
            );
            for &root_hash in neigs.sorted_values().iter() {
                self.cone_coverage_centers_recur(
                    depth_start,
                    root_hash,
                    shs_cone_radius,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    &minmax_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            bmoc_builder.to_bmoc_packing()
        }
    }

    // TODO: make a generic function with cone_coverage_approx_recur or at least cone_coverage_fullin_recur
    #[allow(clippy::too_many_arguments)]
    fn cone_coverage_centers_recur<F>(
        &self,
        depth: u8,                            // cell depth
        hash: u64,                            // cell hash
        shs_cone_radius: f64, // squared_half_segment corresponding to the cone radius
        shs_computer: &F, // function to compute the squared_half_segment between a point and the cone center
        shs_minmax: &[MinMax], // min and max shs at various delta depth
        recur_depth: u8,  // conveniency delta depth index
        bmoc_builder: &mut BMOCBuilderUnsafe, // builder in which cells are appended (flag always set to true!)
    ) where
        F: Fn((f64, f64)) -> f64,
    {
        let center = get(depth).center(hash);
        let shs = shs_computer(center);
        let MinMax { min, max } = shs_minmax[recur_depth as usize];
        if shs <= min {
            bmoc_builder.push(depth, hash, true);
        } else if shs <= max {
            if depth == self.depth {
                if shs <= shs_cone_radius {
                    bmoc_builder.push(depth, hash, true);
                } // else do nothing, recur end here
            } else {
                // because 'min' is a lower bound (approx), we may split cells fully included in the MOC
                // we could check the 4 corers here, we assume it is more costly than post-merging (to be verified!)
                let hash = hash << 2;
                let depth = depth + 1;
                let recur_depth = recur_depth + 1;
                self.cone_coverage_centers_recur(
                    depth,
                    hash,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_centers_recur(
                    depth,
                    hash | 1_u64,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_centers_recur(
                    depth,
                    hash | 2_u64,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_centers_recur(
                    depth,
                    hash | 3_u64,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
            }
        }
    }

    /// Returns a hierarchical view of the list of cells fully covered by the given cone.
    /// All BMOC flags are set to 1.
    ///
    /// # Input
    /// - `cone_lon` the longitude of the center of the cone, in radians
    /// - `cone_lat` the latitude of the center of the cone, in radians
    /// - `cone_radius` the radius of the cone, in radians
    ///
    /// # Output
    /// - the list of cells fully covered by the given cone, in a BMOC (hierarchical view with
    ///   all flags set to 1).
    ///
    /// # Example
    /// ```rust
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 4_u8;
    /// let nested3 = get(depth);
    ///
    /// let lon = 13.158329_f64.to_radians();
    /// let lat = -72.80028_f64.to_radians();
    /// let radius = 5.64323_f64.to_radians();
    ///
    /// let actual_res = nested3.cone_coverage_fullin(lon, lat, radius);
    /// let expected_res: [u64; 2] = [2081, 2082];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///      assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn cone_coverage_fullin(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
        // Special case: the full sky is covered
        if cone_radius >= PI {
            return self.allsky_bmoc_builder().to_bmoc();
        }
        // Common variable
        let cos_cone_lat = cone_lat.cos();
        let shs_cone_radius = to_squared_half_segment(cone_radius);
        // Special case of very large radius: test the 12 base cells
        if !has_best_starting_depth(cone_radius) {
            let distances = largest_center_to_vertex_distances_with_radius(
                0,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius,
            );
            let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius),
            );
            for h in 0..12 {
                self.cone_coverage_fullin_recur(
                    0,
                    h,
                    shs_cone_radius,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    &minmax_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            return bmoc_builder.to_bmoc_packing();
        }
        // Normal case
        let depth_start = best_starting_depth(cone_radius);
        if depth_start >= self.depth {
            let neigs: Vec<u64> = self
                .neighbors(self.hash(cone_lon, cone_lat), true)
                .values_vec()
                .into_iter()
                .filter(|neigh| {
                    let mut fully_in = true;
                    for (lon, lat) in self.vertices(*neigh) {
                        fully_in &= squared_half_segment(
                            lon - cone_lon,
                            lat - cone_lat,
                            lon.cos(),
                            cos_cone_lat,
                        ) <= shs_cone_radius;
                    }
                    fully_in
                })
                .collect();
            let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
            for neig in neigs {
                bmoc_builder.push(self.depth, neig, true);
            }
            bmoc_builder.to_bmoc()
        } else {
            let distances = largest_center_to_vertex_distances_with_radius(
                depth_start,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius,
            );
            let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
            let root_layer = get(depth_start);
            let root_center_hash = root_layer.hash(cone_lon, cone_lat);
            let neigs = root_layer.neighbors(root_center_hash, true);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius),
            );
            for &root_hash in neigs.sorted_values().iter() {
                self.cone_coverage_fullin_recur(
                    depth_start,
                    root_hash,
                    shs_cone_radius,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    &minmax_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            bmoc_builder.to_bmoc_packing()
        }
    }

    // TODO: make a generic function with cone_coverage_approx_recur or at least cone_coverage_centers_recur
    #[allow(clippy::too_many_arguments)]
    fn cone_coverage_fullin_recur<F>(
        &self,
        depth: u8,                            // cell depth
        hash: u64,                            // cell hash
        shs_cone_radius: f64, // squared_half_segment corresponding to the cone radius
        shs_computer: &F, // function to compute the squared_half_segment between a point and the cone center
        shs_minmax: &[MinMax], // min and max shs at various delta depth
        recur_depth: u8,  // conveniency delta depth index
        bmoc_builder: &mut BMOCBuilderUnsafe, // builder in which cells are appended (flag always set to true!)
    ) where
        F: Fn((f64, f64)) -> f64,
    {
        let center = get(depth).center(hash);
        let shs = shs_computer(center);
        let MinMax { min, max } = shs_minmax[recur_depth as usize];
        if shs <= min {
            bmoc_builder.push(depth, hash, true);
        } else if shs <= max {
            if depth == self.depth {
                // Center must be in the cone
                if shs <= shs_cone_radius {
                    // And the 4 vertices must also be in the cone
                    for vertex_coo in self.vertices(hash) {
                        // If vertex out of the cone
                        if shs_computer(vertex_coo) > shs_cone_radius {
                            return; // do nothing, recur end here
                        }
                    }
                    bmoc_builder.push(depth, hash, true);
                } // else do nothing, recur end here
            } else {
                // because 'min' is a lower bound (approx), we may split cells fully included in the MOC
                // we could check the 4 corers here, we assume it is more costly than post-merging (to be verified!)
                let hash = hash << 2;
                let depth = depth + 1;
                let recur_depth = recur_depth + 1;
                self.cone_coverage_fullin_recur(
                    depth,
                    hash,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_fullin_recur(
                    depth,
                    hash | 1_u64,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_fullin_recur(
                    depth,
                    hash | 2_u64,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_fullin_recur(
                    depth,
                    hash | 3_u64,
                    shs_cone_radius,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
            }
        }
    }

    /// Returns a hierarchical view of the list of cells fully covered by the given cone.
    /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
    /// This is approximated: cell tagged as partially overlapped could be fully overlapped
    /// (but we guarantee that cells tagged as fully overlapped are really fully overlapped).
    ///
    /// # Input
    /// - `cone_lon` longitude of the center of the ring, in radians
    /// - `cone_lat` latitude of the center of the ring, in radians
    /// - `cone_radius_int` internal radius of the ring, in radians
    /// - `cone_radius_ext` external radius of the ring, in radians
    ///
    /// # Output
    /// - the list of cells overlapped by the given ring, in a BMOC (hierarchical view also telling
    ///   if a cell is fully or partially covered).
    ///
    /// # Panics
    /// * if `cone_radius_ext > PI`
    /// * if `cone_radius_ext <= cone_radius_int`
    ///
    /// # Usecase
    /// * annulus region provided by the delay between Fermi and INTEGRAL for GRBs
    ///   (asked for Gravitational Waves applications).
    ///
    ///
    /// # Example
    /// ```rust
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 4_u8;
    /// let nested = get(depth);
    ///
    /// let lon = 13.158329_f64.to_radians();
    /// let lat = -72.80028_f64.to_radians();
    /// let radius_int = 5.64323_f64.to_radians();
    /// let radius_ext = 10.0_f64.to_radians();
    ///
    /// let actual_res = nested.ring_coverage_approx(lon, lat, radius_int, radius_ext);
    /// let expected_res: [u64; 40] = [2050, 2051, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061,
    ///   2062, 2063, 2080, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2094, 2176,
    ///   2177, 2178, 2817, 2820, 2821, 2822, 2823, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2880];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///      assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn ring_coverage_approx(
        &self,
        cone_lon: f64,
        cone_lat: f64,
        cone_radius_int: f64,
        cone_radius_ext: f64,
    ) -> BMOC {
        self.ring_coverage_approx_internal(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
            .to_bmoc_packing()
    }

    pub fn ring_coverage_approx_custom(
        &self,
        delta_depth: u8,
        cone_lon: f64,
        cone_lat: f64,
        cone_radius_int: f64,
        cone_radius_ext: f64,
    ) -> BMOC {
        if delta_depth == 0 {
            self.ring_coverage_approx(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
        } else {
            // TODO: change the algo not to put all cell in the MOC before pruning it:
            // - make a second recur function returning the number of sub-cell overlapped by the cone
            get(self.depth + delta_depth)
                .ring_coverage_approx_internal(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
                .to_lower_depth_bmoc_packing(self.depth)
        }
    }

    fn ring_coverage_approx_internal(
        &self,
        cone_lon: f64,
        cone_lat: f64,
        cone_radius_int: f64,
        cone_radius_ext: f64,
    ) -> BMOCBuilderUnsafe {
        assert!(cone_radius_ext <= PI);
        assert!(
            cone_radius_int < cone_radius_ext,
            "{} >= {} ",
            cone_radius_int,
            cone_radius_ext
        );
        // Common variable
        let cos_cone_lat = cone_lat.cos();
        let shs_cone_radius_int = to_squared_half_segment(cone_radius_int);
        // Special case of very large radius: test the 12 base cells
        if !has_best_starting_depth(cone_radius_ext) {
            let distances_int = largest_center_to_vertex_distances_with_radius(
                0,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius_int,
            );
            let minmax_int_array: Box<[MinMax]> =
                to_shs_min_max_array(cone_radius_int, &distances_int);
            let distances_ext = largest_center_to_vertex_distances_with_radius(
                0,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius_ext,
            );
            let minmax_ext_array: Box<[MinMax]> =
                to_shs_min_max_array(cone_radius_ext, &distances_ext);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius_ext),
            );
            for h in 0..12 {
                self.ring_coverage_approx_recur(
                    0,
                    h,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    shs_cone_radius_int,
                    &minmax_int_array,
                    &minmax_ext_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            return bmoc_builder;
        }
        // Normal case
        let depth_start = best_starting_depth(cone_radius_ext);
        if depth_start >= self.depth {
            let shs_max = to_squared_half_segment(
                cone_radius_ext
                    + largest_center_to_vertex_distance_with_radius(
                        depth_start,
                        cone_lon,
                        cone_lat,
                        cone_radius_ext,
                    ),
            );
            let root_layer = get(depth_start);
            let center_hash_at_depth_start = root_layer.hash(cone_lon, cone_lat);
            let mut neigs: Vec<u64> = root_layer
                .neighbors(center_hash_at_depth_start, true)
                .values_vec()
                .iter()
                .filter(|neigh| {
                    for (lon, lat) in self.vertices(**neigh) {
                        if squared_half_segment(
                            lon - cone_lon,
                            lat - cone_lat,
                            lon.cos(),
                            cos_cone_lat,
                        ) > shs_cone_radius_int
                        {
                            return true;
                        }
                    }
                    false
                }) // remove cell fully in the internal raidus
                .map(h_to_h_and_shs(cone_lon, cone_lat, cos_cone_lat, root_layer))
                .filter(shs_lower_than(shs_max))
                .map(self.h_and_shs_to_lower_h(depth_start))
                .collect();
            neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
            neigs.dedup(); // remove duplicates (vector must be sorted first)
            let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
            for neig in neigs {
                bmoc_builder.push(self.depth, neig, false);
            }
            bmoc_builder
        } else {
            let distances_int = largest_center_to_vertex_distances_with_radius(
                depth_start,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius_int,
            );
            let minmax_int_array: Box<[MinMax]> =
                to_shs_min_max_array(cone_radius_int, &distances_int);
            let distances_ext = largest_center_to_vertex_distances_with_radius(
                depth_start,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius_ext,
            );
            let minmax_ext_array: Box<[MinMax]> =
                to_shs_min_max_array(cone_radius_ext, &distances_ext);
            let root_layer = get(depth_start);
            let root_center_hash = root_layer.hash(cone_lon, cone_lat);
            let neigs = root_layer.neighbors(root_center_hash, true);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius_ext),
            );
            for &root_hash in neigs.sorted_values().iter() {
                self.ring_coverage_approx_recur(
                    depth_start,
                    root_hash,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    shs_cone_radius_int,
                    &minmax_int_array,
                    &minmax_ext_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            bmoc_builder
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn ring_coverage_approx_recur<F>(
        &self,
        depth: u8,                // cell depth
        hash: u64,                // cell hash
        shs_computer: &F, // function to compute the squared_half_segment between a point and the cone center
        shs_int_cone_radius: f64, // squared_half_segment corresponding to the internal radius
        // shs_ext_cone_radius: f64,  // squared_half_segment corresponding to the external radius
        shs_int_minmax: &[MinMax], // min and max shs at various delta depth for the internal radius
        shs_ext_minmax: &[MinMax], // min and max shs at various delta depth for the external radius
        recur_depth: u8,           // conveniency delta depth index
        bmoc_builder: &mut BMOCBuilderUnsafe, // builder in which cells are appended (flag always set to true!)
    ) where
        F: Fn((f64, f64)) -> f64,
    {
        let center = get(depth).center(hash);
        let shs = shs_computer(center);
        let MinMax {
            min: min_int,
            max: max_int,
        } = shs_int_minmax[recur_depth as usize];
        let MinMax {
            min: min_ext,
            max: max_ext,
        } = shs_ext_minmax[recur_depth as usize];
        // We recall here that 'min_ext' and 'min_int' are lower bounds
        // while 'max_ext' and 'max_ext' are upper bounds
        if min_ext >= shs && shs >= max_int {
            // is fully in the ring (for sure)
            bmoc_builder.push(depth, hash, true);
        } else if max_ext >= shs && shs >= min_int {
            // may overlap (or not) or be fully in the ring
            if depth == self.depth {
                // If all 4 vertices are in the small radius, reject the cell,
                // i.e. if at least one vertex is not in the small radius, accept the cell
                for vertex_coo in self.vertices(hash) {
                    if shs_computer(vertex_coo) > shs_int_cone_radius {
                        bmoc_builder.push(depth, hash, false); // could be true, but we use a fast test only
                        break;
                    }
                }
            } else {
                //because 'min_int' and 'max_ext' are approximations, we may split cells fully
                // included in the MOC. We could check the 4 corers here, but we assume it is more costly
                // than post-merging (to be verified!)
                let hash = hash << 2;
                let depth = depth + 1;
                let recur_depth = recur_depth + 1;
                self.ring_coverage_approx_recur(
                    depth,
                    hash,
                    shs_computer,
                    shs_int_cone_radius,
                    shs_int_minmax,
                    shs_ext_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.ring_coverage_approx_recur(
                    depth,
                    hash | 1_u64,
                    shs_computer,
                    shs_int_cone_radius,
                    shs_int_minmax,
                    shs_ext_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.ring_coverage_approx_recur(
                    depth,
                    hash | 2_u64,
                    shs_computer,
                    shs_int_cone_radius,
                    shs_int_minmax,
                    shs_ext_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.ring_coverage_approx_recur(
                    depth,
                    hash | 3_u64,
                    shs_computer,
                    shs_int_cone_radius,
                    shs_int_minmax,
                    shs_ext_minmax,
                    recur_depth,
                    bmoc_builder,
                );
            }
        } // else fully out of the ring (for sure)
    }

    /// Returns a hierarchical view of the list of cells overlapped by the given cone.
    /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
    /// This is approximated: cell tagged as partially overlapped could be fully overlapped
    /// (but we guarantee that cells tagged as fully overlapped are really fully overlapped).
    ///
    /// The algorithm is fast but approximated: it may return false positive,
    /// i.e. cells which are near from the cone but do not overlap it.
    /// To control the approximation, see the method
    /// [cone_coverage_approx_custom](#method.cone_coverage_approx_custom)
    ///
    /// # Input
    /// - `cone_lon` the longitude of the center of the cone, in radians
    /// - `cone_lat` the latitude of the center of the cone, in radians
    /// - `cone_radius` the radius of the cone, in radians
    ///
    /// # Output
    /// - the list of cells overlapped by the given cone, in a BMOC (hierarchical view also telling
    ///   if a cell is fully or partially covered).
    ///
    /// # Example
    /// ```rust
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 3_u8;
    /// let nested3 = get(depth);
    ///
    /// let lon = 13.158329_f64.to_radians();
    /// let lat = -72.80028_f64.to_radians();
    /// let radius = 5.64323_f64.to_radians();
    ///
    /// let actual_res = nested3.cone_coverage_approx(lon, lat, radius);
    /// let expected_res: [u64; 10] = [512, 514, 515, 520, 521, 522, 544, 705, 708, 709];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///      assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn cone_coverage_approx(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
        self.cone_coverage_approx_internal(cone_lon, cone_lat, cone_radius)
            .to_bmoc_packing()
    }

    /// Returns a hierarchical view of the list of cells overlapped by the given cone.
    /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
    /// The algorithm is fast but approximated: it may return false positive,
    /// i.e. cells which are near from the cone but do not overlap it.
    /// To control the approximation, you can choose to perform the computations at a deeper depth
    /// using the `delta_depth` parameter.
    ///
    /// # Input
    /// - `delta_depth` the difference between this Layer depth and the depth at which the computations
    ///   are made (should remain quite small).
    /// - `cone_lon` the longitude of the center of the cone, in radians
    /// - `cone_lat` the latitude of the center of the cone, in radians
    /// - `cone_radius` the radius of the cone, in radians
    ///
    /// # Output
    /// - the list of cells overlapped by the given cone, in a BMOC (hierarchical view also telling
    ///   if a cell is fully or partially covered).
    ///
    /// # Panics
    /// If this layer depth + `delta_depth` > the max depth (i.e. 29)
    ///
    /// # Example
    /// ```rust
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 3_u8;
    /// let nested3 = get(depth);
    ///
    /// let lon = 13.158329_f64.to_radians();
    /// let lat = -72.80028_f64.to_radians();
    /// let radius = 5.64323_f64.to_radians();
    ///
    /// let actual_res = nested3.cone_coverage_approx_custom(2, lon, lat, radius);
    /// let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///     assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn cone_coverage_approx_custom(
        &self,
        delta_depth: u8,
        cone_lon: f64,
        cone_lat: f64,
        cone_radius: f64,
    ) -> BMOC {
        if delta_depth == 0 {
            self.cone_coverage_approx(cone_lon, cone_lat, cone_radius)
        } else {
            // TODO: change the algo not to put all cell in the MOC before pruning it
            get(self.depth + delta_depth)
                .cone_coverage_approx_internal(cone_lon, cone_lat, cone_radius)
                //.to_lower_depth_bmoc(self.depth)
                .to_lower_depth_bmoc_packing(self.depth)
        }
    }

    fn cone_coverage_approx_internal(
        &self,
        cone_lon: f64,
        cone_lat: f64,
        cone_radius: f64,
    ) -> BMOCBuilderUnsafe {
        // Special case: the full sky is covered
        if cone_radius >= PI {
            return self.allsky_bmoc_builder();
        }
        // Common variable
        let cos_cone_lat = cone_lat.cos();
        // Special case of very large radius: test the 12 base cells
        if !has_best_starting_depth(cone_radius) {
            let distances = largest_center_to_vertex_distances_with_radius(
                0,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius,
            );
            let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius),
            );
            for h in 0..12 {
                self.cone_coverage_approx_recur(
                    0,
                    h,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    &minmax_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            return bmoc_builder;
        }
        // Normal case
        let depth_start = best_starting_depth(cone_radius);
        if depth_start >= self.depth {
            let shs_max = to_squared_half_segment(
                cone_radius
                    + largest_center_to_vertex_distance_with_radius(
                        depth_start,
                        cone_lon,
                        cone_lat,
                        cone_radius,
                    ),
            );
            let root_layer = get(depth_start);
            let center_hash_at_depth_start = root_layer.hash(cone_lon, cone_lat);
            let mut neigs: Vec<u64> = root_layer
                .neighbors(center_hash_at_depth_start, true)
                .values_vec()
                .iter()
                .map(h_to_h_and_shs(cone_lon, cone_lat, cos_cone_lat, root_layer))
                .filter(shs_lower_than(shs_max))
                .map(self.h_and_shs_to_lower_h(depth_start))
                .collect();
            neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
            neigs.dedup(); // remove duplicates (vector must be sorted first)
            // return BMOC::create_unsafe(self.depth, neigs.into_boxed_slice());
            let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
            for neig in neigs {
                bmoc_builder.push(self.depth, neig, false);
            }
            bmoc_builder
        } else {
            let distances = largest_center_to_vertex_distances_with_radius(
                depth_start,
                self.depth + 1,
                cone_lon,
                cone_lat,
                cone_radius,
            );
            let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
            let root_layer = get(depth_start);
            let root_center_hash = root_layer.hash(cone_lon, cone_lat);
            let neigs = root_layer.neighbors(root_center_hash, true);
            let mut bmoc_builder = BMOCBuilderUnsafe::new(
                self.depth,
                self.n_moc_cell_in_cone_upper_bound(cone_radius),
            );
            for &root_hash in neigs.sorted_values().iter() {
                self.cone_coverage_approx_recur(
                    depth_start,
                    root_hash,
                    &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                    &minmax_array,
                    0,
                    &mut bmoc_builder,
                );
            }
            bmoc_builder
        }
    }

    fn cone_coverage_approx_recur<F>(
        &self,
        depth: u8,                            // cell depth
        hash: u64,                            // cell hash
        shs_computer: &F, // function to compute the squared_half_segment between a point and the cone center
        shs_minmax: &[MinMax], // min and max shs at various delta depth
        recur_depth: u8,  // conveniency delta depth index
        bmoc_builder: &mut BMOCBuilderUnsafe, // builder in which cells are appended
    ) where
        F: Fn((f64, f64)) -> f64,
    {
        let center = get(depth).center(hash);
        let shs = shs_computer(center);
        let MinMax { min, max } = shs_minmax[recur_depth as usize];
        if shs <= min {
            // Fast inclusion test
            bmoc_builder.push(depth, hash, true);
        } else if shs <= max {
            // Fast rejection test
            if depth == self.depth {
                bmoc_builder.push(depth, hash, false);
            } else {
                // because 'min' is a lower bound (approx), we may split cells fully included in the MOC
                // we could check the 4 corers here, we assume it is more costly than post-merging (to be verified!)
                let hash = hash << 2;
                let depth = depth + 1;
                let recur_depth = recur_depth + 1;
                self.cone_coverage_approx_recur(
                    depth,
                    hash,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_approx_recur(
                    depth,
                    hash | 1_u64,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_approx_recur(
                    depth,
                    hash | 2_u64,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
                self.cone_coverage_approx_recur(
                    depth,
                    hash | 3_u64,
                    shs_computer,
                    shs_minmax,
                    recur_depth,
                    bmoc_builder,
                );
            }
        }
    }

    /// cone_radius in radians
    /// TODO: find a better function!!
    #[inline]
    #[allow(dead_code)]
    fn ncell_in_cone_upper_bound(&self, cone_radius: f64) -> usize {
        // cell_area = 4 pi / ncell
        // cone_area = pi r^2
        // => area_ratio = cone_area / cell_area = r^2 * ncell / 4 = r^2 * (ncell >> 2)
        let mut area_ratio = pow2(cone_radius) * ((self.n_hash >> 2) as f64);
        area_ratio += 1.0_f64; // to be sure the minimum = 1
        // We want:
        //  - ratio = 1    --> n = 9
        //  - ratio = +inf --> n = 1.25 (125% => + 25%)
        //  - function shape = a/x + b
        // => a / 1 + b = 9  and a / +inf + b = 1.25
        // => b = 1.25 and a = 7.75
        let correction_factor = 1.25_f64 + 7.75_f64 / area_ratio;
        (correction_factor * area_ratio) as usize
    }

    #[inline]
    fn n_moc_cell_in_cone_upper_bound(&self, cone_radius: f64) -> usize {
        // 	NEW UPPER BOUND: supposedly more robust (and faster to compute)
        // But requires delta_depth to be given in parameter!
        // At lower resolution depth, max 9 cells (partially) overlapped:
        // => grid nside max = 3 * (2^DeltaDepth)
        // Worst case, for each nside max row (or col) at depth max:
        // - 2 (both sides) x (1 cell overllapping externaly + 1 cell overlapping internally + 1 no-fusioned internall cell)
        // - x2 to be conservative
        // 12 << delta_depth (delta_depth = diff between best starting depth and MOC depth max)

        // OLD UPPER BOUND (KEEP DURING THE TRANSITION): fails in rare circumstances,
        //   e.g. depth = 14, radius = 0.001, alpha = 0.002 ,delta = -1.3;
        // cell_area = 4 * pi / ncell = 4 * pi / (3 * 4 * nside^2) = pi / (3 * nside^2) =  pi * r^2
        // cell_radius = r = 1 / (sqrt(3) * nside)
        // As a very simple and naive rule, we take 4x the number of cells needed to cover
        // the cone external annulus
        // Annulus area = 4 pi ((R + r)^2 - R^2) = 4 pi (r^2 + 2rR)
        // N cells = 4 pi (r^2 + 2rR) / 4 pi r^2 = 1 + 2 R/r = 1 + 2 * sqrt(3) * nside * R
        const TWICE_SQRT_3: f64 = 2.0_f64 * 1.732_050_807_568_877_2_f64; // sqrt(3)
        4_usize * (1_usize + (self.nside as f64 * TWICE_SQRT_3 * cone_radius + 0.99_f64) as usize)
    }

    /*pub fn cone_coverage(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
      // Special case: the full sky is covered
      if cone_radius >= PI {
        return self.allsky_bmoc();
      }
      // Common variable
      let cos_cone_lat = cone_lat.cos();
      // Special case of very large radius: test the 12 base cells
      if !has_best_starting_depth(cone_radius) {

      }
      // Normal case

    }*/

    /// Returns a hierarchical view of the list of cells overlapped by the given elliptical cone.
    /// The BMOC also tells if the cell if fully or partially overlapped by the elliptical cone.
    /// The algorithm is approximated: it may return false positive,
    /// i.e. cells which are near from the elliptical cone but do not overlap it.
    /// To control the approximation, see the method
    /// [cone_coverage_approx_custom](#method.elliptical_cone_coverage_custom)
    ///
    /// # Input
    /// - `lon` the longitude of the center of the elliptical cone, in radians
    /// - `lat` the latitude of the center of the elliptical cone, in radians
    /// - `a` the semi-major axis of the elliptical cone, in radians
    /// - `b` the semi-minor axis of the elliptical cone, in radians
    /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
    ///
    /// # Output
    /// - the list of cells overlapped by the given elliptical cone, in a BMOC
    ///   (hierarchical view also telling if a cell is fully or partially covered).
    ///
    /// # Panics
    /// - if the semi-major axis is > PI/2
    ///
    /// # Example
    /// ```rust
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 3_u8;
    /// let nested3 = get(depth);
    ///
    /// let lon = 36.80105218_f64.to_radians();
    /// let lat = 56.78028536_f64.to_radians();
    /// let a = 14.93_f64.to_radians();
    /// let b = 4.93_f64.to_radians();
    /// let pa = 75.0_f64.to_radians();
    ///
    /// let actual_res = nested3.elliptical_cone_coverage(lon, lat, a, b, pa);
    /// let expected_res: [u64; 16] = [27, 30, 39, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 56, 57];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///     assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn elliptical_cone_coverage(&self, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
        self.elliptical_cone_coverage_internal(lon, lat, a, b, pa)
            .to_bmoc_packing()
    }

    /// Returns a hierarchical view of the list of cells overlapped by the given elliptical cone.
    /// The BMOC also tells if the cell if fully or partially overlapped by the elliptical cone.
    /// The algorithm is approximated: it may return false positive,
    /// i.e. cells which are near from the cone but do not overlap it.
    /// To control the approximation, you can choose to perform the computations at a deeper depth
    /// using the `delta_depth` parameter.
    ///
    /// # Input
    /// - `delta_depth` the difference between this Layer depth and the depth at which the computations
    ///   are made (should remain quite small).
    /// - `lon` the longitude of the center of the elliptical cone, in radians
    /// - `lat` the latitude of the center of the elliptical cone, in radians
    /// - `a` the semi-major axis of the elliptical cone, in radians
    /// - `b` the semi-minor axis of the elliptical cone, in radians
    /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
    ///
    /// # Output
    /// - the list of cells overlapped by the given elliptical cone, in a BMOC
    ///   (hierarchical view also telling if a cell is fully or partially covered).
    ///
    /// # Panics
    /// - if the semi-major axis is > PI/2
    /// - if this layer depth + `delta_depth` > the max depth (i.e. 29)
    ///
    pub fn elliptical_cone_coverage_custom(
        &self,
        delta_depth: u8,
        lon: f64,
        lat: f64,
        a: f64,
        b: f64,
        pa: f64,
    ) -> BMOC {
        if delta_depth == 0 {
            self.elliptical_cone_coverage(lon, lat, a, b, pa)
        } else {
            // TODO: change the algo not to put all cell in the MOC and pruning it
            get(self.depth + delta_depth)
                .elliptical_cone_coverage_internal(lon, lat, a, b, pa)
                .to_lower_depth_bmoc_packing(self.depth)
        }
    }

    pub fn elliptical_cone_coverage_internal(
        &self,
        lon: f64,
        lat: f64,
        a: f64,
        b: f64,
        pa: f64,
    ) -> BMOCBuilderUnsafe {
        if a >= FRAC_PI_2 {
            panic!("Unable to handle ellipses with a semi-major axis > PI/2");
        }
        // Common variable
        let sph_ellipse = EllipticalCone::new(lon, lat, a, b, pa);
        // Special case of very large radius: test the 12 base cells
        if !has_best_starting_depth(a) {
            let mut bmoc_builder =
                BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(a));
            let distances: Box<[f64]> =
                largest_center_to_vertex_distances_with_radius(0, self.depth + 1, lon, lat, a);
            for h in 0..12 {
                self.elliptical_cone_coverage_recur(
                    0,
                    h,
                    &sph_ellipse,
                    &distances,
                    0,
                    &mut bmoc_builder,
                );
            }
            return bmoc_builder;
        }
        // Normal case
        let depth_start = best_starting_depth(a);
        let root_layer = get(depth_start);
        let root_center_hash = root_layer.hash(lon, lat);
        // Small ellipse case
        if depth_start >= self.depth {
            let distance = largest_center_to_vertex_distance_with_radius(depth_start, lon, lat, a);
            let mut neigs: Vec<u64> = root_layer
                .neighbors(root_center_hash, true)
                .values_vec()
                .into_iter()
                .filter(|h| {
                    let (l, b) = root_layer.center(*h);
                    sph_ellipse.contains(l, b) || sph_ellipse.overlap_cone(l, b, distance)
                })
                .map(|h| h >> ((depth_start - self.depth) << 1)) // h_to_lower_depth
                .collect();
            neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
            neigs.dedup(); // remove duplicates (vector must be sorted first)
            let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
            for neig in neigs {
                bmoc_builder.push(self.depth, neig, false);
            }
            bmoc_builder
        } else {
            let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
                depth_start,
                self.depth + 1,
                lon,
                lat,
                a,
            );
            let neigs = root_layer.neighbors(root_center_hash, true);
            let mut bmoc_builder =
                BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(a));
            for &root_hash in neigs.sorted_values().iter() {
                self.elliptical_cone_coverage_recur(
                    depth_start,
                    root_hash,
                    &sph_ellipse,
                    &distances,
                    0,
                    &mut bmoc_builder,
                );
            }
            bmoc_builder
        }
    }

    fn elliptical_cone_coverage_recur(
        &self,
        depth: u8,
        hash: u64,
        ellipse: &EllipticalCone,
        distances: &[f64],
        recur_depth: u8,
        bmoc_builder: &mut BMOCBuilderUnsafe,
    ) {
        let (lon, lat) = get(depth).center(hash);
        let distance = distances[recur_depth as usize];
        /*eprintln!("d: {}; h: {}; lon: {}, lat: {}; dist: {}; contains: {}; overlap: {}",
        &depth, &hash, &lon.to_degrees(), &lat.to_degrees(), &distance.to_degrees(),
        &ellipse.contains_cone(lon, lat, distance),
        &ellipse.overlap_cone(lon, lat, distance));*/
        if ellipse.contains_cone(lon, lat, distance) {
            bmoc_builder.push(depth, hash, true);
        } else if ellipse.contains(lon, lat) || ellipse.overlap_cone(lon, lat, distance) {
            if depth == self.depth {
                let mut is_full = true;
                for (lon, lat) in self.vertices(hash).iter() {
                    is_full &= ellipse.contains(*lon, *lat); // Not sure computation not done if is_full==false, to be verfied
                }
                bmoc_builder.push(depth, hash, is_full);
            } else {
                let hash = hash << 2;
                let depth = depth + 1;
                let recur_depth = recur_depth + 1;
                self.elliptical_cone_coverage_recur(
                    depth,
                    hash,
                    ellipse,
                    distances,
                    recur_depth,
                    bmoc_builder,
                );
                self.elliptical_cone_coverage_recur(
                    depth,
                    hash | 1_u64,
                    ellipse,
                    distances,
                    recur_depth,
                    bmoc_builder,
                );
                self.elliptical_cone_coverage_recur(
                    depth,
                    hash | 2_u64,
                    ellipse,
                    distances,
                    recur_depth,
                    bmoc_builder,
                );
                self.elliptical_cone_coverage_recur(
                    depth,
                    hash | 3_u64,
                    ellipse,
                    distances,
                    recur_depth,
                    bmoc_builder,
                );
            }
        }
    }

    /// Returns a hierarchical view of the list of cells overlapped by the given box.
    /// The BMOC also tells if the cell if fully or partially overlapped by the box.
    ///
    /// # Input
    /// - `lon` the longitude of the center of the box, in radians
    /// - `lat` the latitude of the center of the box, in radians
    /// - `a` the semi-major axis of the box (half the box width), in radians
    /// - `b` the semi-minor axis of the box (half the box height), in radians
    /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
    ///
    /// # Output
    /// - the list of cells overlapped by the given box, in a BMOC
    ///   (hierarchical view also telling if a cell is fully or partially covered).
    ///
    /// # Panics
    /// - if `a` not in `]0, pi/2]`
    /// - if `b` not in `]0, a]`
    /// - if `pa` not in `[0, pi[`
    ///
    pub fn box_coverage(&self, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
        let center = Coo3D::from_sph_coo(lon, lat);
        let vertices = Layer::box2polygon(center.lon(), center.lat(), a, b, pa);
        self.custom_polygon_coverage(
            &vertices,
            &ContainsSouthPoleMethod::ControlPointIn(center),
            true,
        )
    }

    // # Panics
    // * if `lon` not in `[0, 2\pi[`
    // * if `lat` not in `[-\pi/2, \pi/2]`
    // * if `a` not in `]0, \pi/2]` or \`]0, \pi/2]`
    // * if `b` not in `]0, a]`
    // * if `pa` not in `[0, \pi[`
    fn box2polygon(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> Vec<(f64, f64)> {
        assert!(
            (0.0..TWICE_PI).contains(&lon),
            "Expected: lon in [0, 2pi[. Actual: {}",
            lon
        );
        assert!(
            (-HALF_PI..=HALF_PI).contains(&lat),
            "Expected: lat in [-pi/2, pi/2]. Actual: {}",
            lat
        );
        assert!(
            0.0 < a && a <= HALF_PI,
            "Expected: a in ]0, pi/2]. Actual: {}",
            a
        );
        assert!(0.0 < b && b <= a, "Expected: b in ]0, a]. Actual: {}", b);
        assert!(
            (0.0..PI).contains(&pa),
            "Expected: pa in [0, pi[. Actual: {}",
            pa
        );
        // Compute spherical coordinates
        let frame_rotation = RefToLocalRotMatrix::from_center(lon, lat);
        // By application of the Thales theorem, the new point has the property:
        //   sin(new_lat) / sin(dlat) = (cos(dlon) * cos(new_lat)) / cos(dlat)
        // With (imagine looking a the triangle in the (x, z) plane:
        // * old_x = cos(lat)
        // * new_x = cos(dlon) * cos(new_lat)
        // * old_z = sin(dlat)
        // * new_z = sin(new_lat)
        // Leading to:
        //   tan(new_lat) = cos(dlon) * tan(dlat)
        let lon = a;
        let (sin_lon, cos_lon) = lon.sin_cos();
        let lat = (cos_lon * b.tan()).atan();
        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_pa, cos_pa) = pa.sin_cos();
        // Rotation by the position angle
        // - upper right (before rotation by PA)
        let (x1, y1, z1) = (cos_lon * cos_lat, sin_lon * cos_lat, sin_lat);
        // - apply rotation (sin and cos are revere since theta = pi/2 - pa)
        let (y2, z2) = (y1 * sin_pa - z1 * cos_pa, y1 * cos_pa + z1 * sin_pa);
        let mut vertices = Vec::with_capacity(4);
        vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
        // - lower right (before rotation by PA) = (y1, -z1)
        let (y2, z2) = (y1 * sin_pa + z1 * cos_pa, y1 * cos_pa - z1 * sin_pa);
        vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
        // - lower left (before rotation by PA) = (-y1, -z1)
        let (y2, z2) = (-y1 * sin_pa + z1 * cos_pa, -y1 * cos_pa - z1 * sin_pa);
        vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
        // - upper left (before rotation by PA) = (-y1, z1)
        let (y2, z2) = (-y1 * sin_pa - z1 * cos_pa, -y1 * cos_pa + z1 * sin_pa);
        vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
        vertices
    }

    /// Returns a hierarchical view of the list of cells overlapped by the given polygon.
    /// Self intersecting polygons are supported.
    /// The BMOC also tells if the cell if fully or partially overlapped by the polygon.
    ///
    /// If you want the complementary solution, apply the NOT operator on the BMOC.
    ///
    /// This method supports both an *exact* (need more tests) and an *approximated* solution.
    /// The second one being faster (TODO: measure and provided perf differences as a function of the
    /// number of vertices in the polygon).
    ///
    /// The approximation is the following one:
    /// > when testing the intersection between a polygon segment and an HEALPix cell edge
    /// > we consider that each edge of the HEALPix cell is on a great-circle arc (which is not
    /// > true, especially a low resolutions).
    ///
    /// For the exact solution:
    /// > for each polygon segment, we first test if the segment contains a 'special point',
    /// > if it is the case, we add it to the list of cell number computed from each polygin vertex
    /// > A 'special point' is a point such that in the HEALPix projection Euclidean plane
    ///
    /// ```math
    /// \mathrm{d}\DeltaX(z) / \mathrm{d}Y = \pm 1
    /// ```
    ///
    /// # Input
    /// - `vertices` the list of vertices (in a slice) coordinates, in radians
    ///   `[(lon, lat), (lon, lat), ..., (lon, lat)]`
    /// - `exact_solution` if set
    ///
    /// # Output
    /// - the list of cells overlapped by the given polygon, in a BMOC (hierarchical view also telling
    ///   if a cell is fully or partially covered).
    ///
    /// # Example
    /// ```rust
    /// use healpix::compass_point::{Direction};
    /// use healpix::nested::{get, Layer};
    ///
    /// let depth = 3_u8;
    /// let nested3 = get(depth);
    ///
    /// let actual_res = nested3.polygon_coverage(&[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)], false);
    /// let expected_res: [u64; 8] = [304, 305, 306, 307, 308, 310, 313, 316];
    /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    ///     assert_eq!(h1, *h2);
    /// }
    /// ```
    pub fn polygon_coverage(&self, vertices: &[(f64, f64)], exact_solution: bool) -> BMOC {
        self.custom_polygon_coverage(vertices, &ContainsSouthPoleMethod::Default, exact_solution)
    }

    /// Same as [polygon_coverage](struct.Layer.html#method.polygon_coverage) but with a custom
    /// Polygon builder (in case the default behaviour is not satisfactory).
    pub fn custom_polygon_coverage(
        &self,
        vertices: &[(f64, f64)],
        south_pole_method: &ContainsSouthPoleMethod,
        exact_solution: bool,
    ) -> BMOC {
        let poly = Polygon::new_custom(
            vertices
                .iter()
                .map(|(lon, lat)| LonLat {
                    lon: *lon,
                    lat: *lat,
                })
                .collect::<Vec<LonLat>>()
                .into_boxed_slice(),
            south_pole_method,
        );

        let bounding_cone: Cone = Cone::bounding_cone(poly.vertices());
        let mut depth_start = 0;
        let neigs: Vec<u64> = if let ContainsSouthPoleMethod::Default = south_pole_method {
            if !has_best_starting_depth(bounding_cone.radius()) {
                // poly.must_contain(bounding_cone.center()); // Opposite of bounding cone out of?
                (0..12).collect()
            } else {
                depth_start = best_starting_depth(bounding_cone.radius()).min(self.depth);
                let root_layer = get(depth_start);
                let LonLat { lon, lat } = bounding_cone.center().lonlat();
                let center_hash = root_layer.hash(lon, lat);
                let mut neigs: Vec<u64> = root_layer.neighbors(center_hash, true).values_vec();
                neigs.sort_unstable();
                neigs
            }
        } else {
            (0..12).collect()
        };
        // Compute and sort the list of cells containing at least one polygon vertex
        let mut sorted_poly_vertices_hash = self.hashs_vec(poly.vertices());
        // Special treatment for the exact solution
        if exact_solution {
            let vertices = poly.vertices();
            let mut left = &vertices[vertices.len() - 1];
            for right in vertices {
                let special_lonlats = arc_special_points(left, right, 1.0e-14, 20);
                sorted_poly_vertices_hash.append(&mut self.hashs_vec(&special_lonlats));
                left = right;
            }
        }
        // Back to general case
        sorted_poly_vertices_hash.sort_unstable();
        sorted_poly_vertices_hash.dedup();
        let sorted_poly_vertices_hash = sorted_poly_vertices_hash.into_boxed_slice();
        // Build the list (removing duplicated) for all deltaDepth?
        let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 10_000_usize);
        for root_hash in neigs {
            self.polygon_coverage_recur(
                &mut bmoc_builder,
                depth_start,
                root_hash,
                &poly,
                &sorted_poly_vertices_hash,
            );
        }
        bmoc_builder.to_bmoc()
    }

    fn polygon_coverage_recur(
        &self,
        moc_builder: &mut BMOCBuilderUnsafe,
        depth: u8,
        hash: u64,
        poly: &Polygon,
        sorted_poly_vertices_hash: &[u64],
    ) {
        if is_in_list(depth, hash, self.depth, sorted_poly_vertices_hash) {
            if depth == self.depth {
                moc_builder.push(depth, hash, false);
            } else {
                let hash = hash << 2;
                let depth = depth + 1;
                self.polygon_coverage_recur(
                    moc_builder,
                    depth,
                    hash,
                    poly,
                    sorted_poly_vertices_hash,
                );
                self.polygon_coverage_recur(
                    moc_builder,
                    depth,
                    hash | 1,
                    poly,
                    sorted_poly_vertices_hash,
                );
                self.polygon_coverage_recur(
                    moc_builder,
                    depth,
                    hash | 2,
                    poly,
                    sorted_poly_vertices_hash,
                );
                self.polygon_coverage_recur(
                    moc_builder,
                    depth,
                    hash | 3,
                    poly,
                    sorted_poly_vertices_hash,
                );
            }
        } else {
            let (n_vertices_in_poly, poly_vertices) = n_vertices_in_poly(depth, hash, poly);
            if (n_vertices_in_poly > 0 && n_vertices_in_poly < 4)
                || has_intersection(poly, poly_vertices)
            {
                if depth == self.depth {
                    moc_builder.push(depth, hash, false);
                } else {
                    // I known, I don't like this code repetition, TODO: see how to remove it
                    let hash = hash << 2;
                    let depth = depth + 1;
                    self.polygon_coverage_recur(
                        moc_builder,
                        depth,
                        hash,
                        poly,
                        sorted_poly_vertices_hash,
                    );
                    self.polygon_coverage_recur(
                        moc_builder,
                        depth,
                        hash | 1,
                        poly,
                        sorted_poly_vertices_hash,
                    );
                    self.polygon_coverage_recur(
                        moc_builder,
                        depth,
                        hash | 2,
                        poly,
                        sorted_poly_vertices_hash,
                    );
                    self.polygon_coverage_recur(
                        moc_builder,
                        depth,
                        hash | 3,
                        poly,
                        sorted_poly_vertices_hash,
                    );
                }
            } else if n_vertices_in_poly == 4 {
                moc_builder.push(depth, hash, true);
            }
        }
    }

    fn hashs_vec<T: LonLatT>(&self, poly_vertices: &[T]) -> Vec<u64> {
        poly_vertices
            .iter()
            .map(|coo| self.hash(coo.lon(), coo.lat()))
            .collect::<Vec<u64>>()
    }

    fn allsky_bmoc_builder(&self) -> BMOCBuilderUnsafe {
        let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 12);
        bmoc_builder.push_all(0_u8, 0_u64, 12_u64, true);
        bmoc_builder
    }
}
