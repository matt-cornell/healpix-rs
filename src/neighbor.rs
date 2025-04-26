//! Functions relating to the neighboring edges.

use crate::dir::PosOrNeg;
use crate::dir::{Cardinal, Direction, Ordinal, map::DirectionMap};
use crate::zoc::{ZOrderCurve, get_zoc};
use crate::{HashBits, HashParts};

/// Returns (mod4 - 1) in [0, 2], and 3 if mod4 == 0 (i.e. the previous value in [0, 3] range)
#[inline]
fn prev(mod4: u8) -> u8 {
    debug_assert!(mod4 < 4);
    (((mod4 as i8) - 1) & 3) as u8
}

/// Returns (mod4 + 1) in [1, 3], and 0 if mod4 == 3 (i.e. the next value in [0, 3] range)
#[inline]
fn next(mod4: u8) -> u8 {
    debug_assert!(mod4 < 4);
    (mod4 + 1) & 3
}

/// Returns (mod4 + 2) in [2, 3], and 0 if mod4 == 2 and 1 if mod4 == 3 (i.e. the opposite value in [0, 3] range)
#[inline]
fn oppo(mod4: u8) -> u8 {
    debug_assert!(mod4 < 4);
    (mod4 + 2) & 3
}

/// Returns the input value: useless, just used to improve code legibility
#[inline]
fn iden(mod4: u8) -> u8 {
    debug_assert!(mod4 < 4);
    mod4
}

/// Compute the base cell from its (i, j) coordinates:
/// - i: index along the longitude axis ( = base_cell modulo 4)
/// - j: index along the latitude axis ( = base_cell / 4)
///   - = 0 for the cells covering the north polar cap
///   - = 1 for the cells with are only in the equatorial region
///   - = 2 for the cells covering the south polar cap
#[inline]
fn base_cell(i: u8, j: u8) -> u8 {
    debug_assert!(i < 4 && j < 3);
    (j << 2) + i
}

impl super::Layer {
    /// Retuns the hash value of the neighbor cell of the cell of given hash, in the given direction.
    /// If the cell do not have a neighbor in the given direction (which is the case of the
    /// eastmost and westmost cells in polar caps base cells and northmost and southmost cells of the
    /// equatorial region base cells), the return Option is None.
    ///
    /// # Input
    /// - `hash` the hash value of the cell we look for the neighbor
    /// - `direction` the direction of the neighbor we look for the cell number
    ///
    /// # Output
    /// - the cell number (hash value) of the neighbor of the given cell in the given direction
    ///   (`None` if their is no neighbor in the given direction) .
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    ///
    /// # Example
    /// ```rust
    /// use healpix::dir::Direction;
    /// use healpix::get;
    ///
    /// let depth = 0u8;
    /// let nested0 = get(depth);
    ///
    /// assert_eq!(5 , nested0.neighbor(4, Direction::E).unwrap());
    /// ```
    #[inline]
    pub fn neighbor(&self, hash: u64, direction: Direction) -> Option<u64> {
        let h_parts: HashParts = self.decode_hash(hash);
        self.neighbor_from_parts(h_parts.d0h, h_parts.i, h_parts.j, direction)
    }

    /// Returns the hash values of all the neighbor cells of the cell of given hash (8-connected).
    /// The given cell itself can be included (setting the `include_center` parameters to `true`).
    ///
    /// # Input
    /// - `hash` the hash value of the cell we look for the neighbors
    /// - `include_center` include (or not) the input cell in the Direction::C key of the returned map
    ///
    /// # Output
    /// - the cell number (hash value) of the neighbor of the given cell in the given direction
    ///   (`None` if their is no neighbor in the given direction) .
    ///
    /// # Panics
    /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
    ///
    /// # Example
    /// ```rust
    /// use healpix::dir::Direction;
    /// use healpix::get;
    ///
    /// let depth = 0u8;
    /// let nested0 = get(depth);
    ///
    /// assert_eq!(5 , *nested0.neighbors(4).get(Direction::E).unwrap());
    /// ```
    pub fn neighbors(&self, hash: u64) -> DirectionMap<u64> {
        self.check_hash(hash);
        let mut result_map = DirectionMap::new();
        let HashBits { d0h, i, j } = self.pull_bits_appart(hash);
        if self.is_in_base_cell_border(i, j) {
            self.edge_cell_neighbors(hash, &mut result_map);
        } else {
            self.inner_cell_neighbors(d0h, i, j, &mut result_map);
        }
        result_map
    }

    /// Returns the hash values of the cells which are on the internal border of a square of
    /// size `(1 + 2k) x (1 + 2k)` cells around the given cell (on the 2D projection plane,
    /// in the non-degenerated cases). This corresponds to an extension of the 8-connected definition.
    ///
    /// If that square does not cross base cells the result will contain `(1 + 2k) × (1 + 2k)` cells,
    /// if it does there will be fewer cells.
    ///
    /// The regular `neighbors` methods correspond to `k=1`, `k=0` returns the input cell.
    ///
    /// # Panics
    /// * if `k` is larger than or equals to `nside`.
    ///
    /// # Warning
    /// The algorithm we use works in the 2D projection plane.
    /// When the corner of a base cell has 2 neighbors instead of 3, the result shape on the sphere
    /// may be strange. Those corners are:
    /// * East and West corners of both north polar cap base cells (cells 0 to 3) and south polar cap
    ///   base cells (cells 8 to 11)
    /// * North and south corners of equatorial base cells (cell 4 to 7)
    pub fn kth_neighbors(&self, hash: u64, k: u32) -> Vec<u64> {
        match k {
            0 => [hash; 1].to_vec(),
            1 => self.neighbors(hash).into_values().collect(),
            k if k <= self.nside => {
                let HashParts { d0h, i, j } = self.decode_hash(hash);
                let mut result = Vec::with_capacity(((k << 1) as usize) << 2);
                self.kth_neighbors_internal(d0h, i, j, k, &mut result);
                result
            }
            _ => panic!(
                "The 'k' parameter is too large. Expected: <{}. Actual: {}.",
                self.nside, k
            ),
        }
    }

    /// Returns the hash values of all cells which are within the internal border of a square of
    /// size `(1 + 2k) x (1 + 2k)` cells around the given cell.
    ///
    /// If that square does not cross base cells the result will contain `(1 + 2k) × (1 + 2k)` cells,
    /// if it does there will be fewer cells.
    ///
    /// The regular `neighbors` methods correspond to `k=1`, `k=0` returns the input cell.
    ///
    /// # Panics
    /// * if `k` is larger than or equals to `nside`.
    ///
    /// # Warning
    /// The algorithm we use works in the 2D projection plane.
    /// When the corner of a base cell has 2 neighbors instead of 3, the result shape on the sphere
    /// may be strange. Those corners are:
    /// * East and West corners of both north polar cap base cells (cells 0 to 3) and south polar cap base cells (cells 8 to 11)
    /// * North and south corners of equatorial base cells (cell 4 to 7)
    pub fn kth_neighborhood(&self, hash: u64, k: u32) -> Vec<u64> {
        match k {
            0 => vec![hash],
            k if k <= self.nside => {
                let HashParts { d0h, i, j } = self.decode_hash(hash);
                let capacity = {
                    let k = k as usize;

                    ((k * (k + 1)) << 2) | 1
                };
                let mut result = Vec::with_capacity(capacity);
                result.push(hash);
                for r in 1..=k {
                    self.kth_neighbors_internal(d0h, i, j, r, &mut result);
                }

                result
            }
            _ => panic!(
                "The 'k' parameter is too large. Expected: <{}. Actual: {}.",
                self.nside, k
            ),
        }
    }

    fn kth_neighbors_internal(&self, d0h: u8, i: u32, j: u32, k: u32, result: &mut Vec<u64>) {
        if i >= k && j >= k && i + k < self.nside && j + k < self.nside {
            let xfrom = i - k;
            let xto = i + k;
            let yfrom = j - k;
            let yto = j + k;
            // S (inclusive) to W (exclusive)
            for y in yfrom..yto {
                result.push(self.build_hash_from_parts(d0h, xfrom, y));
            }
            // W (inclusive) to N (exclusive)
            for x in xfrom..xto {
                result.push(self.build_hash_from_parts(d0h, x, yto));
            }
            // N (inslusive) to E (exclusive)
            for y in (yfrom + 1..=yto).rev() {
                result.push(self.build_hash_from_parts(d0h, xto, y));
            }
            // E (inclusive) to S (exclusive)
            for x in (xfrom + 1..=xto).rev() {
                result.push(self.build_hash_from_parts(d0h, x, yfrom));
            }
        } else {
            let k = k as i32;
            let i = i as i32;
            let j = j as i32;

            let xfrom = i - k;
            let xto = i + k;
            let yfrom = j - k;
            let yto = j + k;

            // In this method, both i and j can be < 0 or >= nside
            let partial_compute =
                move |nside: i32, d0h: u8, i: i32, j: i32, k: i32, res: &mut Vec<u64>| {
                    let xfrom = i - k;
                    let xto = i + k;
                    let yfrom = j - k;
                    let yto = j + k;
                    // S (inclusive) to W (exclusive)
                    if (0..nside).contains(&xfrom) {
                        for y in yfrom.max(0)..yto.min(nside) {
                            res.push(self.build_hash_from_parts(d0h, xfrom as u32, y as u32));
                        }
                    }
                    // W (inclusive) to N (exclusive)
                    if (0..nside).contains(&yto) {
                        for x in xfrom.max(0)..xto.min(nside) {
                            res.push(self.build_hash_from_parts(d0h, x as u32, yto as u32));
                        }
                    }
                    // N (inslusive) to E (exclusive)
                    if (0..nside).contains(&xto) {
                        for y in ((yfrom + 1).max(0)..=yto.min(nside - 1)).rev() {
                            res.push(self.build_hash_from_parts(d0h, xto as u32, y as u32));
                        }
                    }
                    // E (inclusive) to S (exclusive)
                    if (0..nside).contains(&yfrom) {
                        for x in ((xfrom + 1).max(0)..=xto.min(nside - 1)).rev() {
                            res.push(self.build_hash_from_parts(d0h, x as u32, yfrom as u32));
                        }
                    }
                };

            let nside = self.nside as i32;
            let overflow_sw = xfrom < 0;
            let overflow_ne = xto >= nside;
            let overflow_se = yfrom < 0;
            let overflow_nw = yto >= nside;
            let overflow_s = overflow_sw && overflow_se;
            let overflow_w = overflow_sw && overflow_nw;
            let overflow_n = overflow_nw && overflow_ne;
            let overflow_e = overflow_se && overflow_ne;

            if overflow_s {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::S) {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_sw {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::SW)
                {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_w {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::W) {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_nw {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::NW)
                {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_n {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::N) {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_ne {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::NE)
                {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_e {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::E) {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            if overflow_se {
                if let Some((d0h, i, j)) = self.to_neighbor_base_cell_coo(d0h, i, j, Direction::SE)
                {
                    partial_compute(nside, d0h, i, j, k, result);
                }
            }
            partial_compute(nside, d0h, i, j, k, result);
        }
    }

    /// Same method as [neighbors](#method.neighbors) except that neighbors are appended
    /// to the given vector.
    pub fn append_bulk_neighbors(&self, hash: u64, dest: &mut Vec<u64>) {
        self.check_hash(hash);
        let HashBits { d0h, i, j } = self.pull_bits_appart(hash);
        if self.is_in_base_cell_border(i, j) {
            self.append_bulk_edge_cell_neighbors(hash, dest);
        } else {
            self.append_bulk_inner_cell_neighbors(d0h, i, j, dest);
        }
    }

    /// Returns the hash values corresponding to the internal bounds of the given `hash` at
    /// the hash depth + the given `delta_depth`:
    /// - the first quarter contains the southeast border (the z-order curve x-axis with y = 0);
    /// - the second quarter contains the northeast border (the z-order y-axis with x = xmax - 1);
    /// - the third quarter contains the northwest border (the z-order curve x-axis with y = ymax - 1);
    /// - the forth quarter contains the southwest border (the y-axis with x = 0).
    ///
    /// The hashes are ordered consecutively, starting from the south (x=0, y=0) cell in the
    /// anti-clokwise direction.
    ///
    /// # Input
    /// - `hash ` the hash for which we look for the internal bounds
    /// - `delta_depth` difference between the depth of the edge cells and the depth of the given cell  
    ///
    /// # Output
    /// - the cell numbers (hash values) of the given hash inner edge ordered consecutively,
    ///   starting from the south (x=0, y=0) cell in the anti-clokwise direction.
    ///
    pub fn internal_edge(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
        // Compute the x and y part masks for deltaDepth.
        let zoc = get_zoc(delta_depth);
        let twice_dd = delta_depth << 1;
        let x_max_bits = crate::x_mask(delta_depth);
        let y_max_bits = x_max_bits << 1;
        // Prepare hashes of depth of (self.depth + delta_depth), switching hash bits of 2 delta_depth to the left.
        hash <<= twice_dd;
        // Prepare filling the result.
        // am1 stands for a - 1, i.e. nSide - 1, i.e. the index of the last cell along the x or y-axis
        let am1 = (1_u32 << delta_depth) - 1; // 2^deltaDepth - 1
        let mut res: Vec<u64> = Vec::with_capacity((am1 << 2) as usize);
        // Southeast axis
        res.push(hash);
        for k in 1..am1 {
            let x_bits = zoc.i02h(k);
            res.push(hash | x_bits);
        }
        // Northeast axis
        res.push(hash | x_max_bits);
        for k in 1..am1 {
            res.push(hash | zoc.oj2h(k) | x_max_bits);
        }
        // Northwest axis
        res.push(hash | y_max_bits | x_max_bits);
        for k in 1..am1 {
            res.push(hash | y_max_bits | zoc.i02h(am1 - k));
        }
        // Southwest axis
        res.push(hash | y_max_bits);
        for k in 1..am1 {
            res.push(hash | zoc.oj2h(am1 - k));
        }
        res.into_boxed_slice()
    }

    /// Same as method [internal_edge](#method.internal_edge) except that the returned array is sorted.
    pub fn internal_edge_sorted(mut hash: u64, delta_depth: u8) -> Vec<u64> {
        if delta_depth == 0 {
            return vec![hash; 1];
        } else if delta_depth == 1 {
            hash <<= 2;
            return vec![hash, hash | 1, hash | 2, hash | 3];
        }
        // Compute the x and y part masks for deltaDepth.
        let zoc = get_zoc(delta_depth);
        let twice_dd = delta_depth << 1;
        let x_max_bits = crate::x_mask(delta_depth);
        let y_max_bits = x_max_bits << 1;
        // Prepare hashes of depth of (self.depth + delta_depth), switching hash bits of 2 delta_depth to the left.
        hash <<= twice_dd;
        // Set grid size (nSide inside the cell of depth this.depth)
        let nside = 1 << delta_depth;
        let am1 = nside - 1;
        let n_half_side = nside >> 1;
        // South sub-square (dividing in 4 sub-squares)
        let mut x = 1_u32;
        let mut lim = 2_u32;
        let mut k0 = 1_usize;
        let mut k1 = 2_usize;
        let mut k2 = (am1 + n_half_side) as usize;
        let mut k3 = ((am1 << 1) + n_half_side) as usize;
        let size = (am1 << 2) as usize;
        let mut result = vec![0_u64; size];
        // Set South corner (first element)
        result[0] = hash;
        // Set east corner
        result[k2 - 1] = hash | x_max_bits;
        // Set west corner
        result[k3 - 1] = hash | y_max_bits;
        // Set north corner (last element)
        result[size - k0] = hash | y_max_bits | x_max_bits;
        while x < n_half_side {
            // while (k < nHalfSize)
            x += 1;
            let xn = zoc.ij2h(x, nside - x); // x shuffled
            let xs = xn & x_max_bits;
            let xn = xn & y_max_bits;
            // South square, south east part
            result[k0] = hash | xs;
            k0 += 1;
            // South square, south west part
            result[k1] = hash | (xs << 1);
            k1 += 1;
            // East square, north east part
            result[k2] = hash | (xs << 1) | x_max_bits;
            k2 += 1;
            // West square, north west part
            result[k3] = hash | y_max_bits | xs;
            k3 += 1;
            // North square, north west
            result[size - k0] = hash | y_max_bits | (xn >> 1);
            // North square, north east
            result[size - k1] = hash | xn | x_max_bits;
            // West square, north west part
            result[size - k2] = hash | xn;
            // East square, nort east part
            result[size - k3] = hash | (xn >> 1);
            // Change k0, k1 and limit if x== limit.
            // The following lines of code are equivalent to:
            if x == lim {
                k0 = k1;
                k1 += lim as usize; // +2 +4 +8
                lim <<= 1; // 4 8 16 32 ...
            }
            // To be tested if they are faster (since no risk of branch miss-prediction):
            // probably true for small deltaDepth but not for large deltaDepth.
            /* let mut tmp = x & lim;  debug_assert!((x < lim && tmp == 0) || (x == lim && tmp == x));
            k0 += (tmp >> 1) as usize;
            k1 += tmp as usize;
            tmp -= x;           debug_assert!((x < lim            ) || (x == lim && tmp == 0));
            tmp = 1 >> tmp;     debug_assert!((x < lim && tmp == 0) || (x == lim && tmp == 1));
            lim <<= tmp;*/
        }
        result
    }

    /// Provides the list of all cells of depth this layer depth + the given `delta_depth`
    /// surrounding the cell of given hash value.  
    ///
    /// Here the result of both following codes:
    /// ![External edge depth 1, cells 10 and 11, delta_depth = +2](external_edge.png)
    ///
    /// ```rust
    /// use healpix::Layer;
    ///
    /// let depth = 1;
    /// let delta_depth = 2;
    ///
    /// let hash = 10;
    /// let actual_res = Layer::get(depth).external_edge(hash, delta_depth, true);
    /// let expected_res: [u64; 19] = [85, 87, 93, 95, 117, 138, 139, 142, 143, 154, 176, 178, 184, 186, 415, 437, 439, 445, 447];
    /// assert_eq!(actual_res, expected_res);
    ///
    /// let hash = 11;
    /// let actual_res = Layer::get(depth).external_edge(hash, delta_depth, true);
    /// let expected_res: [u64; 20] = [63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255];
    /// assert_eq!(actual_res, expected_res);
    ///
    /// ```
    pub fn external_edge(&self, hash: u64, delta_depth: u8, sorted: bool) -> Vec<u64> {
        let mut dest = Vec::with_capacity((4 + (self.nside << 2)) as usize); // 4 borders (nside) + 4 corners (1)
        self.append_external_edge(hash, delta_depth, sorted, &mut dest);
        dest
    }

    /// Similar to [external_edge](#method.external_edge) except that the result is appended to the
    /// provided  `dest` list.
    pub fn append_external_edge(
        &self,
        hash: u64,
        delta_depth: u8,
        sorted: bool,
        dest: &mut Vec<u64>,
    ) {
        //} -> Box<[u64]> {
        self.check_hash(hash);
        if delta_depth == 0 {
            if sorted {
                let mut tmp: Vec<u64> = Vec::with_capacity(8);
                self.append_bulk_neighbors(hash, &mut tmp);
                tmp.sort_unstable();
                dest.append(&mut tmp);
            } else {
                self.append_bulk_neighbors(hash, dest);
            }
            return;
        }

        // let mut edge = Vec::with_capacity((4 + (self.nside << 2)) as usize); // 4 borders (nside) + 4 corners (1)
        let h_bits: HashBits = self.pull_bits_appart(hash);
        if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
            // Not easy: opposite directions depends on base cell neighbors
            let mut neighbors = DirectionMap::new();
            self.edge_cell_neighbors(hash, &mut neighbors);
            let mut neighbors = neighbors.into_iter().collect::<arrayvec::ArrayVec<_, 8>>();
            if sorted {
                neighbors.sort_unstable_by_key(|v| v.1);
            }
            let h_parts: HashParts = self.decode_hash(hash);
            for (direction, hash_value) in neighbors {
                let dir_from_neig = if h_parts.d0h == self.h_2_d0h(hash_value) {
                    direction.opposite()
                } else if self.depth == 0 {
                    direction_from_neighbor(h_parts.d0h, direction)
                } else {
                    edge_cell_direction_from_neighbor(
                        h_parts.d0h,
                        self.direction_in_base_cell_border(h_bits.i, h_bits.j)
                            .expect("No for current cell"),
                        direction,
                    )
                };
                append_sorted_internal_edge_element(hash_value, delta_depth, dir_from_neig, dest);
            }
        } else {
            // Easy: always use the opposite direction
            let mut neighbors = DirectionMap::new();
            self.inner_cell_neighbors(h_bits.d0h, h_bits.i, h_bits.j, &mut neighbors);
            let mut neighbors = neighbors.into_iter().collect::<arrayvec::ArrayVec<_, 8>>();
            if sorted {
                neighbors.sort_unstable_by_key(|v| v.1);
            }
            for (direction, hash_value) in neighbors {
                append_sorted_internal_edge_element(
                    hash_value,
                    delta_depth,
                    direction.opposite(),
                    dest,
                );
            }
        }
        // edge.into_boxed_slice()
    }

    #[inline]
    fn is_in_base_cell_border(&self, i_in_base_cell_bits: u64, j_in_base_cell_bits: u64) -> bool {
        0_u64 == i_in_base_cell_bits
            || i_in_base_cell_bits == self.x_mask
            || 0_u64 == j_in_base_cell_bits
            || j_in_base_cell_bits == self.y_mask
    }

    #[inline]
    fn direction_in_base_cell_border(
        &self,
        i_in_base_cell_bits: u64,
        j_in_base_cell_bits: u64,
    ) -> Option<Direction> {
        let i = if 0_u64 == i_in_base_cell_bits {
            0
        } else if i_in_base_cell_bits == self.x_mask {
            2
        } else {
            1
        };
        let j = if 0_u64 == j_in_base_cell_bits {
            0
        } else if j_in_base_cell_bits == self.y_mask {
            2
        } else {
            1
        };
        // S SE E
        // SW C NE
        // W NW N
        Direction::try_from_index([4, 3, 2, 5, u8::MAX, 1, 6, 7, 0][3 * j + i])
    }

    fn inner_cell_neighbors(
        &self,
        d0h_bits: u64,
        i_in_d0h_bits: u64,
        j_in_d0h_bits: u64,
        result_map: &mut DirectionMap<u64>,
    ) {
        let ij = self.z_order_curve.h2ij(i_in_d0h_bits | j_in_d0h_bits);
        let i = self.z_order_curve.ij2i(ij);
        let j = self.z_order_curve.ij2j(ij);
        // Compute i-1 and j-1 bits.
        // Depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could result
        // in making fewer operations
        let ij = self.z_order_curve.ij2h(i - 1, j - 1);
        let im1_bits = ij & self.x_mask;
        let jm1_bits = ij & self.y_mask;
        // Compute i+1 and j+1 bits.
        // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
        // result in making fewer operations
        let ij = self.z_order_curve.ij2h(i + 1, j + 1);
        let ip1_bits = ij & self.x_mask;
        let jp1_bits = ij & self.y_mask;
        // Unrolled for loop an Direction enumset
        result_map.insert(Direction::S, d0h_bits | im1_bits | jm1_bits);
        result_map.insert(Direction::SE, d0h_bits | i_in_d0h_bits | jm1_bits);
        result_map.insert(Direction::E, d0h_bits | ip1_bits | jm1_bits);
        result_map.insert(Direction::SW, d0h_bits | im1_bits | j_in_d0h_bits);
        result_map.insert(Direction::NE, d0h_bits | ip1_bits | j_in_d0h_bits);
        result_map.insert(Direction::W, d0h_bits | im1_bits | jp1_bits);
        result_map.insert(Direction::NW, d0h_bits | i_in_d0h_bits | jp1_bits);
        result_map.insert(Direction::N, d0h_bits | ip1_bits | jp1_bits);
    }

    fn append_bulk_inner_cell_neighbors(
        &self,
        d0h_bits: u64,
        i_in_d0h_bits: u64,
        j_in_d0h_bits: u64,
        dest: &mut Vec<u64>,
    ) {
        let ij = self.z_order_curve.h2ij(i_in_d0h_bits | j_in_d0h_bits);
        let i = self.z_order_curve.ij2i(ij);
        let j = self.z_order_curve.ij2j(ij);
        // Compute i-1 and j-1 bits.
        // Depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could result
        // in making fewer operations
        let ij = self.z_order_curve.ij2h(i - 1, j - 1);
        let im1_bits = ij & self.x_mask;
        let jm1_bits = ij & self.y_mask;
        // Compute i+1 and j+1 bits.
        // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
        // result in making fewer operations
        let ij = self.z_order_curve.ij2h(i + 1, j + 1);
        let ip1_bits = ij & self.x_mask;
        let jp1_bits = ij & self.y_mask;
        // Unrolled for loop an Direction enumset
        dest.push(d0h_bits | im1_bits | jm1_bits);
        dest.push(d0h_bits | i_in_d0h_bits | jm1_bits);
        dest.push(d0h_bits | ip1_bits | jm1_bits);
        dest.push(d0h_bits | im1_bits | j_in_d0h_bits);
        dest.push(d0h_bits | ip1_bits | j_in_d0h_bits);
        dest.push(d0h_bits | im1_bits | jp1_bits);
        dest.push(d0h_bits | i_in_d0h_bits | jp1_bits);
        dest.push(d0h_bits | ip1_bits | jp1_bits);
    }

    fn edge_cell_neighbors(&self, hash: u64, result_map: &mut DirectionMap<u64>) {
        // Could have simply been edgeCellNeighbours(hash, EnumSet.allOf(Direction.class) result)
        // but we preferred to unroll the for loop.
        let HashParts { d0h, i, j } = self.decode_hash(hash);
        for dir in Direction::VALUES {
            if let Some(hash) = self.neighbor_from_parts(d0h, i, j, dir) {
                result_map.insert(dir, hash);
            }
        }
    }

    fn append_bulk_edge_cell_neighbors(&self, hash: u64, dest: &mut Vec<u64>) {
        let HashParts { d0h, i, j } = self.decode_hash(hash);
        dest.extend(
            Direction::VALUES
                .iter()
                .filter_map(|&dir| self.neighbor_from_parts(d0h, i, j, dir)),
        );
    }

    fn neighbor_from_parts(&self, d0h: u8, i: u32, j: u32, dir: Direction) -> Option<u64> {
        let i = (i as i32) + (dir.se() as i32);
        let j = (j as i32) + (dir.sw() as i32);
        let d0_neighbor_dir = Direction::from_sesw(
            self.neighbor_base_cell_offset(i),
            self.neighbor_base_cell_offset(j),
        );
        self.neighbor_from_shifted_coos(d0h, i as u32, j as u32, d0_neighbor_dir)
    }

    /*fn neighbor_base_cell_offset(&self, offset: i8, coo: u32) -> i8 {
      if coo == 0_u32 && offset == -1_i8 {
        -1_i8
      } else if coo == self.nside_minus_1 && offset == 1_i8 {
        1_i8
      } else {
        0_i8
      }
    }*/
    /// This method has a single input parameters `coo` which must be in `[-1, nside]`, and returns:
    /// - -1 if `coo` == -1
    /// -  0 if `coo` in `[0, nside[`
    /// -  1 if `coo` == nside
    #[inline]
    fn neighbor_base_cell_offset(&self, coo_in_base_cell: i32) -> PosOrNeg {
        debug_assert!(-1_i32 <= coo_in_base_cell && coo_in_base_cell <= (self.nside as i32));
        let offset = (coo_in_base_cell >> 31 | coo_in_base_cell >> self.depth) as i8;
        debug_assert!(
            (coo_in_base_cell == -1_i32 && offset == -1_i8)
                || (coo_in_base_cell == (self.nside as i32) && offset == 1_i8)
                || (-1_i32 < coo_in_base_cell
                    && coo_in_base_cell < (self.nside as i32)
                    && offset == 0_i8)
        );
        PosOrNeg::from_val(offset)
    }

    #[inline]
    fn neighbor_from_shifted_coos(
        &self,
        d0h: u8,
        i: u32,
        j: u32,
        base_cell_neighbor_dir: Option<Direction>,
    ) -> Option<u64> {
        if let Some(dir) = base_cell_neighbor_dir {
            let d0h_mod_4 = d0h & 0b11;
            match d0h >> 2 {
                0 => self.npc_neighbor(d0h_mod_4, i, j, dir),
                1 => self.eqr_neighbor(d0h_mod_4, i, j, dir),
                2 => self.spc_neighbor(d0h_mod_4, i, j, dir),
                _ => unreachable!("Base cell must be in [0, 12["),
            }
        } else {
            debug_assert!(i < self.nside && j < self.nside);
            Some(self.build_hash_from_parts(d0h, i, j))
        }
    }
    fn npc_neighbor(
        &self,
        d0h_mod_4: u8,
        i: u32,
        j: u32,
        base_cell_neighbor_dir: Direction,
    ) -> Option<u64> {
        let m = self.nside_minus_1;
        Some(match base_cell_neighbor_dir {
            Direction::S => self.build_hash_from_parts(base_cell(iden(d0h_mod_4), 2), m, m),
            Direction::SE => self.build_hash_from_parts(base_cell(next(d0h_mod_4), 1), i, m),
            Direction::SW => self.build_hash_from_parts(base_cell(iden(d0h_mod_4), 1), m, j),
            Direction::NE => self.build_hash_from_parts(base_cell(next(d0h_mod_4), 0), j, m),
            Direction::NW => self.build_hash_from_parts(base_cell(prev(d0h_mod_4), 0), m, i),
            Direction::N => self.build_hash_from_parts(base_cell(oppo(d0h_mod_4), 0), m, m),
            _ => None?,
        })
    }
    fn eqr_neighbor(
        &self,
        d0h_mod_4: u8,
        i: u32,
        j: u32,
        base_cell_neighbor_dir: Direction,
    ) -> Option<u64> {
        let m = self.nside_minus_1;
        Some(match base_cell_neighbor_dir {
            Direction::SE => self.build_hash_from_parts(base_cell(iden(d0h_mod_4), 2), i, m),
            Direction::E => self.build_hash_from_parts(base_cell(next(d0h_mod_4), 1), 0, m),
            Direction::SW => self.build_hash_from_parts(base_cell(prev(d0h_mod_4), 2), m, j),
            Direction::NE => self.build_hash_from_parts(base_cell(iden(d0h_mod_4), 0), 0, j),
            Direction::W => self.build_hash_from_parts(base_cell(prev(d0h_mod_4), 1), m, 0),
            Direction::NW => self.build_hash_from_parts(base_cell(prev(d0h_mod_4), 0), i, 0),
            _ => None?,
        })
    }
    fn spc_neighbor(
        &self,
        d0h_mod_4: u8,
        i: u32,
        j: u32,
        base_cell_neighbor_dir: Direction,
    ) -> Option<u64> {
        Some(match base_cell_neighbor_dir {
            Direction::S => self.build_hash_from_parts(base_cell(oppo(d0h_mod_4), 2), 0, 0),
            Direction::SE => self.build_hash_from_parts(base_cell(next(d0h_mod_4), 2), 0, i),
            Direction::SW => self.build_hash_from_parts(base_cell(prev(d0h_mod_4), 2), j, 0),
            Direction::NE => self.build_hash_from_parts(base_cell(next(d0h_mod_4), 1), 0, j),
            Direction::NW => self.build_hash_from_parts(base_cell(iden(d0h_mod_4), 1), i, 0),
            Direction::N => self.build_hash_from_parts(base_cell(iden(d0h_mod_4), 0), 0, 0),
            _ => None?,
        })
    }

    /// In the 2D projection plane, transform the given `(i, j)` coordinates in the frame defined by:
    /// * origin: South vertex of the input base cell
    /// * x-axis: South vertex to East vertex of the input base cell
    /// * y-axis: South vertex to West vertex of the input base cell
    ///   into coordinates in the frame attached to the neighbor cell of given direction (with respect to
    ///   the input base cell).
    /// # Params:
    /// * `new_base_cell_dir`: direction of the base cell in which we want the new coordinates, with
    ///   respect to the given `d0h` base cell.
    /// # Retuned params:
    /// * returns `None` if the base cell has no neighbor in the given direction.
    /// * `r.0`: base cell number of the cell in which the new coordinates are expressed.
    /// * `r.1`: coordinate along the `S to SE` axis in the new base cell.
    /// * `r.2`: coordinate along the `S to SW` axis in the new base cell.
    /// # Info:
    /// * input `i` and `j` coordinates are supposed to be in `[0, nside[` while output coordinates
    ///   may be negative and or larger that `nside`.
    pub fn to_neighbor_base_cell_coo(
        &self,
        d0h: u8,
        i: i32,
        j: i32,
        new_base_cell_dir: Direction,
    ) -> Option<(u8, i32, i32)> {
        let d0h_mod_4 = d0h & 0b11;
        match d0h >> 2 {
            0 => self.to_neighbor_base_cell_coo_npc(d0h_mod_4, i, j, new_base_cell_dir),
            1 => self.to_neighbor_base_cell_coo_eqr(d0h_mod_4, i, j, new_base_cell_dir),
            2 => self.to_neighbor_base_cell_coo_spc(d0h_mod_4, i, j, new_base_cell_dir),
            _ => unreachable!("Base cell must be in [0, 12["),
        }
    }
    fn to_neighbor_base_cell_coo_npc(
        &self,
        d0h_mod_4: u8,
        i: i32,
        j: i32,
        new_base_cell_dir: Direction,
    ) -> Option<(u8, i32, i32)> {
        let n = self.nside as i32;
        let m = (n << 1) - 1;
        match new_base_cell_dir {
            Direction::S => Some((base_cell(iden(d0h_mod_4), 2), n + i, n + j)),
            Direction::SE => Some((base_cell(next(d0h_mod_4), 1), i, n + j)),
            Direction::SW => Some((base_cell(iden(d0h_mod_4), 1), n + i, j)),
            Direction::NE => Some((base_cell(next(d0h_mod_4), 0), j, m - i)),
            Direction::NW => Some((base_cell(prev(d0h_mod_4), 0), m - j, i)),
            Direction::N => Some((base_cell(oppo(d0h_mod_4), 0), m - i, m - j)),
            _ => None,
        }
    }
    fn to_neighbor_base_cell_coo_eqr(
        &self,
        d0h_mod_4: u8,
        i: i32,
        j: i32,
        new_base_cell_dir: Direction,
    ) -> Option<(u8, i32, i32)> {
        let n = self.nside as i32;
        match new_base_cell_dir {
            Direction::SE => Some((base_cell(iden(d0h_mod_4), 2), i, n + j)),
            Direction::E => Some((base_cell(next(d0h_mod_4), 1), i - n, n + j)),
            Direction::SW => Some((base_cell(prev(d0h_mod_4), 2), n + i, j)),
            Direction::NE => Some((base_cell(iden(d0h_mod_4), 0), i - n, j)),
            Direction::W => Some((base_cell(prev(d0h_mod_4), 1), n + i, j - n)),
            Direction::NW => Some((base_cell(prev(d0h_mod_4), 0), i, j - n)),
            _ => None,
        }
    }
    fn to_neighbor_base_cell_coo_spc(
        &self,
        d0h_mod_4: u8,
        i: i32,
        j: i32,
        new_base_cell_dir: Direction,
    ) -> Option<(u8, i32, i32)> {
        let n = self.nside as i32;
        match new_base_cell_dir {
            Direction::S => Some((base_cell(oppo(d0h_mod_4), 2), -i - 1, -j - 1)),
            Direction::SE => Some((base_cell(next(d0h_mod_4), 2), -j - 1, i)),
            Direction::SW => Some((base_cell(prev(d0h_mod_4), 2), j, -i - 1)),
            Direction::NE => Some((base_cell(next(d0h_mod_4), 1), i - n, j)),
            Direction::NW => Some((base_cell(iden(d0h_mod_4), 1), i, j - n)),
            Direction::N => Some((base_cell(iden(d0h_mod_4), 0), i - n, j - n)),
            _ => None,
        }
    }
}

/// Returns the direction of a cell on the inner edge of the given base cell from its neighbor
/// located at the given direction in a different base cell.
/// # Inputs
/// - `base_cell` the base cell containing the sub-cell we are looking for the direction from its
///   neighbor in the given `neighbor_direction`
/// - `inner_direction` the direction of the sub-cell in the edge of the given base cell
/// - `neighbor_direction` direction of the neighbor of the sub-cell from which we are looking
///   at the direction of the sub-cell
///   
pub fn edge_cell_direction_from_neighbor(
    base_cell: u8,
    inner_direction: Direction,
    neighbor_direction: Direction,
) -> Direction {
    match base_cell >> 2 {
        // <=> basce_cell / 4
        0 => npc_egde_direction_from_neighbor(inner_direction, neighbor_direction),
        1 => eqr_edge_direction_from_neighbor(inner_direction, neighbor_direction),
        2 => spc_edge_direction_from_neighbor(inner_direction, neighbor_direction),
        _ => panic!("Base cell must be in [0, 12["),
    }
}

fn npc_egde_direction_from_neighbor(
    inner_direction: Direction,
    neighbor_direction: Direction,
) -> Direction {
    match neighbor_direction {
        Direction::E => match inner_direction {
            Direction::N | Direction::NE => Direction::N,
            Direction::E => panic!("No neighbor in direction {:?}", neighbor_direction),
            Direction::S | Direction::SE => neighbor_direction.opposite(),
            _ => unreachable!(),
        },
        Direction::W => match inner_direction {
            Direction::N | Direction::NW => Direction::N,
            Direction::W => panic!("No neighbor in direction {:?}", neighbor_direction),
            Direction::S | Direction::SW => neighbor_direction.opposite(),
            _ => unreachable!(),
        },
        Direction::NE => {
            assert!([Direction::N, Direction::E, Direction::NE].contains(&inner_direction));
            Direction::NW
        }
        Direction::NW => {
            assert!([Direction::N, Direction::W, Direction::NW].contains(&inner_direction));
            Direction::NE
        }
        Direction::N => match inner_direction {
            Direction::N => Direction::N,
            Direction::E | Direction::NE => Direction::W,
            Direction::W | Direction::NW => Direction::E,
            _ => unreachable!(),
        },
        _ => neighbor_direction.opposite(),
    }
}

fn eqr_edge_direction_from_neighbor(
    _inner_direction: Direction,
    neighbor_direction: Direction,
) -> Direction {
    neighbor_direction.opposite()
}

fn spc_edge_direction_from_neighbor(
    inner_direction: Direction,
    neighbor_direction: Direction,
) -> Direction {
    match neighbor_direction {
        Direction::E => match inner_direction {
            Direction::S | Direction::SE => Direction::S,
            Direction::E => panic!("No neighbor in direction {:?}", neighbor_direction),
            Direction::N | Direction::NE => neighbor_direction.opposite(),
            _ => unreachable!(),
        },
        Direction::W => match inner_direction {
            Direction::S | Direction::SW => Direction::S,
            Direction::W => panic!("No neighbor in direction {:?}", neighbor_direction),
            Direction::N | Direction::NW => neighbor_direction.opposite(),
            _ => unreachable!(),
        },
        Direction::SE => {
            assert!([Direction::S, Direction::E, Direction::SE].contains(&inner_direction));
            Direction::SW
        }
        Direction::SW => {
            assert!([Direction::S, Direction::W, Direction::SW].contains(&inner_direction));
            Direction::SE
        }
        Direction::S => match inner_direction {
            Direction::S => Direction::S,
            Direction::E | Direction::SE => Direction::W,
            Direction::W | Direction::SW => Direction::E,
            _ => unreachable!(),
        },
        _ => neighbor_direction.opposite(),
    }
}

/// Returns the direction of the given base cell from its neighbor base cell located
/// in the given direction.
/// # Panics
/// If the base cell has no neighbor in the given direction (i.e. N/S for equatorial cells
/// and E/W for polar caps cells)
pub fn direction_from_neighbor(base_cell: u8, neighbor_direction: Direction) -> Direction {
    match base_cell >> 2 {
        // <=> basce_cell / 4
        0 => npc_direction_from_neighbor(neighbor_direction),
        1 => eqr_direction_from_neighbor(neighbor_direction),
        2 => spc_direction_from_neighbor(neighbor_direction),
        _ => panic!("Base cell must be in [0, 12["),
    }
}

fn npc_direction_from_neighbor(neighbor_direction: Direction) -> Direction {
    match neighbor_direction {
        Direction::E | Direction::W => {
            panic!("No neighbor in direction {:?}", neighbor_direction)
        }
        Direction::NE => Direction::NW,
        Direction::NW => Direction::NE,
        Direction::N => Direction::N,
        _ => neighbor_direction.opposite(),
    }
}

fn eqr_direction_from_neighbor(neighbor_direction: Direction) -> Direction {
    match neighbor_direction {
        Direction::S | Direction::N => {
            panic!("No neighbor in direction {:?}", neighbor_direction)
        }
        _ => neighbor_direction.opposite(),
    }
}

fn spc_direction_from_neighbor(neighbor_direction: Direction) -> Direction {
    match neighbor_direction {
        Direction::E | Direction::W => {
            panic!("No neighbor in direction {:?}", neighbor_direction)
        }
        Direction::S => Direction::S,
        Direction::SE => Direction::SW,
        Direction::SW => Direction::SE,
        _ => neighbor_direction.opposite(),
    }
}

fn append_sorted_internal_edge_element(
    hash: u64,
    delta_depth: u8,
    direction: Direction,
    result: &mut Vec<u64>,
) {
    if direction.is_cardinal() {
        result.push(internal_corner(
            hash,
            delta_depth,
            direction.unwrap_cardinal(),
        ));
    } else {
        append_internal_edge_part(hash, delta_depth, direction.unwrap_ordinal(), result);
    }
}

/// Returns the hash values of the cells of depth this layer depth + the given `delta_depth`
/// located in the internal edge of given direction in the given cell.
///
/// # Info
/// The returned vector is sorted.
///
/// ```rust
/// use healpix::dir::Ordinal;
/// use healpix::neighbor::internal_edge_part;
///
///
/// let delta_depth = 1;
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::SE), vec![0, 1]);
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::SW), vec![0, 2]);
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::NE), vec![1, 3]);
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::NW), vec![2, 3]);
///
/// let delta_depth = 2;
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::SE), vec![ 0,  1,  4,  5]);
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::SW), vec![ 0,  2,  8, 10]);
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::NE), vec![ 5,  7, 13, 15]);
/// assert_eq!(internal_edge_part(0, delta_depth, Ordinal::NW), vec![10, 11, 14, 15]);
/// ```
pub fn internal_edge_part(hash: u64, delta_depth: u8, direction: Ordinal) -> Vec<u64> {
    let mut dest = Vec::new();
    append_internal_edge_part(hash, delta_depth, direction, &mut dest);
    dest
}

/// Same as [internal_edge_part](fn.internal_edge_part.html) except that the result is appended
/// to the given vec.
pub fn append_internal_edge_part(
    hash: u64,
    delta_depth: u8,
    direction: Ordinal,
    result: &mut Vec<u64>,
) {
    fn append_internal_edge_southeast(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
        hash <<= delta_depth << 1;
        let nside = 1_u32 << delta_depth; // 2^deltaDepth
        for x in 0..nside {
            result.push(hash | get_zoc(delta_depth).i02h(x));
        }
    }
    fn append_internal_edge_southwest(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
        hash <<= delta_depth << 1;
        let nside = 1_u32 << delta_depth; // 2^deltaDepth
        for y in 0..nside {
            result.push(hash | get_zoc(delta_depth).oj2h(y));
        }
    }
    fn append_internal_edge_northeast(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
        hash <<= delta_depth << 1;
        let nside = 1_u32 << delta_depth; // 2^deltaDepth
        let x_bits = get_zoc(delta_depth).i02h(nside - 1);
        for y in 0..nside {
            result.push(hash | get_zoc(delta_depth).oj2h(y) | x_bits);
        }
    }
    fn append_internal_edge_northwest(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
        hash <<= delta_depth << 1;
        let nside = 1_u32 << delta_depth; // 2^deltaDepth
        let y_bits = get_zoc(delta_depth).oj2h(nside - 1);
        for x in 0..nside {
            result.push(hash | get_zoc(delta_depth).i02h(x) | y_bits);
        }
    }
    match direction {
        Ordinal::SE => append_internal_edge_southeast(hash, delta_depth, result),
        Ordinal::SW => append_internal_edge_southwest(hash, delta_depth, result),
        Ordinal::NE => append_internal_edge_northeast(hash, delta_depth, result),
        Ordinal::NW => append_internal_edge_northwest(hash, delta_depth, result),
    }
}

/// Returns the hash value of the cell of depth this layer depth + the given `delta_depth`
/// located in the corner of given direction in the given cell.
/// ```rust
/// use healpix::dir::Cardinal;
/// use healpix::neighbor::internal_corner;
///
/// let delta_depth = 1;
/// assert_eq!(0 , internal_corner(0, delta_depth, Cardinal::S));
/// assert_eq!(1 , internal_corner(0, delta_depth, Cardinal::E));
/// assert_eq!(2 , internal_corner(0, delta_depth, Cardinal::W));
/// assert_eq!(3 , internal_corner(0, delta_depth, Cardinal::N));
/// ```
pub fn internal_corner(hash: u64, delta_depth: u8, direction: Cardinal) -> u64 {
    match direction {
        Cardinal::S => hash << (delta_depth << 1),
        Cardinal::E => (hash << (delta_depth << 1)) | crate::x_mask(delta_depth),
        Cardinal::N => (hash << (delta_depth << 1)) | crate::xy_mask(delta_depth),
        Cardinal::W => (hash << (delta_depth << 1)) | crate::y_mask(delta_depth),
    }
}
