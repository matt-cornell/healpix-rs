use super::*;
use crate::bmoc::Bmoc;
use crate::coords::Degrees;
use crate::dir::map::DirectionMap;
use Direction::*;
use std::f64::consts::PI;

#[test]
fn testok_hash_d0() {
    let layer = get(0);
    assert_eq!(4_u64, layer.hash(Degrees(0.0, 0.0)));
    assert_eq!(5_u64, layer.hash(Degrees(90.0, 0.0)));
    assert_eq!(6_u64, layer.hash(Degrees(180.0, 0.0)));
    assert_eq!(7_u64, layer.hash(Degrees(270.0, 0.0)));

    assert_eq!(0_u64, layer.hash(Degrees(45.0, 41.0)));
    assert_eq!(1_u64, layer.hash(Degrees(135.0, 41.0)));
    assert_eq!(2_u64, layer.hash(Degrees(225.0, 41.0)));
    assert_eq!(3_u64, layer.hash(Degrees(315.0, 41.0)));

    assert_eq!(8_u64, layer.hash(Degrees(45.0, -41.0)));
    assert_eq!(9_u64, layer.hash(Degrees(135.0, -41.0)));
    assert_eq!(10_u64, layer.hash(Degrees(225.0, -41.0)));
    assert_eq!(11_u64, layer.hash(Degrees(315.0, -41.0)));
}

#[test]
fn testok_hash() {
    let layer = get(3);
    let hash = layer.hash(Degrees(333.5982493968911, -25.919634217871433));
    assert_eq!(735_u64, hash);
}

#[test]
#[allow(clippy::approx_constant)]
fn testok_hash_2() {
    let layer = get(0);
    // ra = 179.99999999999998633839 deg
    // de = -48.13786699999999889561 deg
    let hash = layer.hash([3.141592653589793, -0.8401642740371252]);
    // The difference comes from the fact that:
    //   ra * 4/pi = 4.0 instead of 3.99999999999999969640 due to numerical approximations.
    // In v1, it is a particular case (k=3) in depth0_bits, but we do not handle all of them
    if hash == 9_u64 {
        assert_eq!(9_u64, hash); // with hash_v1
    } else {
        assert_eq!(10_u64, hash); // with hash_v2
    }
}

#[test]
#[allow(clippy::approx_constant)]
fn testok_hash_3() {
    let layer = get(0);
    // ra = 89.99999999999999889877 deg
    // de = -42.68491599999999973256 deg
    let hash = layer.hash([1.5707963267948966, -0.7449923251372079]);
    // The difference comes from the fact that:
    //   ra * 4/pi = 2.0 instead of 1.99999999999999997552 due to numerical approximations.
    // In v1, it is a particular case (k=3) in depth0_bits, but we do not handle all of them
    if hash == 8_u64 {
        assert_eq!(8_u64, hash);
    } else {
        assert_eq!(9_u64, hash);
    }
}

#[test]
fn testok_hash_4() {
    let layer = get(3);
    let ra = 180.0_f64;
    let dec = -45.85_f64;
    let hash = layer.hash(Degrees(ra, dec));
    assert_eq!(682_u64, hash);
}

#[test]
fn testok_neighbor_d0() {
    let layer = get(0);
    // North polar cap
    // - 0
    assert_eq!(2_u64, layer.neighbor(0, N).unwrap());
    assert_eq!(1_u64, layer.neighbor(0, NE).unwrap());
    assert_eq!(3_u64, layer.neighbor(0, NW).unwrap());
    assert_eq!(8_u64, layer.neighbor(0, S).unwrap());
    assert_eq!(5_u64, layer.neighbor(0, SE).unwrap());
    assert_eq!(4_u64, layer.neighbor(0, SW).unwrap());
    assert_eq!(None, layer.neighbor(0, E));
    assert_eq!(None, layer.neighbor(0, W));
    // - 1
    assert_eq!(3_u64, layer.neighbor(1, N).unwrap());
    assert_eq!(2_u64, layer.neighbor(1, NE).unwrap());
    assert_eq!(0_u64, layer.neighbor(1, NW).unwrap());
    assert_eq!(9_u64, layer.neighbor(1, S).unwrap());
    assert_eq!(6_u64, layer.neighbor(1, SE).unwrap());
    assert_eq!(5_u64, layer.neighbor(1, SW).unwrap());
    assert_eq!(None, layer.neighbor(1, E));
    assert_eq!(None, layer.neighbor(1, W));
    // - 2
    assert_eq!(0_u64, layer.neighbor(2, N).unwrap());
    assert_eq!(3_u64, layer.neighbor(2, NE).unwrap());
    assert_eq!(1_u64, layer.neighbor(2, NW).unwrap());
    assert_eq!(10_u64, layer.neighbor(2, S).unwrap());
    assert_eq!(7_u64, layer.neighbor(2, SE).unwrap());
    assert_eq!(6_u64, layer.neighbor(2, SW).unwrap());
    assert_eq!(None, layer.neighbor(2, E));
    assert_eq!(None, layer.neighbor(2, W));
    // - 3
    assert_eq!(1_u64, layer.neighbor(3, N).unwrap());
    assert_eq!(0_u64, layer.neighbor(3, NE).unwrap());
    assert_eq!(2_u64, layer.neighbor(3, NW).unwrap());
    assert_eq!(11_u64, layer.neighbor(3, S).unwrap());
    assert_eq!(4_u64, layer.neighbor(3, SE).unwrap());
    assert_eq!(7_u64, layer.neighbor(3, SW).unwrap());
    assert_eq!(None, layer.neighbor(3, E));
    assert_eq!(None, layer.neighbor(3, W));
    // Equatorial region
    // - 4
    assert_eq!(None, layer.neighbor(4, N));
    assert_eq!(0_u64, layer.neighbor(4, NE).unwrap());
    assert_eq!(3_u64, layer.neighbor(4, NW).unwrap());
    assert_eq!(None, layer.neighbor(4, S));
    assert_eq!(8_u64, layer.neighbor(4, SE).unwrap());
    assert_eq!(11_u64, layer.neighbor(4, SW).unwrap());
    assert_eq!(5_u64, layer.neighbor(4, E).unwrap());
    assert_eq!(7_u64, layer.neighbor(4, W).unwrap());
    // - 5
    assert_eq!(None, layer.neighbor(5, N));
    assert_eq!(1_u64, layer.neighbor(5, NE).unwrap());
    assert_eq!(0_u64, layer.neighbor(5, NW).unwrap());
    assert_eq!(None, layer.neighbor(5, S));
    assert_eq!(9_u64, layer.neighbor(5, SE).unwrap());
    assert_eq!(8_u64, layer.neighbor(5, SW).unwrap());
    assert_eq!(6_u64, layer.neighbor(5, E).unwrap());
    assert_eq!(4_u64, layer.neighbor(5, W).unwrap());
    // - 6
    assert_eq!(None, layer.neighbor(6, N));
    assert_eq!(2_u64, layer.neighbor(6, NE).unwrap());
    assert_eq!(1_u64, layer.neighbor(6, NW).unwrap());
    assert_eq!(None, layer.neighbor(6, S));
    assert_eq!(10_u64, layer.neighbor(6, SE).unwrap());
    assert_eq!(9_u64, layer.neighbor(6, SW).unwrap());
    assert_eq!(7_u64, layer.neighbor(6, E).unwrap());
    assert_eq!(5_u64, layer.neighbor(6, W).unwrap());
    // - 7
    assert_eq!(None, layer.neighbor(7, N));
    assert_eq!(3_u64, layer.neighbor(7, NE).unwrap());
    assert_eq!(2_u64, layer.neighbor(7, NW).unwrap());
    assert_eq!(None, layer.neighbor(7, S));
    assert_eq!(11_u64, layer.neighbor(7, SE).unwrap());
    assert_eq!(10_u64, layer.neighbor(7, SW).unwrap());
    assert_eq!(4_u64, layer.neighbor(7, E).unwrap());
    assert_eq!(6_u64, layer.neighbor(7, W).unwrap());
    //  South polar cap
    // - 8
    assert_eq!(0_u64, layer.neighbor(8, N).unwrap());
    assert_eq!(5_u64, layer.neighbor(8, NE).unwrap());
    assert_eq!(4_u64, layer.neighbor(8, NW).unwrap());
    assert_eq!(10_u64, layer.neighbor(8, S).unwrap());
    assert_eq!(9_u64, layer.neighbor(8, SE).unwrap());
    assert_eq!(11_u64, layer.neighbor(8, SW).unwrap());
    assert_eq!(None, layer.neighbor(8, E));
    assert_eq!(None, layer.neighbor(8, W));
    // - 9
    assert_eq!(1_u64, layer.neighbor(9, N).unwrap());
    assert_eq!(6_u64, layer.neighbor(9, NE).unwrap());
    assert_eq!(5_u64, layer.neighbor(9, NW).unwrap());
    assert_eq!(11_u64, layer.neighbor(9, S).unwrap());
    assert_eq!(10_u64, layer.neighbor(9, SE).unwrap());
    assert_eq!(8_u64, layer.neighbor(9, SW).unwrap());
    assert_eq!(None, layer.neighbor(9, E));
    assert_eq!(None, layer.neighbor(9, W));
    // - 10
    assert_eq!(2_u64, layer.neighbor(10, N).unwrap());
    assert_eq!(7_u64, layer.neighbor(10, NE).unwrap());
    assert_eq!(6_u64, layer.neighbor(10, NW).unwrap());
    assert_eq!(8_u64, layer.neighbor(10, S).unwrap());
    assert_eq!(11_u64, layer.neighbor(10, SE).unwrap());
    assert_eq!(9_u64, layer.neighbor(10, SW).unwrap());
    assert_eq!(None, layer.neighbor(10, E));
    assert_eq!(None, layer.neighbor(10, W));
    // - 11
    assert_eq!(3_u64, layer.neighbor(11, N).unwrap());
    assert_eq!(4_u64, layer.neighbor(11, NE).unwrap());
    assert_eq!(7_u64, layer.neighbor(11, NW).unwrap());
    assert_eq!(9_u64, layer.neighbor(11, S).unwrap());
    assert_eq!(8_u64, layer.neighbor(11, SE).unwrap());
    assert_eq!(10_u64, layer.neighbor(11, SW).unwrap());
    assert_eq!(None, layer.neighbor(11, E));
    assert_eq!(None, layer.neighbor(11, W));
}

#[test]
fn testok_neighbors_d0() {
    let layer = get(0);
    // North polar cap
    check_equals(layer.neighbors(0_u64), [1, 2, 3, 4, 5, 8]);
    check_equals(layer.neighbors(1_u64), [0, 2, 3, 5, 6, 9]);
    check_equals(layer.neighbors(2_u64), [0, 1, 3, 6, 7, 10]);
    check_equals(layer.neighbors(3_u64), [0, 1, 2, 4, 7, 11]);
    // Equatorial region
    check_equals(layer.neighbors(4_u64), [0, 3, 5, 7, 8, 11]);
    check_equals(layer.neighbors(5_u64), [0, 1, 4, 6, 8, 9]);
    check_equals(layer.neighbors(6_u64), [1, 2, 5, 7, 9, 10]);
    check_equals(layer.neighbors(7_u64), [2, 3, 4, 6, 10, 11]);
    //  South polar cap
    check_equals(layer.neighbors(8_u64), [0, 4, 5, 9, 10, 11]);
    check_equals(layer.neighbors(9_u64), [1, 5, 6, 8, 10, 11]);
    check_equals(layer.neighbors(10_u64), [2, 6, 7, 8, 9, 11]);
    check_equals(layer.neighbors(11_u64), [3, 4, 7, 8, 9, 10]);
}

#[test]
fn test_ok_neighbors_t1() {
    let depth = 2_u8;
    let hash = 130;

    let layer = get(depth);
    check_equals_all(
        layer.neighbors(hash),
        [128, 129, 131, 136, 137, 176, 177, 180],
    );
}

fn check_equals(map: DirectionMap<u64>, array: [u64; 6]) {
    let mut values = map.into_values().collect::<Vec<_>>();
    values.sort_unstable();
    assert_eq!(array.len(), values.len(), "array lengths differ");
    assert_eq!(values, array);
    array
        .iter()
        .zip(&values)
        .for_each(|(h1, h2)| assert_eq!(h1, h2));
}

fn check_equals_all(map: DirectionMap<u64>, array: [u64; 8]) {
    let mut values = map.into_values().collect::<Vec<_>>();
    values.sort_unstable();
    assert_eq!(array.len(), values.len(), "array lengths differ");
    assert_eq!(values, array);
    array
        .iter()
        .zip(&values)
        .for_each(|(h1, h2)| assert_eq!(h1, h2));
}

/*#[test]
fn testok_cone_flat() {
  let res = cone_overlap_flat(3, 13.158329_f64.to_radians(), -72.80028_f64.to_radians(), 5.64323_f64.to_radians());
  // println!("@@@@@@@@@@@@@@ {:?}", &res);
}*/

#[test]
fn testok_external_edge_struct_v1() {
    let depth = 1;
    let hash = 10;
    let delta_depth = 2;
    // draw moc 3/117,138,139,142,143,437,439,445,447,176, 178, 184, 186,85, 87, 93, 95,415,154
    // let e = external_edge_struct(depth, hash, delta_depth);
    // println!("{:?}", &e);
    let actual_res = get(depth).external_edge(hash, delta_depth, true);
    let expected_res: [u64; 19] = [
        85, 87, 93, 95, 117, 138, 139, 142, 143, 154, 176, 178, 184, 186, 415, 437, 439, 445, 447,
    ];
    assert_eq!(actual_res, expected_res);
    // println!("{:?}", &actual_res);
    // assert_eq!(expected_res.len(), actual_res.len());
    // for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
    //     assert_eq!(h1, h2);
    // }
}

#[test]
fn testok_external_edge_struct_v2() {
    let depth = 1;
    let hash = 11;
    let delta_depth = 2;
    // draw moc 3/63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255
    let actual_res = get(depth).external_edge(hash, delta_depth, true);
    let expected_res: [u64; 20] = [
        63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251,
        254, 255,
    ];
    assert_eq!(actual_res, expected_res);
    // println!("{:?}", &actual_res);
    // for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
    //     assert_eq!(h1, h2);
    // }
    // assert_eq!(expected_res.len(), actual_res.len());
}

#[test]
fn testok_external_edge_struct_v3() {
    let depth = 0;
    let hash = 0;
    let delta_depth = 2;
    // draw moc 3/63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255
    let actual_res = get(depth).external_edge(hash, delta_depth, true);
    let expected_res: [u64; 18] = [
        26, 27, 30, 31, 47, 53, 55, 61, 63, 69, 71, 77, 79, 90, 91, 94, 95, 143,
    ];
    // let expected_res: [u64; 20] = [63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255];
    assert_eq!(actual_res, expected_res);
    //println!("{:?}", &actual_res);
    // for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
    //     assert_eq!(h1, h2);
    // }
    // assert_eq!(expected_res.len(), actual_res.len());
}

#[test]
fn testok_cone_approx_bmoc() {
    // let res = cone_overlap_approx(5, 0.01, 0.02, 0.05);
    // let res = cone_overlap_approx(6, 160.771389_f64.to_radians(), 64.3813_f64.to_radians(), 0.8962_f64.to_radians());
    let actual_res =
        get(3).cone_coverage_approx(Degrees(13.158329, -72.80028), 5.64323_f64.to_radians());
    let expected_res: [u64; 10] = [512, 514, 515, 520, 521, 522, 544, 705, 708, 709];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
    /*for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
}

#[test]
fn testok_cone_approx_custom_bmoc_2() {
    // let res = cone_overlap_approx(5, 0.01, 0.02, 0.05);
    // let res = cone_overlap_approx(6, 160.771389_f64.to_radians(), 64.3813_f64.to_radians(), 0.8962_f64.to_radians());
    let actual_res = get(3).cone_coverage_approx_custom(
        2,
        Degrees(36.80105218, 56.78028536),
        14.93_f64.to_radians(),
    );
    let expected_res: [u64; 22] = [
        26, 27, 30, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 56, 57, 58, 59, 60,
    ];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
    /*println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
}

#[test]
fn testok_cone_approx_custom_bmoc() {
    let actual_res = get(3).cone_coverage_approx_custom(
        2,
        Degrees(13.158329, -72.80028),
        5.64323_f64.to_radians(),
    );
    let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_cone_approx_custom_bmoc_v2() {
    let actual_res = get(6).cone_coverage_approx_custom(
        3,
        Degrees(8.401978, 84.675171),
        0.0008_f64.to_radians(),
    );
    /*for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    let expected_res: [u64; 1] = [4075];
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }

    println!("----------");

    let actual_res = get(6).cone_coverage_approx_custom(
        3,
        Degrees(8.401978, 84.675171),
        0.0004_f64.to_radians(),
    );
    /*for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_cone_approx_custom_bmoc_dbg() {
    let actual_res =
        get(2).cone_coverage_approx_custom(1, Degrees(20.0, 0.0), 50.0_f64.to_radians());
    let expected_res: [u64; 50] = [
        0, 1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 49, 52, 53, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
        75, 76, 77, 78, 79, 82, 88, 89, 90, 91, 94, 131, 134, 135, 136, 137, 138, 139, 140, 141,
        142, 143, 181, 183, 189,
    ];
    assert_eq!(
        actual_res.into_flat_iter().collect::<Vec<_>>(),
        expected_res
    );
    // assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    // for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_conecenter_bmoc_dbg() {
    let depth = 4_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius = 5.64323_f64;
    let actual_res = get(depth).cone_coverage_centers(Degrees(lon, lat), radius.to_radians());
    let expected_res: [u64; 7] = [2058, 2059, 2080, 2081, 2082, 2083, 2088];
    // println!("draw red circle({} {} {}deg)", lon, lat, radius);
    // to_aladin_moc(&actual_res);
    println!("entries: {:?}", actual_res.entries());
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_conefullin_bmoc_dbg() {
    let depth = 4_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius = 5.64323_f64;
    let actual_res = get(depth).cone_coverage_fullin(Degrees(lon, lat), radius.to_radians());
    let expected_res: [u64; 2] = [2081, 2082];
    println!("draw red circle({lon} {lat} {radius}deg)");
    to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_ring_bmoc_dbg() {
    let depth = 4_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius_int = 5.64323_f64;
    let radius_ext = 10.0_f64;
    let actual_res = get(depth).ring_coverage_approx(
        Degrees(lon, lat),
        radius_int.to_radians(),
        radius_ext.to_radians(),
    );
    let expected_res: [u64; 40] = [
        2050, 2051, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061, 2062, 2063, 2080, 2083, 2084,
        2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2094, 2176, 2177, 2178, 2817, 2820, 2821,
        2822, 2823, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2880,
    ];
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_int);
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_ext);
    // to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_ring_bmoc() {
    let depth = 5_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius_int = 5.64323_f64;
    let radius_ext = 10.0_f64;
    let actual_res = get(depth).ring_coverage_approx_custom(
        2,
        Degrees(lon, lat),
        radius_int.to_radians(),
        radius_ext.to_radians(),
    );
    let expected_res: [u64; 99] = [
        8202, 8203, 8206, 8207, 8218, 8224, 8225, 8226, 8227, 8228, 8229, 8230, 8231, 8232, 8233,
        8234, 8236, 8237, 8239, 8240, 8241, 8242, 8243, 8246, 8248, 8249, 8250, 8251, 8252, 8254,
        8320, 8333, 8335, 8336, 8337, 8338, 8339, 8340, 8342, 8344, 8345, 8346, 8347, 8348, 8355,
        8356, 8357, 8358, 8359, 8360, 8361, 8362, 8363, 8364, 8365, 8366, 8367, 8368, 8369, 8370,
        8376, 8704, 8705, 8706, 8707, 8708, 11280, 11281, 11283, 11284, 11285, 11286, 11287, 11292,
        11293, 11328, 11329, 11330, 11331, 11332, 11333, 11334, 11335, 11336, 11337, 11340, 11341,
        11344, 11345, 11346, 11347, 11348, 11349, 11350, 11351, 11352, 11353, 11520, 11521,
    ];
    // For visual verification (cargo test testok_ring_bmoc -- --nocapture)
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_int);
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_ext);
    // to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    assert_eq!(actual_res.flat_iter().collect::<Vec<_>>(), expected_res);
    // for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

fn to_aladin_moc(bmoc: &dyn Bmoc) {
    print!("draw moc {}/", bmoc.max_depth());
    for cell in bmoc.flat_iter() {
        print!("{cell}, ");
    }
}

#[test]
fn testok_polygone_approx() {
    let actual_res = get(3).polygon_coverage(&[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)], false);
    /*let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }*/
    println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
        println!("@@@@@ cell a: {cell:?}");
    }
}

#[test]
fn testok_polygone_exact_npc() {
    // Aladin: draw polygon(65.11781779000003, 85.012424, 89.70533626000001, 87.06130188, 60.23667431000001, 85.609882)
    let depth = 6;
    let mut vertices = [
        (65.11781779000003, 85.012424),
        (89.70533626000001, 87.06130188),
        (60.23667431000001, 85.609882),
    ];
    let expected_res_approx: [u64; 2] = [4062, 4063];
    let expected_res_exact: [u64; 4] = [4062, 4063, 4084, 4085];

    to_radians(&mut vertices);

    let actual_res_approx = get(depth).polygon_coverage(&vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx
        .flat_iter()
        .zip(expected_res_approx.iter())
    {
        assert_eq!(h1, *h2);
    }

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_npc_2() {
    // Aladin: draw polygon(359.70533626,+87.06130188, 330.23667431,+85.60988200, 335.11781779,+85.01242400)
    let depth = 6;
    let mut vertices = [
        (359.70533626, 87.06130188),
        (330.23667431, 85.60988200),
        (335.11781779, 85.01242400),
    ];
    let expected_res_approx: [u64; 2] = [16350, 16351];
    let expected_res_exact: [u64; 4] = [16350, 16351, 16372, 16373];

    to_radians(&mut vertices);

    let actual_res_approx = get(depth).polygon_coverage(&vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx
        .flat_iter()
        .zip(expected_res_approx.iter())
    {
        assert_eq!(h1, *h2);
    }

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_npc_3() {
    // Aladin:  draw polygon(224.86211710,+78.10924662, 176.91129363 +83.92878811, 135.81578643,+78.24840426, 200.73574863,+73.58038790)
    let depth = 3;
    let mut vertices = [
        (224.86211710, 78.10924662),
        (176.91129363, 83.92878811),
        (135.81578643, 78.24840426),
        (200.73574863, 73.58038790),
    ];
    let expected_res_approx: [u64; 5] = [119, 125, 187, 188, 190];
    let expected_res_exact: [u64; 7] = [119, 125, 127, 187, 188, 190, 191];

    // Pb pour cell 127:
    // from (180, +83) --> (224, +78)=> 191
    // from (135, +78) --> (176, +83) => 127

    to_radians(&mut vertices);

    let actual_res_approx = get(depth).polygon_coverage(&vertices, false);
    // println!("draw moc 3/ {:?}", actual_res_approx.to_flat_array());
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx
        .flat_iter()
        .zip(expected_res_approx.iter())
    {
        assert_eq!(h1, *h2);
    }

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    // println!("draw moc 3/ {:?}", actual_res_exact.to_flat_array());
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_spc() {
    // Aladin: draw polygon(359.70533626,-87.06130188, 330.23667431,-85.60988200, 335.11781779,-85.01242400)
    let depth = 6;
    let mut vertices = [
        (359.70533626, -87.06130188),
        (330.23667431, -85.60988200),
        (335.11781779, -85.01242400),
    ];
    let expected_res_approx: [u64; 2] = [45072, 45074];
    let expected_res_exact: [u64; 4] = [45061, 45063, 45072, 45074];

    to_radians(&mut vertices);

    let actual_res_approx = get(depth).polygon_coverage(&vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx
        .flat_iter()
        .zip(expected_res_approx.iter())
    {
        assert_eq!(h1, *h2);
    }

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_eqr() {
    // In Aladin: draw polygon(180.08758393,-41.43289179, 191.00310758,-29.99207687, 181.59160475,-34.21976170)
    let depth = 3;
    let mut vertices = [
        (180.08758393, -41.43289179),
        (191.00310758, -29.99207687),
        (181.59160475, -34.219761700),
    ];
    let expected_res_approx: [u64; 2] = [384, 385];
    let expected_res_exact: [u64; 4] = [384, 385, 682, 683];

    to_radians(&mut vertices);

    let actual_res_approx = get(depth).polygon_coverage(&vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx
        .flat_iter()
        .zip(expected_res_approx.iter())
    {
        assert_eq!(h1, *h2);
    }

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_2() {
    // In Aladin: draw polygon(174.75937396073138, -49.16744206799886, 185.24062603926856, -49.16744206799887, 184.63292896369916, -42.32049830486584, 175.3670710363009, -42.32049830486584)
    let depth = 3; // 10

    let mut vertices = [
        (174.75937396073138, -49.16744206799886),
        (185.24062603926856, -49.16744206799887),
        (184.63292896369916, -42.32049830486584),
        (175.3670710363009, -42.32049830486584),
    ];
    let expected_res_exact: [u64; 4] = [596, 597, 680, 682];

    to_radians(&mut vertices);

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);

    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res_exact.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/

    assert!(actual_res_exact.deep_size() > 0);

    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_3() {
    // In Aladin: draw polygon(359.05552155,+04.00678949, 003.48359255,-04.16855897, 000.85415475,-04.30100561)
    let depth = 7; // 7

    let mut vertices = [
        (359.05552155, 4.00678949),
        (3.48359255, -4.16855897),
        (0.85415475, -4.30100561),
    ];
    let expected_res_exact: [u64; 86] = [
        69479, 69484, 69485, 69486, 69487, 69488, 69489, 69490, 69491, 69492, 69494, 69496, 69497,
        69498, 69499, 69500, 69501, 69502, 69503, 69572, 69573, 69574, 69575, 69581, 69584, 69585,
        69586, 69587, 69588, 69589, 69590, 69591, 69592, 69593, 69594, 69595, 69596, 69597, 69598,
        69599, 69617, 69619, 69620, 69621, 69622, 69623, 69628, 69629, 69631, 72322, 72328, 72330,
        72352, 72353, 72354, 72355, 72360, 72361, 72362, 72363, 72364, 72366, 75093, 77824, 77825,
        77826, 77827, 77828, 77830, 77831, 77833, 77835, 77836, 77837, 77838, 77839, 77860, 77861,
        77863, 77869, 77872, 77874, 77880, 77882, 77883, 77969,
    ];

    to_radians(&mut vertices);

    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);

    to_aladin_moc(&actual_res_exact);

    println!("@@@@@ FLAT VIEW");
    for cell in actual_res_exact.flat_iter() {
        println!("@@@@@ cell a: {cell:?}");
    }

    assert!(actual_res_exact.deep_size() > 0);

    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_4() {
    // In Aladin: draw polygon(353.8156714, -56.33202193, 6.1843286, -56.33202193, 5.27558041, -49.49378172, 354.72441959, -49.49378172)
    let depth = 3;
    let mut vertices = [
        (353.8156714, -56.33202193),
        (6.1843286, -56.33202193),
        (5.27558041, -49.49378172),
        (354.72441959, -49.49378172),
    ];
    let expected_res_exact: [u64; 4] = [546, 552, 721, 724];
    to_radians(&mut vertices);
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_polygone_exact_5() {
    // In Aladin: draw polygon(,269.0489904481256, 50.93107840310317,256.0142359254377,
    // 13.990475156963315,243.27946433651516, -7.489608522084967,224.32882610807906,
    // -24.613170006428422,181.0522052001408, -32.62716788567398,196.13865549198817,
    // -48.194574021576855,204.11410185161034, -61.679032617555755,204.22079738688976,
    // -75.90259740652621,452.0599233625658, -81.01038799438487,291.81039030880413,
    // -59.94991219586425,298.63694405367926, -35.41386619268488,307.97031962597856,
    // -11.966387215315455,330.9614415037472, 19.1371582387119,310.5212786376428,
    // 20.26388923574062,296.2930605411064, 25.849021337885926,283.4739050648204,
    // 34.80711164805169,269.0489904481256, 50.93107840310317)
    let depth = 6;
    let vertices = [
        (4.695790732486566, 0.8889150097255261),
        (4.4682913488764395, 0.24417985540748033),
        (4.2460276551603116, -0.1307183283958091),
        (3.915276622719796, -0.42958085596528983),
        (3.159957098738856, -0.5694515052059677),
        (3.4232653287700523, -0.8411539982726408),
        (3.5624631270616547, -1.0765021986213243),
        (3.5643253154494583, -1.3247502355595913),
        (7.889934078990009, -1.4138979988201017),
        (5.093052102418385, -1.0463233540993349),
        (5.212197941830804, -0.6180885659230597),
        (5.375096075892637, -0.2088528564758103),
        (5.776366851387001, 0.33400642074068165),
        (5.4196187097295985, 0.35367158642311125),
        (5.171289457253199, 0.45115053076437905),
        (4.947552986866946, 0.6074987013677716),
    ];
    // let expected_res_exact: [u64; 4] =[546, 552, 721, 724];
    // to_radians(&mut vertices);
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    to_aladin_moc(&actual_res_exact);
    /*assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }*/
    /*to_degrees(&mut vertices);
    print!("\ndraw polygon(");
    for (lon, lat) in vertices.iter() {
      print!(",{}, {}", lon, lat);
    }
    println!(")");*/
}

#[test]
fn testok_polygone_exact_6() {
    // In Aladin: draw polygon(268.84102386, -45.73283624, 299.40278164, 38.47909742, 66.0951825, 38.76894431, 96.66953469, -45.43470273)
    let depth = 0;
    let mut vertices = [
        (268.84102386, -45.73283624),
        (299.40278164, 38.47909742),
        (66.0951825, 38.76894431),
        (96.66953469, -45.43470273),
    ];
    let expected_res_exact: [u64; 9] = [0, 3, 4, 5, 7, 8, 9, 10, 11];
    to_radians(&mut vertices);
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
    /*to_degrees(&mut vertices);
    print!("\ndraw polygon(");
    for (lon, lat) in vertices.iter() {
      print!(",{}, {}", lon, lat);
    }
    println!(")");*/
}

#[test]
fn testok_polygone_exact_fxp() {
    let depth = 10;
    let mut vertices = [
        (000.9160848095, 01.0736381331),
        (001.5961114529, -00.7062969568),
        (359.6079412529, 01.1296198985),
        (358.7836886856, -01.3663552354),
        (001.8720201899, 00.4097184220),
        (358.4159783831, 00.2376811155),
        (359.8319515193, -01.2824324848),
        (358.7798765255, 00.9896544935),
        (001.5440798843, 00.8056786162),
    ];
    // 10/4806091 10/4806094
    to_radians(&mut vertices);
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    // to_aladin_moc(&actual_res_exact);
    for h in actual_res_exact.flat_iter() {
        assert!(
            h != 4806091 || h != 4806094,
            "Contains 10/4806091 or 10/4806094"
        );
    }
}

#[test]
fn testok_polygone_aladinlite() {
    use std::f64::consts::TAU;
    let depth = 9;
    let vertices = [
        (-1.3220905656465063 + TAU, -1.2576425633152113),
        (-1.322092139407076 + TAU, -1.2576425633127768),
        (-1.3220921394018748 + TAU, -1.2576422856351774),
        (-1.3220905656426547 + TAU, -1.257642285637612),
    ];
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    to_aladin_moc(&actual_res_exact);
    /*for h in actual_res_exact.flat_iter() {
      if h == 4806091 || h == 4806094 {
        assert!(false, "Contains 10/4806091 or 10/4806094")
      }
    }*/
}

/*
Just used to test visually, everything was OK
#[test]
fn testok_polygon_aladinlitev3(){
  let depth = 2;
  let vertices = [
    (2.762765958009174, -1.1558461939389928), (1.8511822776529407, -1.0933557485463217),
    (1.2223938385258124, -0.9257482282529843), (0.7433555354994513, -0.7588156605217202),
    (0.3057973131669505, -0.6330807489082455), (-0.13178353851370578, -0.5615526603403407),
    (-0.5856691165836956, -0.5400327446156468), (-0.4015158856736461, -0.255385904573062),
    (-0.36447295312493533, 0.055767062790012194), (-0.46118312909272624, 0.3439014598749883),
    (-0.7077120187325645, 0.5509825177593931), (-1.0770630524112423, 0.6032198101231719),
    (-1.4301394011313524, 0.46923928877261123), (-1.4878998343615297, 0.8530569092752573),
    (-1.9282441740007021, 1.1389659504677638), (-2.7514745722416727, 1.0909408969626568),
    (-3.043766621461535, 0.7764645583461347), (-3.0294959096529364, 0.4089232125719977),
    (-2.878305764522211, 0.05125802264111997), (3.0734739115502143, 0.057021963147195695),
    (2.7821939201916965, -0.05879878325187485), (2.5647762800185414, -0.27363879379223127),
    (2.4394176374317205, -0.5545388864809319), (2.444623710270669, -0.8678251267471937)
  ];
  let mut s = String::new();
  for coo in vertices.iter().map(|(a, b)| [*a, *b]).flatten() {
    s.push(',');
    s.push_str(&format!("{}", (coo as f64).to_degrees()));
  }
  println!("draw polygon({})", &s.as_str()[1..]);
  // In Aladin:
  // draw polygon(2.762765958009174, -1.1558461939389928, 1.8511822776529407, -1.0933557485463217, 1.2223938385258124, -0.9257482282529843, 0.7433555354994513, -0.7588156605217202,  0.3057973131669505, -0.6330807489082455, -0.13178353851370578, -0.5615526603403407,  -0.5856691165836956, -0.5400327446156468, -0.4015158856736461, -0.255385904573062,  -0.36447295312493533, 0.055767062790012194, -0.46118312909272624, 0.3439014598749883,  -0.7077120187325645, 0.5509825177593931, -1.0770630524112423, 0.6032198101231719,  -1.4301394011313524, 0.46923928877261123, -1.4878998343615297, 0.8530569092752573,  -1.9282441740007021, 1.1389659504677638, -2.7514745722416727, 1.0909408969626568,  -3.043766621461535, 0.7764645583461347, -3.0294959096529364, 0.4089232125719977,  -2.878305764522211, 0.05125802264111997, 3.0734739115502143, 0.057021963147195695, 2.7821939201916965, -0.05879878325187485, 2.5647762800185414, -0.27363879379223127, 2.4394176374317205, -0.5545388864809319, 2.444623710270669, -0.8678251267471937)
  // to_radians(&mut vertices);
  let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
  to_aladin_moc(&actual_res_exact);
  let hash = hash(depth, 2.762765958009174, -1.1558461939389928);
  println!("hash: {}", hash);
}*/

#[test]
fn testok_polygone_exact_markus() {
    // In Aladin:
    // draw polygon(...)
    let depth = 12;
    let coos = [
        1.8513602782977292_f64,
        -0.18592998414051692,
        1.8512258422891248,
        -0.18579783541887313,
        1.8511638051632155,
        -0.185675857268537,
        1.8511017577699869,
        -0.1856453564601204,
        1.8510397110839012,
        -0.1856148549652769,
        1.8509776667949591,
        -0.18557418856024166,
        1.8509052799363221,
        -0.18554368392668075,
        1.8508432374924337,
        -0.18550301603535885,
        1.8507708544432935,
        -0.1854623454440125,
        1.8507088117613895,
        -0.18543184028926155,
        1.8506467721028483,
        -0.18539117022619508,
        1.8505847358196907,
        -0.18534033525596078,
        1.8505330419248691,
        -0.18528950207341022,
        1.8504710184894149,
        -0.18519800896602515,
        1.850419329785655,
        -0.18513701052017761,
        1.8503883238565562,
        -0.18507601709824864,
        1.8503573186288944,
        -0.18501502350745627,
        1.8503263170225652,
        -0.1849438655291929,
        1.8502849728262747,
        -0.1848828686900894,
        1.8502643132238532,
        -0.184811713286104,
        1.8502539912630134,
        -0.18475072501341952,
        1.8502333294756885,
        -0.18468973371793343,
        1.850233354215554,
        -0.18460841998467828,
        1.8502333758556728,
        -0.18453727047097696,
        1.850243746325784,
        -0.18443563131867075,
        1.8502437739504332,
        -0.18434415338072863,
        1.850243801563952,
        -0.18425267544669602,
        1.850233483954607,
        -0.18418152293981005,
        1.8502335055540653,
        -0.18411037344041434,
        1.850223185237512,
        -0.18404938513277874,
        1.8501818522075997,
        -0.18397822335057393,
        1.8501405202635473,
        -0.18390706126594134,
        1.8500991860698157,
        -0.18384606309162968,
        1.8500475110573917,
        -0.1837952254787585,
        1.8499958370161058,
        -0.1837443873907134,
        1.8499441675707649,
        -0.18368338461546654,
        1.84990284099673,
        -0.18361222078765688,
        1.8498718533556298,
        -0.18354106034805112,
        1.8498305249747435,
        -0.18348006020220078,
        1.8498098766287663,
        -0.1834089031712837,
        1.8497995779550107,
        -0.18330725725741462,
        1.8497892640149673,
        -0.1832462681693668,
        1.8497892913795262,
        -0.18317511869547,
        1.8498099966304176,
        -0.18309381264003968,
        1.8498410342777623,
        -0.18302267446611675,
        1.8498927334119273,
        -0.18298203607383598,
        1.8499651015965932,
        -0.1829515684700408,
        1.850037472401915,
        -0.18291093572064188,
        1.8500995055507055,
        -0.18287029878505212,
        1.8501408712721863,
        -0.18280932641072814,
        1.8501718997230883,
        -0.1827483506134699,
        1.850182255036083,
        -0.18268736846694722,
        1.8501305967663695,
        -0.18261620321629635,
        1.8500582543810962,
        -0.18258568769748984,
        1.8499962452637886,
        -0.18256533891575163,
        1.8499238933743258,
        -0.18256531429383693,
        1.8498721909927536,
        -0.18262628139511278,
        1.8498308160127352,
        -0.1827075802042478,
        1.8497894475832763,
        -0.18276855028745648,
        1.8497377337273002,
        -0.18284984456644895,
        1.8496860346994894,
        -0.1828904815280938,
        1.8496239939937291,
        -0.1829412780938848,
        1.8495412828905007,
        -0.18298190121966892,
        1.8494792487450675,
        -0.1830123677623018,
        1.8494172231098334,
        -0.1830225051988025,
        1.8493448651057844,
        -0.1830224730958927,
        1.8492828585330476,
        -0.18299195220809547,
        1.8492828925303078,
        -0.18292080274479733,
        1.8492829216622935,
        -0.18285981749099187,
        1.8492829556396726,
        -0.1827886680286644,
        1.8492829944578515,
        -0.18270735435791027,
        1.8492830284122925,
        -0.18263620489631682,
        1.8492933935110967,
        -0.1825752244003239,
        1.849324430143685,
        -0.1825142533027064,
        1.8493347945288308,
        -0.182453272729678,
        1.849365839361594,
        -0.18237197298389143,
        1.849396873991587,
        -0.1823110014833154,
        1.8494382433000762,
        -0.18225003429912384,
        1.8495106086441788,
        -0.18220940835074173,
        1.849572633374496,
        -0.18217894145457372,
        1.849655331928133,
        -0.18213831785758774,
        1.8497276901408222,
        -0.18210785331326518,
        1.849800051409277,
        -0.1820672236239383,
        1.8498724078572921,
        -0.18203675720966933,
        1.8499344285723411,
        -0.1820062863060658,
        1.8500171252356468,
        -0.1819554931554334,
        1.850079147551529,
        -0.18191485643812733,
        1.8501825076068366,
        -0.18187423161174024,
        1.8502445277525188,
        -0.18183359306249244,
        1.8503168814810491,
        -0.18179295669421697,
        1.8503788941242136,
        -0.18177264507944257,
        1.8504512486801108,
        -0.18172184276349443,
        1.8505029333185445,
        -0.1816710345197491,
        1.8505442852757363,
        -0.1816100591568629,
        1.850585633907972,
        -0.1815592476998868,
        1.8506476565420458,
        -0.1814779478344996,
        1.8506786715997467,
        -0.18141696921824987,
        1.850709688119604,
        -0.18134582621623202,
        1.8507303681750746,
        -0.18128484516779256,
        1.8507510498462858,
        -0.18121369982844085,
        1.8507407311398347,
        -0.18114254826722087,
        1.8507304105793165,
        -0.18108156089920815,
        1.8506994263298109,
        -0.1810104050299821,
        1.8506477725552963,
        -0.18095957295895188,
        1.8505857869062527,
        -0.18090873805888288,
        1.850554803069368,
        -0.18084774560122668,
        1.8505858179298913,
        -0.1807766032700415,
        1.8506168296497307,
        -0.18071562497989227,
        1.850647840673332,
        -0.18065464651625232,
        1.8507201778206104,
        -0.180603840703177,
        1.8507821733848713,
        -0.18059368882855678,
        1.8508338422657715,
        -0.18055304172508826,
        1.8508648499351492,
        -0.1804920620550578,
        1.850864864843558,
        -0.1804107483219105,
        1.8508648797465577,
        -0.1803294345853517,
        1.850864892782244,
        -0.18025828506295855,
        1.8508752354417124,
        -0.18019730159328218,
        1.8508959092420871,
        -0.1801363199016534,
        1.8509269155834893,
        -0.18006517565476338,
        1.8509992407064624,
        -0.18002453044035827,
        1.8510612321725082,
        -0.1799940470354819,
        1.8511335526348274,
        -0.17997372852525145,
        1.851195539039006,
        -0.1799737362930964,
        1.851278187577929,
        -0.1799737455840484,
        1.851422822521087,
        -0.17997375891123377,
        1.8514951406828215,
        -0.17996359995564287,
        1.8515674580202341,
        -0.17996360428695535,
        1.8516707685022664,
        -0.179963608856252,
        1.851743085839695,
        -0.17996361092195487,
        1.8518050721289219,
        -0.17996361195004668,
        1.8518670584181505,
        -0.17995344807278466,
        1.8519393756213955,
        -0.17995344760633805,
        1.8520116928246395,
        -0.17995344620699838,
        1.8520840104304412,
        -0.1799636080947026,
        1.85217699158969,
        -0.17999409638466657,
        1.8522596433499237,
        -0.18003474808418177,
        1.8523216338197017,
        -0.18007540027750576,
        1.8523836232924833,
        -0.1800957233463537,
        1.8524559467486013,
        -0.18013637309113859,
        1.8525282700518246,
        -0.18016685768320762,
        1.8525902599855844,
        -0.1801770142472329,
        1.852662584575065,
        -0.1802074971054239,
        1.8527349067467758,
        -0.18021765059379727,
        1.8528175607341215,
        -0.18022780140668926,
        1.8529105465946833,
        -0.18023794908706292,
        1.8529311978668475,
        -0.18017695989513224,
        1.8529001920025943,
        -0.1801159803832753,
        1.8528485225687163,
        -0.1800448401341806,
        1.8527865285329002,
        -0.1800041937674299,
        1.8527245370097851,
        -0.17997371093336756,
        1.8526625447015321,
        -0.17993306319323382,
        1.8525902237418577,
        -0.17990258032858317,
        1.852528235689744,
        -0.17988225954154966,
        1.8525902170349364,
        -0.1798517592273769,
        1.8526522020591731,
        -0.17985175088889574,
        1.852745179595501,
        -0.17985173709613428,
        1.8528071646197006,
        -0.17985172704426694,
        1.8528691533611787,
        -0.17987204474727547,
        1.8529414794564334,
        -0.17992285245321354,
        1.8530034737194898,
        -0.179963497108399,
        1.8530757955426747,
        -0.17998381042047257,
        1.853137789136592,
        -0.1800142893692613,
        1.8532307820893295,
        -0.18006508861459566,
        1.8532927774692713,
        -0.1800955658471284,
        1.8533444476726812,
        -0.180146373561122,
        1.8533857958125843,
        -0.1802276762574201,
        1.8534271393909703,
        -0.18028865021083082,
        1.8534684809191593,
        -0.18033945964033768,
        1.8534994974343275,
        -0.18041060021371805,
        1.8535305147629872,
        -0.1804817406124282,
        1.8535718618956585,
        -0.18054271348844708,
        1.8535925454408997,
        -0.1806036924563041,
        1.8535615675880976,
        -0.180664687201364,
        1.8534995824661735,
        -0.180695198243738,
        1.8534375908057812,
        -0.18070538017125248,
        1.853365263656628,
        -0.1807053999696568,
        1.8532929365074355,
        -0.18070541883480354,
        1.8532309443469994,
        -0.1807155984771237,
        1.8531689592273146,
        -0.180756270078698,
        1.853168976197942,
        -0.18082741958169457,
        1.8531999916204283,
        -0.18089856179674516,
        1.8532516711771596,
        -0.18095953455695846,
        1.8533240071244552,
        -0.18097984465025554,
        1.8533963464617944,
        -0.18101031802316253,
        1.8534583475915447,
        -0.18102046497809307,
        1.8535203520344468,
        -0.1810407754602831,
        1.8535926932809386,
        -0.18107124629742702,
        1.853551379553357,
        -0.1811322441480593,
        1.853468719459801,
        -0.18115259680598478,
        1.8533860559168442,
        -0.1811627840320468,
        1.853313725401153,
        -0.18117296738149677,
        1.8532413971792996,
        -0.1811933140107519,
        1.8531897401786994,
        -0.18123398329392268,
        1.853158754286251,
        -0.18129497580183623,
        1.8531587711626412,
        -0.18136612529137458,
        1.8531794557183354,
        -0.18143726997806436,
        1.853210472282553,
        -0.1814982479071053,
        1.8532724917514882,
        -0.18155921797476143,
        1.8533345126147671,
        -0.18162018735461194,
        1.8533758634520165,
        -0.18167099743739729,
        1.853458540674906,
        -0.18168113878211153,
        1.8535308838873166,
        -0.1816912819845117,
        1.8536135584443616,
        -0.18169125683075096,
        1.8536755677428784,
        -0.18170140137622234,
        1.8537375632869104,
        -0.18167088839128523,
        1.8537995509138552,
        -0.1816200462987067,
        1.8538615448392612,
        -0.18158953194282895,
        1.8539028623273504,
        -0.1815386960249918,
        1.8539545090822382,
        -0.18147769175202094,
        1.853995824861937,
        -0.18142685514865847,
        1.8540578073165581,
        -0.18137601020034316,
        1.85410945917895,
        -0.18133533292176704,
        1.85418177740647,
        -0.18129464667012765,
        1.8542541168289186,
        -0.18130478054090335,
        1.8542851446671542,
        -0.18136575251690873,
        1.854274842432865,
        -0.18143690644015867,
        1.8542438675549624,
        -0.18149790493589874,
        1.854212891989415,
        -0.1815589032600803,
        1.8541819157359403,
        -0.18161990141274395,
        1.854140604494005,
        -0.18168090359572794,
        1.854088957915552,
        -0.18174190958079806,
        1.854026967573646,
        -0.18178259066192382,
        1.8539649723831393,
        -0.1818131068473544,
        1.853902976498352,
        -0.1818436223473159,
        1.8538409873050585,
        -0.1818944655828114,
        1.853789335355931,
        -0.18195546880818889,
        1.8537583521286545,
        -0.1820164646202317,
        1.8537583733736718,
        -0.1820774498829631,
        1.8537583981667958,
        -0.18214859935568553,
        1.8537274204772618,
        -0.18222992341714173,
        1.8536860961411439,
        -0.18228075804079932,
        1.8536344321829328,
        -0.18232143141876903,
        1.8535724222652281,
        -0.18233161484398933,
        1.8535000734963125,
        -0.18233163639387198,
        1.8534483839488172,
        -0.1822909943717122,
        1.8534173601355401,
        -0.18223001777030146,
        1.8534173427097982,
        -0.1821690325050502,
        1.8534069901763455,
        -0.18210805008898873,
        1.8533759678871828,
        -0.182047073256815,
        1.8532209657591807,
        -0.1821385907538852,
        1.8532209657591807,
        -0.1821385907538852,
        1.8531589595755087,
        -0.18215893381474696,
        1.853117630451771,
        -0.18220976424895258,
        1.8530659629067465,
        -0.18225043238787383,
        1.8530246322191117,
        -0.18230126213656264,
        1.8529833049442839,
        -0.18237242000328946,
        1.8529316368446482,
        -0.1824232511157904,
        1.8528799677781902,
        -0.18247408175220095,
        1.8528282977445811,
        -0.182524911912516,
        1.8527456170839731,
        -0.18256558233245482,
        1.852693944920749,
        -0.1826164112546596,
        1.8526215951840501,
        -0.18263674988737813,
        1.8525492436239446,
        -0.18264692337438293,
        1.8524768883584832,
        -0.18262660329194685,
        1.8523838605427814,
        -0.18259612000798175,
        1.8522908361603814,
        -0.18259612781505855,
        1.8522184838629452,
        -0.18259613281959508,
        1.852104788287783,
        -0.18261646722026986,
        1.8520531092781092,
        -0.18265712602218745,
        1.8520427750984183,
        -0.18271811162760823,
        1.8521151290255806,
        -0.1827181088916202,
        1.8521771466774282,
        -0.18271810580290904,
        1.8522495006045756,
        -0.1827181013319043,
        1.8523218562395611,
        -0.18273842435046514,
        1.8523632070863645,
        -0.18279940611373796,
        1.8524355642360182,
        -0.18281972766392637,
        1.8524975830524693,
        -0.18281972102856434,
        1.8525906112771218,
        -0.18281970978849002,
        1.8526526286180132,
        -0.18280953722531323,
        1.8527146457259827,
        -0.18279936397582477,
        1.8527869955569074,
        -0.182768859929351,
        1.8528490119468843,
        -0.1827586851930074,
        1.8529110241835405,
        -0.18272818134741206,
        1.8529730357223029,
        -0.1826976768158609,
        1.8530350465633145,
        -0.18266717159837986,
        1.8530970567067129,
        -0.18263666569499454,
        1.8531590661526396,
        -0.18260615910573078,
        1.8531900594433808,
        -0.1825451666039551,
        1.8532003803817598,
        -0.18248417888771545,
        1.8531693557982265,
        -0.18241303670180406,
        1.853221024111172,
        -0.18237236760955008,
        1.8532933734227546,
        -0.18237234966205038,
        1.8533553899086845,
        -0.18238249774591547,
        1.8534174124410632,
        -0.1824129735649136,
        1.8534484373148825,
        -0.1824739501651511,
        1.853448458080047,
        -0.18254509964029003,
        1.8534484758839718,
        -0.18260608490476157,
        1.8535001684529953,
        -0.18264672692445688,
        1.8535621948350693,
        -0.1826772011390782,
        1.8536241988278674,
        -0.18263652519325926,
        1.853655187053907,
        -0.18257553012162162,
        1.8536655028126596,
        -0.18251454155155494,
        1.8536758183442161,
        -0.18245355296267193,
        1.8537171406126367,
        -0.18239255420394243,
        1.8537998077329572,
        -0.1823417052451024,
        1.8538514706437519,
        -0.18230103034319703,
        1.8538721153030797,
        -0.18222987351480946,
        1.853872092765795,
        -0.18216888825415814,
        1.8538720664800044,
        -0.18209773878303329,
        1.8539237226908758,
        -0.1820468990060636,
        1.8539857287020418,
        -0.18203671169602123,
        1.8540580694387256,
        -0.1820265196690583,
        1.854140749099136,
        -0.18202648673473754,
        1.8542130938019286,
        -0.1820264569165996,
        1.854337108668429,
        -0.18201623941775788,
        1.8544094484765554,
        -0.1820060428557937,
        1.8544507543028441,
        -0.181934874485717,
        1.8544610551126304,
        -0.18186372024756248,
        1.85449202468799,
        -0.18179255634714783,
        1.8545229984280474,
        -0.18173155648472142,
        1.8545746400423422,
        -0.18167054650397377,
        1.854626285662072,
        -0.18161970025692645,
        1.8546882749522515,
        -0.18158917676023903,
        1.8547295847169856,
        -0.1815383347544417,
        1.8547398865501494,
        -0.1814773442116756,
        1.854750182777349,
        -0.18140618943984393,
        1.8547708123902436,
        -0.18133502930381967,
        1.8547811134037422,
        -0.18127403868276232,
        1.854812081019679,
        -0.18121303722084287,
        1.8548843980213114,
        -0.18118250611807843,
        1.8549567315037854,
        -0.18118246671147994,
        1.8549980941266584,
        -0.18123326482089308,
        1.8550394575272688,
        -0.18128406262345753,
        1.8551118226373557,
        -0.1813348422576304,
        1.855142872491283,
        -0.18141613789159422,
        1.8551635831794189,
        -0.18148727523052058,
        1.8551739541751902,
        -0.18154825439023753,
        1.855173997427532,
        -0.18161940384446773,
        1.8551740345118601,
        -0.1816803890900812,
        1.8551844122249557,
        -0.18175153243476497,
        1.8551844680618457,
        -0.18184301030035063,
        1.855194846250852,
        -0.18191415362325747,
        1.8552052247243542,
        -0.18198529692593465,
        1.85521559721333,
        -0.18204627600163992,
        1.8552259699446727,
        -0.18210725505748945,
        1.8552570133751127,
        -0.18216822162789617,
        1.8552983928572575,
        -0.18222918170647745,
        1.8553604311072789,
        -0.18226980021820402,
        1.8554328057910092,
        -0.1823104114751146,
        1.8555051679248793,
        -0.18233069338463606,
        1.8555878661791823,
        -0.18235096750615315,
        1.8556498867316906,
        -0.18236109019332603,
        1.8557325932425572,
        -0.1823915263827604,
        1.855794629261471,
        -0.18242197587644818,
        1.8558773224810987,
        -0.18243208151878457,
        1.855949672603792,
        -0.18243202927598628,
        1.855980726089584,
        -0.18249299182360623,
        1.8559394208942452,
        -0.18254384304099072,
        1.8558877788897787,
        -0.1825947013807476,
        1.8558361359239577,
        -0.18264555924593065,
        1.8557741484255865,
        -0.1826862596482374,
        1.8557018381076122,
        -0.1827472948592403,
        1.8556398484486467,
        -0.1827879937788979,
        1.8555675213744636,
        -0.18282869884964942,
        1.8554951932178876,
        -0.1828694029882902,
        1.8554332140249616,
        -0.18293042803810927,
        1.8553919136685315,
        -0.18300160364397636,
        1.855381623102711,
        -0.18307275958212196,
        1.8553816759809545,
        -0.18315407323467503,
        1.8553920595113906,
        -0.18322521618783746,
        1.855392099311757,
        -0.18328620142857147,
        1.8553818082597895,
        -0.18335735737122666,
        1.8553508416788267,
        -0.18342852618988437,
        1.8553302055429384,
        -0.18348952425176293,
        1.8553095689502415,
        -0.18355052223910293,
        1.8552682623384746,
        -0.18362169694954694,
        1.8552372864379694,
        -0.1836827009405365,
        1.8552166480855876,
        -0.18374369859095774,
        1.8552166986250527,
        -0.18382501226104841,
        1.8552167428619664,
        -0.18389616172419765,
        1.855206442141301,
        -0.18395715315164085,
        1.8552064862673985,
        -0.18402830261842207,
        1.8552168693296875,
        -0.1840994459148106,
        1.8552375853711334,
        -0.18416041877281975,
        1.855237629947007,
        -0.18423156824474607,
        1.855237680907982,
        -0.18431288192966777,
        1.8552377318872193,
        -0.18439419561762005,
        1.8552274305273724,
        -0.18445518709732775,
        1.8552171352810847,
        -0.18452634277208088,
        1.8552068397748254,
        -0.18459749843080292,
        1.855206884025601,
        -0.1846686479184675,
        1.8551965818312721,
        -0.1847296393511224,
        1.8551759642309624,
        -0.18483129373905474,
        1.8551760081158577,
        -0.1849024432368502,
        1.8551657051572494,
        -0.18496343462191608,
        1.8551347639410862,
        -0.18509558760038614,
        1.8551244601711432,
        -0.18515657891799897,
        1.8550728040399498,
        -0.18523792247843815,
        1.8550418228732382,
        -0.1853090896669299,
        1.8549901581162065,
        -0.1853802682587192,
        1.854948821818579,
        -0.18543111226157105,
        1.854876454267602,
        -0.18547180849782238,
        1.8548247857687108,
        -0.18554298557650512,
        1.8547420794541367,
        -0.18559384975419665,
        1.8546593662655284,
        -0.18563454849317843,
        1.8545662944694563,
        -0.1856447583426744,
        1.8545042480799074,
        -0.18565495212383804,
        1.8544318595148424,
        -0.18566514996322853,
        1.854359475413111,
        -0.1856855110893831,
        1.8542767670714466,
        -0.18574653263011986,
        1.8542250834943699,
        -0.18580753997754113,
        1.8541837366695544,
        -0.18585837835974803,
        1.8541527401739855,
        -0.18592954067736409,
        1.8541217385690594,
        -0.18599053860317263,
        1.8540803893603461,
        -0.18604137623222627,
        1.8540287089531733,
        -0.1861227102396788,
        1.8539873578459334,
        -0.18617354718929408,
        1.8539460218020507,
        -0.18626504074818293,
        1.853915015919087,
        -0.1863260375513494,
        1.8538736698208658,
        -0.18639720213192207,
        1.8538219791307038,
        -0.18646837003308542,
        1.8537806343526355,
        -0.18654969816950612,
        1.8537392847482848,
        -0.1866208617760621,
        1.8536772461860471,
        -0.18669203180853772,
        1.8536048585853184,
        -0.18675304014573832,
        1.8535324629682506,
        -0.18679371908332434,
        1.8534600662469602,
        -0.18683439708737642,
        1.8533876568005647,
        -0.18683441721708252,
        1.8533566071093006,
        -0.1867734401494385,
        1.8533565871984605,
        -0.18670229051176154,
        1.8533565672938108,
        -0.18663114087943428,
        1.8533668824835334,
        -0.18652949579605713,
        1.8533772059318696,
        -0.18645834339686834,
        1.8533771828928605,
        -0.18637702955217636,
        1.853387505920691,
        -0.18630587714490543,
        1.8533978286788215,
        -0.18623472472359145,
        1.85339781117922,
        -0.18617373935273773,
        1.8533977849385794,
        -0.18608226130288735,
        1.8533977645364603,
        -0.18601111171378606,
        1.8533977354016076,
        -0.18590946945138145,
        1.8533977150149308,
        -0.18583831987299143,
        1.8533976946346133,
        -0.18576717029885748,
        1.8533976742606526,
        -0.18569602072888974,
        1.8533562921093778,
        -0.18564521076173426,
        1.853294235280328,
        -0.185624898469215,
        1.8532321815180122,
        -0.1856147497127729,
        1.8531701304682395,
        -0.18561476449242317,
        1.8531184331532222,
        -0.1856655974032131,
        1.8530563905581923,
        -0.18570626782034258,
        1.85298400067569,
        -0.18572661117418918,
        1.8529116062702775,
        -0.18572662514312013,
        1.852828869806907,
        -0.18572663996183605,
        1.8527254475954322,
        -0.18571649254160064,
        1.8526427083397266,
        -0.18569617616057874,
        1.8525806522761148,
        -0.18566569173581843,
        1.8524875724573289,
        -0.18564537437166742,
        1.8524255210534533,
        -0.18564538090231686,
        1.8523427858482646,
        -0.18564538854050341,
        1.8522600498960486,
        -0.18563523073199883,
        1.852177314258372,
        -0.1856250717014105,
        1.852084235855293,
        -0.18558441918638036,
        1.8520118425508667,
        -0.18555392885326194,
        1.8519497911080103,
        -0.18548278054767153,
        1.851908424134568,
        -0.1854014672210403,
        1.8518773997294453,
        -0.1853404820287538,
        1.8518463760311155,
        -0.18527949666729704,
        1.8518050118457905,
        -0.1852286752531498,
        1.8517739901487216,
        -0.18513719683245522,
        1.8517326282580837,
        -0.1850660464477075,
        1.8516705858039069,
        -0.1850050592936848,
        1.851598203195129,
        -0.18496439940759357,
        1.851525821036894,
        -0.18493390280685654,
        1.8514637786409567,
        -0.18492373446428506,
        1.8514017364801523,
        -0.18491356543464932,
        1.8513396935551185,
        -0.18491355993727504,
        1.851277648396113,
        -0.18493388219148216,
        1.851225939686831,
        -0.18498469761055747,
        1.8511845707188266,
        -0.18503551382315714,
        1.8511431968466974,
        -0.18511682239303032,
        1.8511328476108682,
        -0.18517780637161896,
        1.8511431803673817,
        -0.18523879304451857,
        1.851163854589929,
        -0.18529978100878403,
        1.8511948706319659,
        -0.185360770150041,
        1.8512362276491314,
        -0.18543192451942578,
        1.8512775880031034,
        -0.1854827501407481,
        1.85132929184779,
        -0.18552341223756463,
        1.851360310101894,
        -0.18559456469755148,
        1.8513913300837215,
        -0.18565555276513232,
        1.8514326928101248,
        -0.18571654147556824,
        1.851443030010248,
        -0.18577752762071573,
        1.8514326828964742,
        -0.18583851218390054,
        1.8513602782977292,
        -0.18592998414051692,
    ];
    let mut vertices: Vec<(f64, f64)> = Vec::with_capacity(coos.len() >> 1);
    for i in (0..coos.len()).step_by(2) {
        vertices.push((coos[i], coos[i + 1]))
    }
    // let expected_res_exact: [u64; 9] =...
    // Problem was an error in Newtom-Raphson method
    let _actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    // So far we are happy if it does not panic :)
    // to_aladin_moc(&actual_res_exact);
}

#[test]
#[allow(clippy::excessive_precision)]
fn testok_polygone_exact_mt_pb_with_java_lib() {
    // In Aladin:
    // draw polygon(...)
    let depth = 12;
    // First polygone
    let coos = [
        -1.5702949547333407,
        -0.7295093151415473,
        -1.5702171673769950,
        -0.7295093016804524,
        -1.5701393800214274,
        -0.7295092852142693,
        -1.5700615926667945,
        -0.7295092657429985,
    ];
    let mut vertices: Vec<(f64, f64)> = Vec::with_capacity(coos.len() >> 1);
    for i in (0..coos.len()).step_by(2) {
        vertices.push((coos[i], coos[i + 1]))
    }
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    to_aladin_moc(&actual_res_exact);
    let mut s = String::new();
    for coo in vertices.iter().flat_map(|(a, b)| [*a, *b]) {
        s.push(',');
        s.push_str(&coo.to_degrees().to_string());
    }
    println!("draw polygon({})", &s.as_str()[1..]);

    // Second polygon
    let coos = [
        -1.5706045044233712,
        -0.7295105218483977,
        -1.5705168372776197,
        -0.7295105199399403,
        -1.5704291701320918,
        -0.7295105142145686,
        -1.5703415029870114,
        -0.7295105046722821,
    ];
    let mut vertices: Vec<(f64, f64)> = Vec::with_capacity(coos.len() >> 1);
    for i in (0..coos.len()).step_by(2) {
        vertices.push((coos[i], coos[i + 1]))
    }
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    to_aladin_moc(&actual_res_exact);
    let mut s = String::new();
    for coo in vertices.iter().flat_map(|(a, b)| [*a, *b]) {
        s.push(',');
        s.push_str(&coo.to_degrees().to_string());
    }
    println!("draw polygon({})", &s.as_str()[1..]);
}

#[test]
fn testok_polygone_exact_7() {
    // In Aladin: draw polygon(ra1, dec1, ra2, dec2, ...)
    let depth = 4;
    let vertices = [(PI, 1.000), (PI + 1e-3, 1.001), (PI - 1e-3, 1.002)];
    let expected_res_exact: [u64; 3] = [373, 375, 698];
    let actual_res_exact = get(depth).polygon_coverage(&vertices, true);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
        assert_eq!(h1, *h2);
    }
}

#[test]
fn testok_zone_1() {
    let depth = 5;
    let (lon_min, lat_min, lon_max, lat_max) = (
        0.0_f64.to_radians(),
        0.0_f64.to_radians(),
        10.0_f64.to_radians(),
        10.0_f64.to_radians(),
    );
    let expected_res_exact: [u64; 44] = [
        4517, 4518, 4519, 4521, 4522, 4523, 4524, 4525, 4526, 4527, 4528, 4530, 4531, 4536, 4537,
        4538, 4539, 4540, 4542, 4543, 4864, 4865, 4867, 4868, 4869, 4870, 4871, 4876, 4877, 4879,
        4880, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4912,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    // assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_2() {
    let depth = 2;
    let (lon_min, lat_min, lon_max, lat_max) = (
        284.2_f64.to_radians(),
        -25.0_f64.to_radians(),
        26.7_f64.to_radians(),
        85.0_f64.to_radians(),
    );
    let expected_res_exact: [u64; 53] = [
        2, 8, 9, 10, 11, 14, 15, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
        64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 113, 116, 117, 118, 119,
        125, 139, 142, 183, 187, 188, 189, 190, 191,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_3() {
    let depth = 4;
    let (lon_min, lat_min, lon_max, lat_max) = (
        12.7_f64.to_radians(),
        64.5_f64.to_radians(),
        350.2_f64.to_radians(),
        65.0_f64.to_radians(),
    );
    // println!("Test zone (12.7, 64.5) -> (350.2, 65.0)");
    let expected_res_exact: [u64; 65] = [
        127, 207, 211, 212, 213, 214, 216, 217, 218, 227, 228, 229, 230, 232, 233, 383, 447, 463,
        467, 468, 469, 470, 472, 473, 474, 483, 484, 485, 486, 488, 489, 490, 639, 703, 719, 723,
        724, 725, 726, 728, 729, 730, 739, 740, 741, 742, 744, 745, 746, 959, 975, 979, 980, 981,
        982, 984, 985, 986, 995, 996, 997, 998, 1000, 1001, 1002,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_4() {
    let depth = 5;
    let (lon_min, lat_min, lon_max, lat_max) = (
        75.7_f64.to_radians(),
        -75.0_f64.to_radians(),
        12.2_f64.to_radians(),
        -74.5_f64.to_radians(),
    );
    let expected_res_exact: [u64; 71] = [
        8257, 8258, 8259, 8260, 8321, 8322, 8323, 8328, 9245, 9246, 9247, 9261, 9262, 9263, 9265,
        9266, 9267, 9268, 9272, 9281, 9282, 9283, 9284, 9288, 9345, 9346, 9347, 9348, 9352, 10269,
        10270, 10271, 10285, 10286, 10287, 10289, 10290, 10291, 10292, 10296, 10305, 10306, 10307,
        10308, 10312, 10369, 10370, 10371, 10372, 10376, 11293, 11294, 11295, 11309, 11310, 11311,
        11313, 11314, 11315, 11316, 11320, 11329, 11330, 11331, 11332, 11336, 11393, 11394, 11395,
        11396, 11400,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_5() {
    let depth = 3;
    let (lon_min, lat_min, lon_max, lat_max) = (
        180_f64.to_radians(),
        (30_f64 + 1e-14_f64).to_radians(),
        360_f64.to_radians(),
        50_f64.to_radians(),
    );
    let expected_res_exact: [u64; 75] = [
        141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 161,
        162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 176, 177, 178, 205, 206,
        207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 225, 226, 227,
        228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 240, 241, 242, 318, 319, 445, 447,
        509, 510, 511,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    // assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_6() {
    let depth = 3;
    let (lon_min, lat_min, lon_max, lat_max) = (
        0_f64.to_radians(),
        -90.0_f64.to_radians(),
        180.0_f64.to_radians(),
        0_f64.to_radians(),
    );

    let expected_res_exact: [u64; 204] = [
        256, 257, 259, 260, 261, 262, 263, 268, 269, 271, 272, 273, 274, 275, 276, 277, 278, 280,
        281, 282, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335,
        336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 352, 353, 354, 355, 356, 357, 358, 360,
        361, 362, 384, 386, 387, 392, 393, 394, 395, 396, 398, 399, 416, 417, 418, 419, 420, 421,
        422, 424, 425, 426, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525,
        526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543,
        544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561,
        562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579,
        580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597,
        598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615,
        616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633,
        634, 635, 636, 637, 638, 639,
    ];

    // println!("Zone: (0, -90) -> (180, 0)");
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_7() {
    let depth = 5;
    let (lon_min, lat_min, lon_max, lat_max) = (
        20_f64.to_radians(),
        0_f64.to_radians(),
        30.0_f64.to_radians(),
        89_f64.to_radians(),
    );

    let expected_res_exact: [u64; 228] = [
        136, 138, 139, 160, 161, 162, 163, 164, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        184, 186, 187, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 523, 524, 525, 526, 527,
        528, 529, 530, 531, 532, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 548, 549, 551,
        560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 571, 572, 573, 574, 575, 584, 586, 587,
        608, 609, 610, 611, 612, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 632, 633, 634,
        635, 638, 639, 660, 661, 663, 704, 705, 706, 707, 708, 709, 710, 711, 717, 720, 721, 722,
        723, 724, 725, 726, 727, 728, 729, 732, 733, 734, 735, 896, 897, 898, 899, 902, 903, 904,
        905, 906, 907, 908, 909, 910, 911, 920, 921, 922, 923, 926, 927, 932, 933, 944, 945, 947,
        948, 949, 950, 951, 992, 993, 994, 995, 998, 999, 1001, 1004, 1005, 1016, 1017, 1018, 1019,
        1022, 1023, 4454, 4455, 4457, 4458, 4459, 4460, 4461, 4462, 4463, 4472, 4474, 4475, 4501,
        4503, 4544, 4545, 4546, 4547, 4548, 4549, 4550, 4551, 4552, 4553, 4555, 4556, 4557, 4558,
        4559, 4560, 4561, 4562, 4563, 4564, 4566, 4567, 4568, 4569, 4570, 4571, 4572, 4573, 4574,
        4575, 4580, 4581, 4583, 4592, 4593, 4594, 4595, 4596, 4597, 4598, 4599, 4600, 4601, 4603,
        4604, 4605, 4606, 4607, 4948, 4949, 4951,
    ];

    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);

    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_8() {
    let depth = 5;
    let (lon_min, lat_min, lon_max, lat_max) = (
        30_f64.to_radians(),
        63_f64.to_radians(),
        20.0_f64.to_radians(),
        89_f64.to_radians(),
    );

    let expected_res_exact: [u64; 755] = [
        503, 507, 508, 509, 510, 511, 759, 763, 764, 765, 766, 767, 823, 827, 828, 829, 830, 831,
        839, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859,
        860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877,
        878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895,
        903, 907, 909, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 923, 924, 925, 926, 927,
        928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945,
        946, 947, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964,
        965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982,
        983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000,
        1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015,
        1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1527, 1531, 1532, 1533, 1534, 1535, 1783,
        1787, 1788, 1789, 1790, 1791, 1847, 1851, 1852, 1853, 1854, 1855, 1863, 1867, 1868, 1869,
        1870, 1871, 1872, 1873, 1874, 1875, 1876, 1877, 1878, 1879, 1880, 1881, 1882, 1883, 1884,
        1885, 1886, 1887, 1888, 1889, 1890, 1891, 1892, 1893, 1894, 1895, 1896, 1897, 1898, 1899,
        1900, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 1911, 1912, 1913, 1914,
        1915, 1916, 1917, 1918, 1919, 1927, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939,
        1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1950, 1951, 1952, 1953, 1954,
        1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969,
        1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984,
        1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999,
        2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014,
        2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029,
        2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 2043, 2044,
        2045, 2046, 2047, 2551, 2555, 2556, 2557, 2558, 2559, 2807, 2811, 2812, 2813, 2814, 2815,
        2871, 2875, 2876, 2877, 2878, 2879, 2887, 2891, 2892, 2893, 2894, 2895, 2896, 2897, 2898,
        2899, 2900, 2901, 2902, 2903, 2904, 2905, 2906, 2907, 2908, 2909, 2910, 2911, 2912, 2913,
        2914, 2915, 2916, 2917, 2918, 2919, 2920, 2921, 2922, 2923, 2924, 2925, 2926, 2927, 2928,
        2929, 2930, 2931, 2932, 2933, 2934, 2935, 2936, 2937, 2938, 2939, 2940, 2941, 2942, 2943,
        2951, 2955, 2956, 2957, 2958, 2959, 2960, 2961, 2962, 2963, 2964, 2965, 2966, 2967, 2968,
        2969, 2970, 2971, 2972, 2973, 2974, 2975, 2976, 2977, 2978, 2979, 2980, 2981, 2982, 2983,
        2984, 2985, 2986, 2987, 2988, 2989, 2990, 2991, 2992, 2993, 2994, 2995, 2996, 2997, 2998,
        2999, 3000, 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013,
        3014, 3015, 3016, 3017, 3018, 3019, 3020, 3021, 3022, 3023, 3024, 3025, 3026, 3027, 3028,
        3029, 3030, 3031, 3032, 3033, 3034, 3035, 3036, 3037, 3038, 3039, 3040, 3041, 3042, 3043,
        3044, 3045, 3046, 3047, 3048, 3049, 3050, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058,
        3059, 3060, 3061, 3062, 3063, 3064, 3065, 3066, 3067, 3068, 3069, 3070, 3071, 3575, 3579,
        3580, 3581, 3582, 3583, 3831, 3835, 3836, 3837, 3838, 3839, 3895, 3899, 3900, 3901, 3902,
        3903, 3911, 3915, 3916, 3917, 3918, 3919, 3920, 3921, 3922, 3923, 3924, 3925, 3926, 3927,
        3928, 3929, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 3939, 3940, 3941, 3942,
        3943, 3944, 3945, 3946, 3947, 3948, 3949, 3950, 3951, 3952, 3953, 3954, 3955, 3956, 3957,
        3958, 3959, 3960, 3961, 3962, 3963, 3964, 3965, 3966, 3967, 3975, 3979, 3980, 3981, 3982,
        3983, 3984, 3985, 3986, 3987, 3988, 3989, 3990, 3991, 3992, 3993, 3994, 3995, 3996, 3997,
        3998, 3999, 4000, 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009, 4010, 4011, 4012,
        4013, 4014, 4015, 4016, 4017, 4018, 4019, 4020, 4021, 4022, 4023, 4024, 4025, 4026, 4027,
        4028, 4029, 4030, 4031, 4032, 4033, 4034, 4035, 4036, 4037, 4038, 4039, 4040, 4041, 4042,
        4043, 4044, 4045, 4046, 4047, 4048, 4049, 4050, 4051, 4052, 4053, 4054, 4055, 4056, 4057,
        4058, 4059, 4060, 4061, 4062, 4063, 4064, 4065, 4066, 4067, 4068, 4069, 4070, 4071, 4072,
        4073, 4074, 4075, 4076, 4077, 4078, 4079, 4080, 4081, 4082, 4083, 4084, 4085, 4086, 4087,
        4088, 4089, 4090, 4091, 4092, 4093, 4094, 4095,
    ];

    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);

    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_9() {
    let depth = 5;

    let (lon_min, lat_min, lon_max, lat_max) = (
        355_f64.to_radians(),
        -90_f64.to_radians(),
        5.0_f64.to_radians(),
        90_f64.to_radians(),
    );

    let expected_res_exact: [u64; 468] = [
        672, 674, 675, 680, 681, 682, 683, 684, 685, 686, 687, 696, 697, 698, 699, 700, 701, 702,
        703, 744, 745, 746, 747, 748, 749, 750, 751, 760, 761, 762, 763, 764, 766, 767, 938, 939,
        942, 943, 954, 955, 958, 959, 1002, 1003, 1006, 1007, 1018, 1019, 1022, 1023, 3408, 3409,
        3411, 3412, 3413, 3414, 3415, 3420, 3421, 3422, 3423, 3444, 3445, 3446, 3447, 3452, 3453,
        3454, 3455, 3540, 3541, 3542, 3543, 3548, 3549, 3550, 3551, 3572, 3573, 3574, 3575, 3580,
        3581, 3583, 3925, 3927, 3933, 3935, 3957, 3959, 3965, 3967, 4053, 4055, 4061, 4063, 4085,
        4087, 4093, 4095, 4096, 4097, 4098, 4099, 4100, 4101, 4102, 4103, 4104, 4105, 4106, 4107,
        4108, 4109, 4110, 4111, 4112, 4114, 4115, 4120, 4121, 4122, 4123, 4124, 4126, 4127, 4128,
        4129, 4131, 4132, 4133, 4134, 4135, 4140, 4141, 4143, 4144, 4145, 4146, 4147, 4148, 4149,
        4150, 4151, 4152, 4153, 4154, 4155, 4156, 4157, 4158, 4159, 4192, 4194, 4195, 4200, 4201,
        4202, 4203, 4204, 4206, 4207, 4240, 4241, 4243, 4244, 4245, 4246, 4247, 4252, 4253, 4255,
        4288, 4289, 4290, 4291, 4292, 4293, 4294, 4295, 4296, 4297, 4298, 4299, 4300, 4301, 4302,
        4303, 4304, 4306, 4307, 4312, 4313, 4314, 4315, 4316, 4318, 4319, 4320, 4321, 4323, 4324,
        4325, 4326, 4327, 4332, 4333, 4335, 4336, 4337, 4338, 4339, 4340, 4341, 4342, 4343, 4344,
        4345, 4346, 4347, 4348, 4349, 4350, 4351, 4512, 4514, 4515, 4520, 4521, 4522, 4523, 4524,
        4526, 4527, 4688, 4689, 4691, 4692, 4693, 4694, 4695, 4700, 4701, 4703, 4864, 4865, 4866,
        4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4882,
        4883, 4888, 4889, 4890, 4891, 4892, 4894, 4895, 4896, 4897, 4899, 4900, 4901, 4902, 4903,
        4908, 4909, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923,
        4924, 4925, 4926, 4927, 4960, 4962, 4963, 4968, 4969, 4970, 4971, 4972, 4974, 4975, 5008,
        5009, 5011, 5012, 5013, 5014, 5015, 5020, 5021, 5023, 5056, 5057, 5058, 5059, 5060, 5061,
        5062, 5063, 5064, 5065, 5066, 5067, 5068, 5069, 5070, 5071, 5072, 5074, 5075, 5080, 5081,
        5082, 5083, 5084, 5086, 5087, 5088, 5089, 5091, 5092, 5093, 5094, 5095, 5100, 5101, 5103,
        5104, 5105, 5106, 5107, 5108, 5109, 5110, 5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118,
        5119, 8192, 8194, 8200, 8202, 8224, 8226, 8232, 8234, 8320, 8322, 8328, 8330, 8352, 8354,
        8360, 8362, 8704, 8706, 8707, 8712, 8713, 8714, 8715, 8736, 8737, 8738, 8739, 8744, 8745,
        8746, 8747, 8832, 8833, 8834, 8835, 8840, 8841, 8842, 8843, 8864, 8865, 8866, 8867, 8872,
        8873, 8874, 8875, 8876, 8878, 8879, 11264, 11265, 11268, 11269, 11280, 11281, 11284, 11285,
        11328, 11329, 11332, 11333, 11344, 11345, 11348, 11349, 11520, 11521, 11523, 11524, 11525,
        11526, 11527, 11536, 11537, 11538, 11539, 11540, 11541, 11542, 11543, 11584, 11585, 11586,
        11587, 11588, 11589, 11590, 11591, 11600, 11601, 11602, 11603, 11604, 11605, 11606, 11607,
        11612, 11613, 11615,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);

    // to_aladin_moc(&actual_res_exact);

    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_10() {
    let depth = 5;

    let (lon_min, lat_min, lon_max, lat_max) = (
        175_f64.to_radians(),
        -90_f64.to_radians(),
        185.0_f64.to_radians(),
        90_f64.to_radians(),
    );

    let expected_res_exact: [u64; 468] = [
        1360, 1361, 1363, 1364, 1365, 1366, 1367, 1372, 1373, 1374, 1375, 1396, 1397, 1398, 1399,
        1404, 1405, 1406, 1407, 1492, 1493, 1494, 1495, 1500, 1501, 1502, 1503, 1524, 1525, 1526,
        1527, 1532, 1533, 1535, 1877, 1879, 1885, 1887, 1909, 1911, 1917, 1919, 2005, 2007, 2013,
        2015, 2037, 2039, 2045, 2047, 2720, 2722, 2723, 2728, 2729, 2730, 2731, 2732, 2733, 2734,
        2735, 2744, 2745, 2746, 2747, 2748, 2749, 2750, 2751, 2792, 2793, 2794, 2795, 2796, 2797,
        2798, 2799, 2808, 2809, 2810, 2811, 2812, 2814, 2815, 2986, 2987, 2990, 2991, 3002, 3003,
        3006, 3007, 3050, 3051, 3054, 3055, 3066, 3067, 3070, 3071, 6144, 6145, 6146, 6147, 6148,
        6149, 6150, 6151, 6152, 6153, 6154, 6155, 6156, 6157, 6158, 6159, 6160, 6162, 6163, 6168,
        6169, 6170, 6171, 6172, 6174, 6175, 6176, 6177, 6179, 6180, 6181, 6182, 6183, 6188, 6189,
        6191, 6192, 6193, 6194, 6195, 6196, 6197, 6198, 6199, 6200, 6201, 6202, 6203, 6204, 6205,
        6206, 6207, 6240, 6242, 6243, 6248, 6249, 6250, 6251, 6252, 6254, 6255, 6288, 6289, 6291,
        6292, 6293, 6294, 6295, 6300, 6301, 6303, 6336, 6337, 6338, 6339, 6340, 6341, 6342, 6343,
        6344, 6345, 6346, 6347, 6348, 6349, 6350, 6351, 6352, 6354, 6355, 6360, 6361, 6362, 6363,
        6364, 6366, 6367, 6368, 6369, 6371, 6372, 6373, 6374, 6375, 6380, 6381, 6383, 6384, 6385,
        6386, 6387, 6388, 6389, 6390, 6391, 6392, 6393, 6394, 6395, 6396, 6397, 6398, 6399, 6560,
        6562, 6563, 6568, 6569, 6570, 6571, 6572, 6574, 6575, 6736, 6737, 6739, 6740, 6741, 6742,
        6743, 6748, 6749, 6751, 6912, 6913, 6914, 6915, 6916, 6917, 6918, 6919, 6920, 6921, 6922,
        6923, 6924, 6925, 6926, 6927, 6928, 6930, 6931, 6936, 6937, 6938, 6939, 6940, 6942, 6943,
        6944, 6945, 6947, 6948, 6949, 6950, 6951, 6956, 6957, 6959, 6960, 6961, 6962, 6963, 6964,
        6965, 6966, 6967, 6968, 6969, 6970, 6971, 6972, 6973, 6974, 6975, 7008, 7010, 7011, 7016,
        7017, 7018, 7019, 7020, 7022, 7023, 7056, 7057, 7059, 7060, 7061, 7062, 7063, 7068, 7069,
        7071, 7104, 7105, 7106, 7107, 7108, 7109, 7110, 7111, 7112, 7113, 7114, 7115, 7116, 7117,
        7118, 7119, 7120, 7122, 7123, 7128, 7129, 7130, 7131, 7132, 7134, 7135, 7136, 7137, 7139,
        7140, 7141, 7142, 7143, 7148, 7149, 7151, 7152, 7153, 7154, 7155, 7156, 7157, 7158, 7159,
        7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 9216, 9217, 9220, 9221, 9232, 9233, 9236,
        9237, 9280, 9281, 9284, 9285, 9296, 9297, 9300, 9301, 9472, 9473, 9475, 9476, 9477, 9478,
        9479, 9488, 9489, 9490, 9491, 9492, 9493, 9494, 9495, 9536, 9537, 9538, 9539, 9540, 9541,
        9542, 9543, 9552, 9553, 9554, 9555, 9556, 9557, 9558, 9559, 9564, 9565, 9567, 10240, 10242,
        10248, 10250, 10272, 10274, 10280, 10282, 10368, 10370, 10376, 10378, 10400, 10402, 10408,
        10410, 10752, 10754, 10755, 10760, 10761, 10762, 10763, 10784, 10785, 10786, 10787, 10792,
        10793, 10794, 10795, 10880, 10881, 10882, 10883, 10888, 10889, 10890, 10891, 10912, 10913,
        10914, 10915, 10920, 10921, 10922, 10923, 10924, 10926, 10927,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);

    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

#[test]
fn testok_zone_11() {
    let depth = 5;

    let (lon_min, lat_min, lon_max, lat_max) = (
        0_f64.to_radians(),
        -5_f64.to_radians(),
        360.0_f64.to_radians(),
        5_f64.to_radians(),
    );

    let expected_res_exact: [u64; 1408] = [
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 32, 1024, 1025, 1026, 1027, 1028, 1029, 1030,
        1031, 1032, 1033, 1034, 1035, 1036, 1040, 1056, 2048, 2049, 2050, 2051, 2052, 2053, 2054,
        2055, 2056, 2057, 2058, 2059, 2060, 2064, 2080, 3072, 3073, 3074, 3075, 3076, 3077, 3078,
        3079, 3080, 3081, 3082, 3083, 3084, 3088, 3104, 4319, 4335, 4339, 4340, 4341, 4342, 4343,
        4344, 4345, 4346, 4347, 4348, 4349, 4350, 4351, 4383, 4399, 4403, 4404, 4405, 4406, 4407,
        4408, 4409, 4410, 4411, 4412, 4413, 4414, 4415, 4419, 4420, 4421, 4422, 4423, 4424, 4425,
        4426, 4427, 4428, 4429, 4430, 4431, 4432, 4433, 4434, 4435, 4436, 4437, 4438, 4439, 4440,
        4441, 4442, 4443, 4444, 4445, 4446, 4447, 4448, 4449, 4450, 4451, 4452, 4453, 4454, 4455,
        4456, 4457, 4458, 4459, 4460, 4461, 4462, 4463, 4464, 4465, 4466, 4467, 4468, 4469, 4470,
        4471, 4472, 4473, 4474, 4475, 4476, 4483, 4484, 4485, 4486, 4487, 4488, 4489, 4490, 4491,
        4492, 4493, 4494, 4495, 4496, 4497, 4498, 4499, 4500, 4501, 4502, 4503, 4504, 4505, 4506,
        4507, 4508, 4509, 4510, 4511, 4512, 4513, 4514, 4515, 4516, 4517, 4518, 4519, 4520, 4521,
        4522, 4523, 4524, 4525, 4526, 4527, 4528, 4529, 4530, 4531, 4532, 4533, 4534, 4535, 4536,
        4537, 4538, 4539, 4540, 4544, 4545, 4546, 4547, 4548, 4549, 4550, 4551, 4552, 4553, 4554,
        4555, 4556, 4560, 4576, 4639, 4655, 4659, 4660, 4661, 4662, 4663, 4664, 4665, 4666, 4667,
        4668, 4669, 4670, 4671, 4675, 4676, 4677, 4678, 4679, 4680, 4681, 4682, 4683, 4684, 4685,
        4686, 4687, 4688, 4689, 4690, 4691, 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699, 4700,
        4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715,
        4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4729, 4730,
        4731, 4732, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751,
        4752, 4753, 4754, 4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4765, 4766,
        4767, 4768, 4769, 4770, 4771, 4772, 4773, 4774, 4775, 4776, 4777, 4778, 4779, 4780, 4781,
        4782, 4783, 4784, 4785, 4786, 4787, 4788, 4789, 4790, 4791, 4792, 4793, 4794, 4795, 4796,
        4800, 4801, 4802, 4803, 4804, 4805, 4806, 4807, 4808, 4809, 4810, 4811, 4812, 4816, 4832,
        4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4880, 4896,
        5343, 5359, 5363, 5364, 5365, 5366, 5367, 5368, 5369, 5370, 5371, 5372, 5373, 5374, 5375,
        5407, 5423, 5427, 5428, 5429, 5430, 5431, 5432, 5433, 5434, 5435, 5436, 5437, 5438, 5439,
        5443, 5444, 5445, 5446, 5447, 5448, 5449, 5450, 5451, 5452, 5453, 5454, 5455, 5456, 5457,
        5458, 5459, 5460, 5461, 5462, 5463, 5464, 5465, 5466, 5467, 5468, 5469, 5470, 5471, 5472,
        5473, 5474, 5475, 5476, 5477, 5478, 5479, 5480, 5481, 5482, 5483, 5484, 5485, 5486, 5487,
        5488, 5489, 5490, 5491, 5492, 5493, 5494, 5495, 5496, 5497, 5498, 5499, 5500, 5507, 5508,
        5509, 5510, 5511, 5512, 5513, 5514, 5515, 5516, 5517, 5518, 5519, 5520, 5521, 5522, 5523,
        5524, 5525, 5526, 5527, 5528, 5529, 5530, 5531, 5532, 5533, 5534, 5535, 5536, 5537, 5538,
        5539, 5540, 5541, 5542, 5543, 5544, 5545, 5546, 5547, 5548, 5549, 5550, 5551, 5552, 5553,
        5554, 5555, 5556, 5557, 5558, 5559, 5560, 5561, 5562, 5563, 5564, 5568, 5569, 5570, 5571,
        5572, 5573, 5574, 5575, 5576, 5577, 5578, 5579, 5580, 5584, 5600, 5663, 5679, 5683, 5684,
        5685, 5686, 5687, 5688, 5689, 5690, 5691, 5692, 5693, 5694, 5695, 5699, 5700, 5701, 5702,
        5703, 5704, 5705, 5706, 5707, 5708, 5709, 5710, 5711, 5712, 5713, 5714, 5715, 5716, 5717,
        5718, 5719, 5720, 5721, 5722, 5723, 5724, 5725, 5726, 5727, 5728, 5729, 5730, 5731, 5732,
        5733, 5734, 5735, 5736, 5737, 5738, 5739, 5740, 5741, 5742, 5743, 5744, 5745, 5746, 5747,
        5748, 5749, 5750, 5751, 5752, 5753, 5754, 5755, 5756, 5763, 5764, 5765, 5766, 5767, 5768,
        5769, 5770, 5771, 5772, 5773, 5774, 5775, 5776, 5777, 5778, 5779, 5780, 5781, 5782, 5783,
        5784, 5785, 5786, 5787, 5788, 5789, 5790, 5791, 5792, 5793, 5794, 5795, 5796, 5797, 5798,
        5799, 5800, 5801, 5802, 5803, 5804, 5805, 5806, 5807, 5808, 5809, 5810, 5811, 5812, 5813,
        5814, 5815, 5816, 5817, 5818, 5819, 5820, 5824, 5825, 5826, 5827, 5828, 5829, 5830, 5831,
        5832, 5833, 5834, 5835, 5836, 5840, 5856, 5888, 5889, 5890, 5891, 5892, 5893, 5894, 5895,
        5896, 5897, 5898, 5899, 5900, 5904, 5920, 6367, 6383, 6387, 6388, 6389, 6390, 6391, 6392,
        6393, 6394, 6395, 6396, 6397, 6398, 6399, 6431, 6447, 6451, 6452, 6453, 6454, 6455, 6456,
        6457, 6458, 6459, 6460, 6461, 6462, 6463, 6467, 6468, 6469, 6470, 6471, 6472, 6473, 6474,
        6475, 6476, 6477, 6478, 6479, 6480, 6481, 6482, 6483, 6484, 6485, 6486, 6487, 6488, 6489,
        6490, 6491, 6492, 6493, 6494, 6495, 6496, 6497, 6498, 6499, 6500, 6501, 6502, 6503, 6504,
        6505, 6506, 6507, 6508, 6509, 6510, 6511, 6512, 6513, 6514, 6515, 6516, 6517, 6518, 6519,
        6520, 6521, 6522, 6523, 6524, 6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6539, 6540,
        6541, 6542, 6543, 6544, 6545, 6546, 6547, 6548, 6549, 6550, 6551, 6552, 6553, 6554, 6555,
        6556, 6557, 6558, 6559, 6560, 6561, 6562, 6563, 6564, 6565, 6566, 6567, 6568, 6569, 6570,
        6571, 6572, 6573, 6574, 6575, 6576, 6577, 6578, 6579, 6580, 6581, 6582, 6583, 6584, 6585,
        6586, 6587, 6588, 6592, 6593, 6594, 6595, 6596, 6597, 6598, 6599, 6600, 6601, 6602, 6603,
        6604, 6608, 6624, 6687, 6703, 6707, 6708, 6709, 6710, 6711, 6712, 6713, 6714, 6715, 6716,
        6717, 6718, 6719, 6723, 6724, 6725, 6726, 6727, 6728, 6729, 6730, 6731, 6732, 6733, 6734,
        6735, 6736, 6737, 6738, 6739, 6740, 6741, 6742, 6743, 6744, 6745, 6746, 6747, 6748, 6749,
        6750, 6751, 6752, 6753, 6754, 6755, 6756, 6757, 6758, 6759, 6760, 6761, 6762, 6763, 6764,
        6765, 6766, 6767, 6768, 6769, 6770, 6771, 6772, 6773, 6774, 6775, 6776, 6777, 6778, 6779,
        6780, 6787, 6788, 6789, 6790, 6791, 6792, 6793, 6794, 6795, 6796, 6797, 6798, 6799, 6800,
        6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815,
        6816, 6817, 6818, 6819, 6820, 6821, 6822, 6823, 6824, 6825, 6826, 6827, 6828, 6829, 6830,
        6831, 6832, 6833, 6834, 6835, 6836, 6837, 6838, 6839, 6840, 6841, 6842, 6843, 6844, 6848,
        6849, 6850, 6851, 6852, 6853, 6854, 6855, 6856, 6857, 6858, 6859, 6860, 6864, 6880, 6912,
        6913, 6914, 6915, 6916, 6917, 6918, 6919, 6920, 6921, 6922, 6923, 6924, 6928, 6944, 7391,
        7407, 7411, 7412, 7413, 7414, 7415, 7416, 7417, 7418, 7419, 7420, 7421, 7422, 7423, 7455,
        7471, 7475, 7476, 7477, 7478, 7479, 7480, 7481, 7482, 7483, 7484, 7485, 7486, 7487, 7491,
        7492, 7493, 7494, 7495, 7496, 7497, 7498, 7499, 7500, 7501, 7502, 7503, 7504, 7505, 7506,
        7507, 7508, 7509, 7510, 7511, 7512, 7513, 7514, 7515, 7516, 7517, 7518, 7519, 7520, 7521,
        7522, 7523, 7524, 7525, 7526, 7527, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7535, 7536,
        7537, 7538, 7539, 7540, 7541, 7542, 7543, 7544, 7545, 7546, 7547, 7548, 7555, 7556, 7557,
        7558, 7559, 7560, 7561, 7562, 7563, 7564, 7565, 7566, 7567, 7568, 7569, 7570, 7571, 7572,
        7573, 7574, 7575, 7576, 7577, 7578, 7579, 7580, 7581, 7582, 7583, 7584, 7585, 7586, 7587,
        7588, 7589, 7590, 7591, 7592, 7593, 7594, 7595, 7596, 7597, 7598, 7599, 7600, 7601, 7602,
        7603, 7604, 7605, 7606, 7607, 7608, 7609, 7610, 7611, 7612, 7616, 7617, 7618, 7619, 7620,
        7621, 7622, 7623, 7624, 7625, 7626, 7627, 7628, 7632, 7648, 7711, 7727, 7731, 7732, 7733,
        7734, 7735, 7736, 7737, 7738, 7739, 7740, 7741, 7742, 7743, 7747, 7748, 7749, 7750, 7751,
        7752, 7753, 7754, 7755, 7756, 7757, 7758, 7759, 7760, 7761, 7762, 7763, 7764, 7765, 7766,
        7767, 7768, 7769, 7770, 7771, 7772, 7773, 7774, 7775, 7776, 7777, 7778, 7779, 7780, 7781,
        7782, 7783, 7784, 7785, 7786, 7787, 7788, 7789, 7790, 7791, 7792, 7793, 7794, 7795, 7796,
        7797, 7798, 7799, 7800, 7801, 7802, 7803, 7804, 7811, 7812, 7813, 7814, 7815, 7816, 7817,
        7818, 7819, 7820, 7821, 7822, 7823, 7824, 7825, 7826, 7827, 7828, 7829, 7830, 7831, 7832,
        7833, 7834, 7835, 7836, 7837, 7838, 7839, 7840, 7841, 7842, 7843, 7844, 7845, 7846, 7847,
        7848, 7849, 7850, 7851, 7852, 7853, 7854, 7855, 7856, 7857, 7858, 7859, 7860, 7861, 7862,
        7863, 7864, 7865, 7866, 7867, 7868, 7872, 7873, 7874, 7875, 7876, 7877, 7878, 7879, 7880,
        7881, 7882, 7883, 7884, 7888, 7904, 7936, 7937, 7938, 7939, 7940, 7941, 7942, 7943, 7944,
        7945, 7946, 7947, 7948, 7952, 7968, 9183, 9199, 9203, 9204, 9205, 9206, 9207, 9208, 9209,
        9210, 9211, 9212, 9213, 9214, 9215, 10207, 10223, 10227, 10228, 10229, 10230, 10231, 10232,
        10233, 10234, 10235, 10236, 10237, 10238, 10239, 11231, 11247, 11251, 11252, 11253, 11254,
        11255, 11256, 11257, 11258, 11259, 11260, 11261, 11262, 11263, 12255, 12271, 12275, 12276,
        12277, 12278, 12279, 12280, 12281, 12282, 12283, 12284, 12285, 12286, 12287,
    ];
    let actual_res_exact = get(depth).zone_coverage([lon_min, lat_min], [lon_max, lat_max]);
    // to_aladin_moc(&actual_res_exact);

    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    assert_eq!(
        actual_res_exact.flat_iter().collect::<Vec<_>>(),
        expected_res_exact
    );
    // for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
    //     assert_eq!(h1, *h2);
    // }
}

fn to_radians(lonlats: &mut [(f64, f64)]) {
    for (lon, lat) in lonlats.iter_mut() {
        *lon = lon.to_radians();
        *lat = lat.to_radians();
    }
}

#[allow(dead_code)]
fn to_degrees(lonlats: &mut [(f64, f64)]) {
    for (lon, lat) in lonlats.iter_mut() {
        *lon = lon.to_degrees();
        *lat = lat.to_degrees();
    }
}

#[test]
fn testok_bmoc_not() {
    let actual_res = get(3).cone_coverage_approx_custom(
        4,
        Degrees(36.80105218, 56.78028536),
        14.93_f64.to_radians(),
    );
    /*println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    println!("@@@@@ HIERARCH VIEW, NOT");*/
    let complement = actual_res.not();
    /*for cell in complement.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    let org = complement.not();
    assert_eq!(actual_res, org);
    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
}

#[test]
fn test_prec_1() {
    let lon_deg = 179.99999999999997_f64;
    let lat_deg = 41.813964843754924_f64;
    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get(0);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    let layer_0 = get(0);
    assert_eq!(1, layer_0.hash(Degrees(lon_deg, lat_deg)));
}

#[test]
fn test_prec_2() {
    let lon_deg = 359.99999999999994_f64;
    let lat_deg = 41.81031502783791_f64;
    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_1 = get(1);
    layer_1.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_1.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);*/
    let layer_1 = get(1);
    assert_eq!(13, layer_1.hash(Degrees(lon_deg, lat_deg)));
}

#[test]
fn test_prec_3() {
    let lon_deg = 359.99999999999994_f64; //359.99999999999994_f64;
    let lat_deg = 41.81031489577861_f64; //41.81031489577857_f64;

    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get(6);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);/*
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    println!("hash: {}", layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));*/

    let layer_6 = get(6);
    assert_eq!(13653, layer_6.hash(Degrees(lon_deg, lat_deg)));
}

/*
#[test]
fn test_prec_4() {
  let lon_deg = 292.49999999999994_f64; //359.99999999999994_f64;
  let lat_deg = 41.810314895778546_f64; //41.81031489577857_f64;

  let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
  // println!("proj_x: {:.17}, proj_y: {:.17}", xy.0, xy.1);
  xy.0 = ensures_x_is_positive(xy.0);
  let layer_0 = get(0);
  layer_0.shift_rotate_scale(&mut xy);
  // println!("x: {}, y: {}", xy.0, xy.1);
  let ij = discretize(xy);
  // println!("i: {}, j: {}", ij.0, ij.1);
  let ij_d0c = layer_0.base_cell_coos(&ij);
  // println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);/*
  let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1, ij, xy/*, lon, lat*/);
  // println!("d0h_bits: {}", d0h_bits);*/
  // println!("hash: {}", layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));

  let layer_2 = get(2);
  assert_eq!(56, layer_2.hash(lon_deg.to_radians(), lat_deg.to_radians()));
}
 */

#[test]
fn test_bilinear_interpolation() {
    let lon_deg = 89.18473162_f64; // 322.99297784_f64;// 324.8778822_f64
    let lat_deg = -28.04159707_f64; // 39.9302924_f64;// -41.08635508_f64
    let res = get(1).bilinear_interpolation(Degrees(lon_deg, lat_deg));
    // println!("{:?}", res);
    // // Result with previous version of hash_dxdy_v1 !
    // assert_eq!(
    //     res,
    //     [
    //         (20, 0.0),
    //         (38, 0.1661686383097217),
    //         (33, 0.2024027885319438),
    //         (20, 0.6314285731583344)
    //     ]
    // );
    // Result with previous version of hash_dxdy_v2 !
    // assert_eq!(
    //     res,
    //     [
    //         (20, 0.0),
    //         (38, 0.1661686383097213),
    //         (33, 0.20240278853194385),
    //         (20, 0.6314285731583348)
    //     ]
    // );
    assert!(
        res.iter()
            .zip(&[
                (20, 0.0),
                (38, 0.1661686383097213),
                (33, 0.20240278853194385),
                (20, 0.6314285731583348)
            ])
            .all(|((ai, af), (bi, bf))| ai == bi && (af - bf).abs() < 0.0000001),
        "left: {res:?}, right: [(20, 0.0), (38, 0.1661686383097213), (33, 0.20240278853194385), (20, 0.6314285731583348)]"
    );
}

#[test]
fn test_bilinear_interpolation_2() {
    let lon_deg = 83.633478_f64;
    let lat_deg = 22.015110_f64;
    let res = get(18).bilinear_interpolation(Degrees(lon_deg, lat_deg));
    // println!("{:?}", res);
    // // Result with previous version of hash_dxdy_v1 !
    // assert_eq!(
    //     res,
    //     [
    //         (405766747916, 0.5757471135241182),
    //         (405766747917, 0.3604806280107034),
    //         (405766747918, 0.039217694696856834),
    //         (405766747919, 0.024554563768321474)
    //     ]
    // );
    assert!(
        res.iter()
            .zip(&[
                (405766747916, 0.5757471135599139),
                (405766747917, 0.3604806280331154),
                (405766747918, 0.039217694661061175),
                (405766747919, 0.024554563745909478)
            ])
            .all(|((ai, af), (bi, bf))| ai == bi && (af - bf).abs() < 0.0000001),
        "left: {res:?}, right: [(405766747916, 0.5757471135599139), (405766747917, 0.3604806280331154), (405766747918, 0.039217694661061175), (405766747919, 0.024554563745909478)]"
    );
}

#[test]
#[allow(clippy::approx_constant)]
fn test_bilinear_interpolation_3() {
    let lon_rad = [0.17453293_f64, 0.43633231_f64, 0.0_f64];
    let lat_rad = [0.08726646_f64, 0.17453293_f64, 0.78539816_f64];
    let depth = 5;
    for (lon, lat) in lon_rad.into_iter().zip(lat_rad) {
        let res = get(depth).bilinear_interpolation(Degrees(lon, lat));
        println!("{res:?}");
    }
}

#[test]
fn test_gen_file() -> std::io::Result<()> {
    use super::get;
    use std::io::prelude::*;

    let depth = 8;
    let layer = get(depth);
    let n_cells = layer.n_hash();

    let mut s = String::with_capacity(8 * 1024);
    s.push_str("i,ra,dec\n");
    for h in 0..n_cells {
        let LonLat { lon, lat } = layer.center(h);
        s.push_str(&format!(
            "{},{},{}\n",
            &h,
            &lon.to_degrees(),
            &lat.to_degrees()
        ));
    }

    let mut file = std::fs::File::create(format!("hpx.{depth}.csv"))?;
    file.write_all(s.as_bytes())?;
    Ok(())
}

#[test]
fn test_to_uniq() {
    for depth in 0..8 {
        for idx in 0..checked::n_hash(depth) {
            assert_eq!((depth, idx), from_uniq(checked::to_uniq(depth, idx)));
        }
    }
}

#[test]
fn test_to_uniq_ivoa() {
    for depth in 0..8 {
        for idx in 0..checked::n_hash(depth) {
            assert_eq!(
                (depth, idx),
                from_uniq_ivoa(checked::to_uniq_ivoa(depth, idx))
            );
        }
    }
}

#[test]
fn test_kth_neighbors() {
    let depth = 8;
    let nside = checked::nside(depth);
    let k = 4;
    let layer = get(depth);

    let expected = [
        // NPC
        [
            vec![
                589817, 589819, 589821, 589820, 283985, 283987, 283993, 283995, 284017, 284019,
                284025, 284028, 284029, 371387, 371385, 371384, 371373, 371372, 371369, 371368, 40,
                41, 44, 45, 56, 57, 51, 49, 27, 25, 19, 17,
            ],
            vec![
                109246, 109244, 109238, 109235, 109234, 109223, 109222, 109219, 109218, 393196,
                393198, 393213, 393212, 393209, 393208, 393197, 21828, 21830, 21836, 21838, 21860,
                21862, 21868, 21869, 21880, 21881, 21884, 21885,
            ],
            vec![
                262108, 262110, 262132, 262134, 262140, 262142, 262109, 196601, 196603, 196605,
                196604, 131043, 131049, 131051, 131063, 131062, 131059, 131058, 131047, 131046,
                65478, 65484, 65486, 65508, 65510, 65516, 65518, 65495, 65494, 65491, 65490, 65479,
            ],
            vec![
                327635, 327641, 327643, 327665, 327667, 327673, 327675, 327639, 327638, 218452,
                218454, 218460, 218462, 218484, 218486, 218487, 43707, 43705, 43699, 43697, 43675,
                43673, 43667, 43666, 43655, 43654, 43651, 43650,
            ],
        ],
        [
            vec![
                655353, 655355, 655357, 655356, 349521, 349523, 349529, 349531, 349553, 349555,
                349561, 349564, 349565, 436923, 436921, 436920, 436909, 436908, 436905, 436904,
                65576, 65577, 65580, 65581, 65592, 65593, 65587, 65585, 65563, 65561, 65555, 65553,
            ],
            vec![
                174782, 174780, 174774, 174771, 174770, 174759, 174758, 174755, 174754, 458732,
                458734, 458749, 458748, 458745, 458744, 458733, 87364, 87366, 87372, 87374, 87396,
                87398, 87404, 87405, 87416, 87417, 87420, 87421,
            ],
            vec![
                65500, 65502, 65524, 65526, 65532, 65534, 65501, 262137, 262139, 262141, 262140,
                196579, 196585, 196587, 196599, 196598, 196595, 196594, 196583, 196582, 131014,
                131020, 131022, 131044, 131046, 131052, 131054, 131031, 131030, 131027, 131026,
                131015,
            ],
            vec![
                393171, 393177, 393179, 393201, 393203, 393209, 393211, 393175, 393174, 21844,
                21846, 21852, 21854, 21876, 21878, 21879, 109243, 109241, 109235, 109233, 109211,
                109209, 109203, 109202, 109191, 109190, 109187, 109186,
            ],
        ],
        [
            vec![
                720889, 720891, 720893, 720892, 415057, 415059, 415065, 415067, 415089, 415091,
                415097, 415100, 415101, 502459, 502457, 502456, 502445, 502444, 502441, 502440,
                131112, 131113, 131116, 131117, 131128, 131129, 131123, 131121, 131099, 131097,
                131091, 131089,
            ],
            vec![
                240318, 240316, 240310, 240307, 240306, 240295, 240294, 240291, 240290, 524268,
                524270, 524285, 524284, 524281, 524280, 524269, 152900, 152902, 152908, 152910,
                152932, 152934, 152940, 152941, 152952, 152953, 152956, 152957,
            ],
            vec![
                131036, 131038, 131060, 131062, 131068, 131070, 131037, 65529, 65531, 65533, 65532,
                262115, 262121, 262123, 262135, 262134, 262131, 262130, 262119, 262118, 196550,
                196556, 196558, 196580, 196582, 196588, 196590, 196567, 196566, 196563, 196562,
                196551,
            ],
            vec![
                458707, 458713, 458715, 458737, 458739, 458745, 458747, 458711, 458710, 87380,
                87382, 87388, 87390, 87412, 87414, 87415, 174779, 174777, 174771, 174769, 174747,
                174745, 174739, 174738, 174727, 174726, 174723, 174722,
            ],
        ],
        [
            vec![
                786425, 786427, 786429, 786428, 480593, 480595, 480601, 480603, 480625, 480627,
                480633, 480636, 480637, 305851, 305849, 305848, 305837, 305836, 305833, 305832,
                196648, 196649, 196652, 196653, 196664, 196665, 196659, 196657, 196635, 196633,
                196627, 196625,
            ],
            vec![
                43710, 43708, 43702, 43699, 43698, 43687, 43686, 43683, 43682, 327660, 327662,
                327677, 327676, 327673, 327672, 327661, 218436, 218438, 218444, 218446, 218468,
                218470, 218476, 218477, 218488, 218489, 218492, 218493,
            ],
            vec![
                196572, 196574, 196596, 196598, 196604, 196606, 196573, 131065, 131067, 131069,
                131068, 65507, 65513, 65515, 65527, 65526, 65523, 65522, 65511, 65510, 262086,
                262092, 262094, 262116, 262118, 262124, 262126, 262103, 262102, 262099, 262098,
                262087,
            ],
            vec![
                524243, 524249, 524251, 524273, 524275, 524281, 524283, 524247, 524246, 152916,
                152918, 152924, 152926, 152948, 152950, 152951, 240315, 240313, 240307, 240305,
                240283, 240281, 240275, 240274, 240263, 240262, 240259, 240258,
            ],
        ],
        // EQR
        [
            vec![
                742737, 742739, 742745, 742747, 742769, 742771, 742777, 742780, 742781, 567995,
                567993, 567992, 567981, 567980, 567977, 567976, 262184, 262185, 262188, 262189,
                262200, 262201, 262195, 262193, 262171, 262169, 262163, 262161,
            ],
            vec![
                40, 41, 44, 38, 36, 14, 12, 6, 4, 371374, 371372, 371369, 371368, 589804, 589806,
                589821, 589820, 589817, 589816, 589805, 283972, 283974, 283980, 283982, 284004,
                284006, 284012, 284013, 284024, 284025, 284028, 284029,
            ],
            vec![
                218436, 218438, 218439, 218450, 218451, 218454, 218455, 43694, 43692, 43686, 43684,
                43662, 43660, 43654, 43651, 43650, 327622, 327628, 327630, 327652, 327654, 327660,
                327662, 327639, 327638, 327635, 327634, 327623,
            ],
            vec![
                786387, 786393, 786395, 786417, 786419, 786425, 786427, 786391, 786390, 480593,
                480595, 480598, 480599, 196610, 196611, 196614, 196615, 196626, 196627, 196625,
                305851, 305849, 305843, 305841, 305819, 305817, 305811, 305810, 305799, 305798,
                305795, 305794,
            ],
        ],
        [
            vec![
                546129, 546131, 546137, 546139, 546161, 546163, 546169, 546172, 546173, 633531,
                633529, 633528, 633517, 633516, 633513, 633512, 327720, 327721, 327724, 327725,
                327736, 327737, 327731, 327729, 327707, 327705, 327699, 327697,
            ],
            vec![
                65576, 65577, 65580, 65574, 65572, 65550, 65548, 65542, 65540, 436910, 436908,
                436905, 436904, 655340, 655342, 655357, 655356, 655353, 655352, 655341, 349508,
                349510, 349516, 349518, 349540, 349542, 349548, 349549, 349560, 349561, 349564,
                349565,
            ],
            vec![
                21828, 21830, 21831, 21842, 21843, 21846, 21847, 109230, 109228, 109222, 109220,
                109198, 109196, 109190, 109187, 109186, 393158, 393164, 393166, 393188, 393190,
                393196, 393198, 393175, 393174, 393171, 393170, 393159,
            ],
            vec![
                589779, 589785, 589787, 589809, 589811, 589817, 589819, 589783, 589782, 283985,
                283987, 283990, 283991, 2, 3, 6, 7, 18, 19, 17, 371387, 371385, 371379, 371377,
                371355, 371353, 371347, 371346, 371335, 371334, 371331, 371330,
            ],
        ],
        [
            vec![
                611665, 611667, 611673, 611675, 611697, 611699, 611705, 611708, 611709, 699067,
                699065, 699064, 699053, 699052, 699049, 699048, 393256, 393257, 393260, 393261,
                393272, 393273, 393267, 393265, 393243, 393241, 393235, 393233,
            ],
            vec![
                131112, 131113, 131116, 131110, 131108, 131086, 131084, 131078, 131076, 502446,
                502444, 502441, 502440, 720876, 720878, 720893, 720892, 720889, 720888, 720877,
                415044, 415046, 415052, 415054, 415076, 415078, 415084, 415085, 415096, 415097,
                415100, 415101,
            ],
            vec![
                87364, 87366, 87367, 87378, 87379, 87382, 87383, 174766, 174764, 174758, 174756,
                174734, 174732, 174726, 174723, 174722, 458694, 458700, 458702, 458724, 458726,
                458732, 458734, 458711, 458710, 458707, 458706, 458695,
            ],
            vec![
                655315, 655321, 655323, 655345, 655347, 655353, 655355, 655319, 655318, 349521,
                349523, 349526, 349527, 65538, 65539, 65542, 65543, 65554, 65555, 65553, 436923,
                436921, 436915, 436913, 436891, 436889, 436883, 436882, 436871, 436870, 436867,
                436866,
            ],
        ],
        [
            vec![
                677201, 677203, 677209, 677211, 677233, 677235, 677241, 677244, 677245, 764603,
                764601, 764600, 764589, 764588, 764585, 764584, 458792, 458793, 458796, 458797,
                458808, 458809, 458803, 458801, 458779, 458777, 458771, 458769,
            ],
            vec![
                196648, 196649, 196652, 196646, 196644, 196622, 196620, 196614, 196612, 305838,
                305836, 305833, 305832, 786412, 786414, 786429, 786428, 786425, 786424, 786413,
                480580, 480582, 480588, 480590, 480612, 480614, 480620, 480621, 480632, 480633,
                480636, 480637,
            ],
            vec![
                152900, 152902, 152903, 152914, 152915, 152918, 152919, 240302, 240300, 240294,
                240292, 240270, 240268, 240262, 240259, 240258, 524230, 524236, 524238, 524260,
                524262, 524268, 524270, 524247, 524246, 524243, 524242, 524231,
            ],
            vec![
                720851, 720857, 720859, 720881, 720883, 720889, 720891, 720855, 720854, 415057,
                415059, 415062, 415063, 131074, 131075, 131078, 131079, 131090, 131091, 131089,
                502459, 502457, 502451, 502449, 502427, 502425, 502419, 502418, 502407, 502406,
                502403, 502402,
            ],
        ],
        // SPC
        [
            vec![
                655362, 655363, 655366, 655364, 720904, 720905, 720908, 720909, 720920, 720921,
                720924, 720918, 720916, 589858, 589859, 589857, 589835, 589833, 589827, 589825,
                524328, 524329, 524332, 524333, 524344, 524345, 524339, 524337, 524315, 524313,
                524307, 524305,
            ],
            vec![
                327720, 327721, 327724, 327718, 327716, 327694, 327692, 327686, 327684, 633515,
                633513, 633507, 633505, 633483, 633481, 633480, 546116, 546118, 546124, 546126,
                546148, 546150, 546156, 546157, 546168, 546169, 546172, 546173,
            ],
            vec![
                283972, 283974, 283975, 283986, 283987, 283990, 283991, 2, 3, 6, 4, 371374, 371372,
                371366, 371364, 371342, 371340, 371334, 371331, 371330, 589766, 589772, 589774,
                589796, 589798, 589804, 589806, 589783, 589782, 589779, 589778, 589767,
            ],
            vec![
                742721, 742723, 742729, 742732, 742733, 742744, 742745, 742748, 742749, 262146,
                262147, 262150, 262151, 262162, 262163, 262161, 567995, 567993, 567987, 567985,
                567963, 567961, 567955, 567954, 567943, 567942, 567939, 567938,
            ],
        ],
        [
            vec![
                720898, 720899, 720902, 720900, 524296, 524297, 524300, 524301, 524312, 524313,
                524316, 524310, 524308, 655394, 655395, 655393, 655371, 655369, 655363, 655361,
                589864, 589865, 589868, 589869, 589880, 589881, 589875, 589873, 589851, 589849,
                589843, 589841,
            ],
            vec![
                393256, 393257, 393260, 393254, 393252, 393230, 393228, 393222, 393220, 699051,
                699049, 699043, 699041, 699019, 699017, 699016, 611652, 611654, 611660, 611662,
                611684, 611686, 611692, 611693, 611704, 611705, 611708, 611709,
            ],
            vec![
                349508, 349510, 349511, 349522, 349523, 349526, 349527, 65538, 65539, 65542, 65540,
                436910, 436908, 436902, 436900, 436878, 436876, 436870, 436867, 436866, 655302,
                655308, 655310, 655332, 655334, 655340, 655342, 655319, 655318, 655315, 655314,
                655303,
            ],
            vec![
                546113, 546115, 546121, 546124, 546125, 546136, 546137, 546140, 546141, 327682,
                327683, 327686, 327687, 327698, 327699, 327697, 633531, 633529, 633523, 633521,
                633499, 633497, 633491, 633490, 633479, 633478, 633475, 633474,
            ],
        ],
        [
            vec![
                524290, 524291, 524294, 524292, 589832, 589833, 589836, 589837, 589848, 589849,
                589852, 589846, 589844, 720930, 720931, 720929, 720907, 720905, 720899, 720897,
                655400, 655401, 655404, 655405, 655416, 655417, 655411, 655409, 655387, 655385,
                655379, 655377,
            ],
            vec![
                458792, 458793, 458796, 458790, 458788, 458766, 458764, 458758, 458756, 764587,
                764585, 764579, 764577, 764555, 764553, 764552, 677188, 677190, 677196, 677198,
                677220, 677222, 677228, 677229, 677240, 677241, 677244, 677245,
            ],
            vec![
                415044, 415046, 415047, 415058, 415059, 415062, 415063, 131074, 131075, 131078,
                131076, 502446, 502444, 502438, 502436, 502414, 502412, 502406, 502403, 502402,
                720838, 720844, 720846, 720868, 720870, 720876, 720878, 720855, 720854, 720851,
                720850, 720839,
            ],
            vec![
                611649, 611651, 611657, 611660, 611661, 611672, 611673, 611676, 611677, 393218,
                393219, 393222, 393223, 393234, 393235, 393233, 699067, 699065, 699059, 699057,
                699035, 699033, 699027, 699026, 699015, 699014, 699011, 699010,
            ],
        ],
        [
            vec![
                589826, 589827, 589830, 589828, 655368, 655369, 655372, 655373, 655384, 655385,
                655388, 655382, 655380, 524322, 524323, 524321, 524299, 524297, 524291, 524289,
                720936, 720937, 720940, 720941, 720952, 720953, 720947, 720945, 720923, 720921,
                720915, 720913,
            ],
            vec![
                262184, 262185, 262188, 262182, 262180, 262158, 262156, 262150, 262148, 567979,
                567977, 567971, 567969, 567947, 567945, 567944, 742724, 742726, 742732, 742734,
                742756, 742758, 742764, 742765, 742776, 742777, 742780, 742781,
            ],
            vec![
                480580, 480582, 480583, 480594, 480595, 480598, 480599, 196610, 196611, 196614,
                196612, 305838, 305836, 305830, 305828, 305806, 305804, 305798, 305795, 305794,
                786374, 786380, 786382, 786404, 786406, 786412, 786414, 786391, 786390, 786387,
                786386, 786375,
            ],
            vec![
                677185, 677187, 677193, 677196, 677197, 677208, 677209, 677212, 677213, 458754,
                458755, 458758, 458759, 458770, 458771, 458769, 764603, 764601, 764595, 764593,
                764571, 764569, 764563, 764562, 764551, 764550, 764547, 764546,
            ],
        ],
    ];

    for d0h in 0..12 {
        let expected_array = &expected[d0h as usize];

        let hs = layer.build_hash_from_parts(d0h, 1, 2);
        let neig = layer.kth_neighbors(hs, k);
        //to_aladin(depth, hs,neig.as_slice());
        let expected = &expected_array[0_usize];
        //if !expected.is_empty() {
        assert_eq!(&neig, expected);
        //}

        let he = layer.build_hash_from_parts(d0h, nside - 2, 2);
        let neig = layer.kth_neighbors(he, k);
        //to_aladin(depth, he,neig.as_slice());
        let expected = &expected_array[1_usize];
        //if !expected.is_empty() {
        assert_eq!(&neig, expected);
        //}

        let hn = layer.build_hash_from_parts(d0h, nside - 2, nside - 3);
        let neig = layer.kth_neighbors(hn, k);
        //to_aladin(depth, hn, neig.as_slice());
        let expected = &expected_array[2_usize];
        //if !expected.is_empty() {
        assert_eq!(&neig, expected);
        //}

        let hw = layer.build_hash_from_parts(d0h, 1, nside - 3);
        let neig = layer.kth_neighbors(hw, k);
        //to_aladin(depth, hw,neig.as_slice());
        let expected = &expected_array[3_usize];
        //if !expected.is_empty() {
        assert_eq!(&neig, expected);
        //}
    }
}

/*#[test]
fn test_xpm1_and_q () {
  let a = Layer::d0h_lh_in_d0c(std::f64::NAN, 0.0);
  eprintln!("{:?}", &a);
  let a = Layer::d0h_lh_in_d0c(std::f64::INFINITY, 0.0);
  eprintln!("{:?}", &a);
}*/
