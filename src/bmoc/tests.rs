use super::*;
use crate::coords::Degrees;
use crate::get;

#[test]
fn test_ok_allsky_and_empty_bmoc() {
    let bmoc_allsky = MutableBmoc::new_allsky(18);
    assert_eq!(bmoc_allsky.entries.len(), 12);
    let bmoc_empty = MutableBmoc::new_empty(18);
    assert_eq!(bmoc_empty.entries.len(), 0);
    assert_eq!(bmoc_allsky.not(), bmoc_empty);
    assert_eq!(bmoc_allsky, bmoc_empty.not());
}

#[test]
fn test_ok_bmoc_not_round_trip() {
    let mut poly_vertices_deg: [f64; 14] = [
        272.511081, -19.487278, 272.515300, -19.486595, 272.517029, -19.471442, 272.511714,
        -19.458837, 272.506430, -19.459001, 272.496401, -19.474322, 272.504821, -19.484924,
    ];
    for v in poly_vertices_deg.iter_mut() {
        *v = (*v).to_radians();
    }
    let vertices: Vec<(f64, f64)> = poly_vertices_deg
        .iter()
        .copied()
        .step_by(2)
        .zip(poly_vertices_deg.iter().copied().skip(1).step_by(2))
        .collect();
    let moc_org = get(14).polygon_coverage(vertices.as_slice(), true);
    let moc_not = moc_org.not();
    let moc_out = moc_not.not();
    println!("len: {}", moc_org.size());
    println!("len: {}", moc_not.size());
    println!("len: {}", moc_out.size());
    assert_eq!(moc_org, moc_out);
}

#[test]
fn test_ok_bmoc_union_and_not() {
    let mut poly1_vertices_deg: [f64; 14] = [
        272.511081, -19.487278, 272.515300, -19.486595, 272.517029, -19.471442, 272.511714,
        -19.458837, 272.506430, -19.459001, 272.496401, -19.474322, 272.504821, -19.484924,
    ];
    for v in poly1_vertices_deg.iter_mut() {
        *v = (*v).to_radians();
    }
    let mut poly2_vertices_deg: [f64; 8] = [
        272.630446, -19.234210, 272.637274, -19.248542, 272.638942, -19.231476, 272.630868,
        -19.226364,
    ];
    for v in poly2_vertices_deg.iter_mut() {
        *v = (*v).to_radians();
    }
    let v1: Vec<(f64, f64)> = poly1_vertices_deg
        .iter()
        .copied()
        .step_by(2)
        .zip(poly1_vertices_deg.iter().copied().skip(1).step_by(2))
        .collect();
    let v2: Vec<(f64, f64)> = poly2_vertices_deg
        .iter()
        .copied()
        .step_by(2)
        .zip(poly2_vertices_deg.iter().copied().skip(1).step_by(2))
        .collect();
    let moc1 = get(14).polygon_coverage(v1.as_slice(), true);
    let moc2 = get(14).polygon_coverage(v2.as_slice(), true);
    let not_moc1 = moc1.not();
    let not_moc2 = moc2.not();
    let union = moc1.or(&moc2);
    let not_inter = not_moc1.and(&not_moc2);
    let union_bis = not_inter.not();
    //to_aladin_moc(&moc1);
    //println!("\n");
    //to_aladin_moc(&moc2);
    //println!("\n");
    // to_aladin_moc(&union);
    //println!("\n");
    // to_aladin_moc(&union_bis);
    //println!("\n");
    assert_eq!(union, union_bis);
}

#[test]
fn test_ok_bmoc_xor_minus_coherency() {
    // No overlapping parts, so we do no test thoroughly XOR and MINUS!!
    let mut poly1_vertices_deg: [f64; 14] = [
        272.511081, -19.487278, 272.515300, -19.486595, 272.517029, -19.471442, 272.511714,
        -19.458837, 272.506430, -19.459001, 272.496401, -19.474322, 272.504821, -19.484924,
    ];
    for v in poly1_vertices_deg.iter_mut() {
        *v = (*v).to_radians();
    }
    let mut poly2_vertices_deg: [f64; 8] = [
        272.630446, -19.234210, 272.637274, -19.248542, 272.638942, -19.231476, 272.630868,
        -19.226364,
    ];
    for v in poly2_vertices_deg.iter_mut() {
        *v = (*v).to_radians();
    }
    let mut poly3_vertices_deg: [f64; 178] = [
        272.536719, -19.461249, 272.542612, -19.476380, 272.537389, -19.491509, 272.540192,
        -19.499823, 272.535455, -19.505218, 272.528024, -19.505216, 272.523437, -19.500298,
        272.514082, -19.503376, 272.502271, -19.500966, 272.488647, -19.490390, 272.481932,
        -19.490913, 272.476737, -19.486589, 272.487633, -19.455645, 272.500386, -19.444996,
        272.503003, -19.437557, 272.512303, -19.432436, 272.514132, -19.423973, 272.522103,
        -19.421523, 272.524511, -19.413250, 272.541021, -19.400024, 272.566264, -19.397500,
        272.564202, -19.389111, 272.569055, -19.383210, 272.588186, -19.386539, 272.593376,
        -19.381832, 272.596327, -19.370541, 272.624911, -19.358915, 272.629256, -19.347842,
        272.642277, -19.341020, 272.651322, -19.330424, 272.653174, -19.325079, 272.648903,
        -19.313708, 272.639616, -19.311098, 272.638128, -19.303083, 272.632705, -19.299839,
        272.627971, -19.289408, 272.628226, -19.276293, 272.633750, -19.270590, 272.615109,
        -19.241810, 272.614704, -19.221196, 272.618224, -19.215572, 272.630809, -19.209945,
        272.633540, -19.198681, 272.640711, -19.195292, 272.643028, -19.186751, 272.651477,
        -19.182729, 272.649821, -19.174859, 272.656782, -19.169272, 272.658933, -19.161883,
        272.678012, -19.159481, 272.689173, -19.176982, 272.689395, -19.183512, 272.678006,
        -19.204016, 272.671112, -19.206598, 272.664854, -19.203523, 272.662760, -19.211156,
        272.654435, -19.214434, 272.652969, -19.222085, 272.656724, -19.242136, 272.650071,
        -19.265092, 272.652868, -19.274296, 272.660871, -19.249462, 272.670041, -19.247807,
        272.675533, -19.254935, 272.673291, -19.273917, 272.668710, -19.279245, 272.671460,
        -19.287043, 272.667507, -19.293933, 272.669261, -19.300601, 272.663969, -19.307130,
        272.672626, -19.308954, 272.675225, -19.316490, 272.657188, -19.349105, 272.657638,
        -19.367455, 272.662447, -19.372035, 272.662232, -19.378566, 272.652479, -19.386871,
        272.645819, -19.387933, 272.642279, -19.398277, 272.629282, -19.402739, 272.621487,
        -19.398197, 272.611782, -19.405716, 272.603367, -19.404667, 272.586162, -19.422703,
        272.561792, -19.420008, 272.555815, -19.413012, 272.546500, -19.415611, 272.537427,
        -19.424213, 272.533081, -19.441402,
    ];
    for v in poly3_vertices_deg.iter_mut() {
        *v = (*v).to_radians();
    }
    let v1: Vec<(f64, f64)> = poly1_vertices_deg
        .iter()
        .copied()
        .step_by(2)
        .zip(poly1_vertices_deg.iter().copied().skip(1).step_by(2))
        .collect();
    let v2: Vec<(f64, f64)> = poly2_vertices_deg
        .iter()
        .copied()
        .step_by(2)
        .zip(poly2_vertices_deg.iter().copied().skip(1).step_by(2))
        .collect();
    let v3: Vec<(f64, f64)> = poly3_vertices_deg
        .iter()
        .copied()
        .step_by(2)
        .zip(poly3_vertices_deg.iter().copied().skip(1).step_by(2))
        .collect();
    let layer = get(18);
    let moc1 = layer.polygon_coverage(v1.as_slice(), true);
    let moc2 = layer.polygon_coverage(v2.as_slice(), true);
    let moc3 = layer.polygon_coverage(v3.as_slice(), true);

    let union12 = moc1.or(&moc2);
    let res1 = moc3.xor(&union12);
    let res2 = moc3.minus(&union12);
    let res3 = moc3.and(&union12.not());

    assert_eq!(res1, res2);
    assert_eq!(res1, res3);
}

#[test]
fn test_ok_bmoc_not() {
    let lon = 13.158329;
    let lat = -72.80028;
    let radius = 5.64323_f64.to_radians();
    let depth = 6;
    let delta_depth = 5;

    let bmoc = get(depth).cone_coverage_approx_custom(delta_depth, Degrees(lon, lat), radius);

    let not_bmoc = bmoc.not();

    println!("BMOC");
    for (flag, range) in bmoc.flagged_ranges() {
        println!("flag: {flag}; range: {range:?}");
    }
    println!("NOT BMOC");
    for (flag, range) in not_bmoc.flagged_ranges() {
        println!("flag: {flag}; range: {range:?}");
    }
    // Asssert that the first range has the flag 'true'.
    assert!(not_bmoc.flagged_ranges().next().unwrap().0);
}
