extern crate peroxide;
use peroxide::*;

#[allow(non_snake_case)]
use std::f64::consts::PI;
use NSCool::tov::eos::*;
use NSCool::consts::consts::*;

fn main() {
//    let (rho_eos, p_eos, nbar_eos) = import_eos("data/APR_EOS_Acc_Fe.dat", UnitSystem::CGS);
//
//    rho_eos.print();
    (1f64 / 1.122e6).print();
    (CGS.c.powi(2) * cgs_to_geom(1f64, Density)).print();
}

//
//pub fn main() {
//    let p_c = K*RHO0C.wpowf(GAMMAF);
//    let init_val = c!(0, 0, p_c);
//    let mut records = solve_with_condition(tov_rhs, init_val, (0f64, 16f64), 1e-3, RK4, |xs| xs[2] >= 0f64);
//    records.print();
//    records.write_with_header("data/tov_rk4.csv", vec!["r", "m", "phi", "p"]).expect("Can't write file");
//}

