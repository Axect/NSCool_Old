extern crate peroxide;
use peroxide::*;

#[allow(non_snake_case)]
use std::f64::consts::PI;
use NSCool::tov::eos::*;
use NSCool::consts::consts::*;

fn main() {
    let (rho_eos, p_eos, nbar_eos) = import_eos("data/APR_EOS_Acc_Fe.dat", UnitSystem::Geometrized);

    let solar = 2E+33;
    let electron = 9.109E-28;
    cgs_to_GeV(solar, Mass).print();
    cgs_to_MeV(electron, Mass).print(); // 0.51MeV
    cgs_to_natural(CGS.c, Velocity).print();
    cgs_to_natural(CGS.hbar, AngularMomentum).print();
}


//pub fn main() {
//    let p_c = K*RHO0C.powf(GAMMAF);
//    let init_val = c!(0, 0, p_c);
//    let mut records = solve_with_condition(tov_rhs, init_val, (0f64, 16f64), 1e-3, GL4(1e-15), |xs| xs[2] >= 0f64);
//    records.print();
//    records.write_with_header("data/tov_gl4.csv", vec!["r", "m", "phi", "p"]).expect("Can't write file");
//}

