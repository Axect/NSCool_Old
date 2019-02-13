extern crate peroxide;
use peroxide::*;

#[allow(non_snake_case)]
use NSCool::tov::eos::*;

use std::f64::consts::PI;


pub fn main() {
    let p_c = K*RHO0C.powf(GAMMAF);
    let init_val = c!(0, 0, p_c);
    let mut records = solve_with_condition(tov_rhs, init_val, (0f64, 16f64), 1e-3, RK4, |xs| xs[2] >= 0f64);
    //let mut records = solve(tov_rhs, init_val, (0f64, 16f64), 1e-3, RK4);
    records.write_with_header("data/tov_rk4.csv", vec!["r", "m", "phi", "p"]).expect("Can't write file");
}

fn tov_derivative_dual(r: Dual, rs: Vec<Dual>) -> Vec<Dual> {
    let r2 = r.pow(2);
    let r3 = r.pow(3);

    let m_old = rs[0];
    let phi_old = rs[1];
    let rho_old = rs[2];

    let rho2 = rho_old.pow(2);
    let rho_gamma = rho_old.pow(GAMMA);

    let m_new = 4f64 * PI * r2 * rho_old;
    let phi_new = (4f64 * PI * K * r3 * rho_gamma + m_old) / (r2 - 2f64 * r * m_old);
    let rho_new = -(4f64 * PI * K * K * r3 * rho_old.pow(GAMMA + 1)
        + m_old
        + (4f64 * PI * K * r3 * rho2 + K * m_old * rho_old))
        / (GAMMAF * K * r2 - 2f64 * GAMMAF * K * r * m_old);

    vec![m_new, phi_new, rho_new]
}
