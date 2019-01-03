extern crate peroxide;

use peroxide::*;

#[allow(non_snake_case)]
use NSCool::tov::ode::*;
use NSCool::tov::rk4::*;

use std::f64::consts::PI;

const K: f64 = 100f64;
const GAMMA: usize = 2;
const GAMMAF: f64 = 2f64;

#[allow(unused_must_use_variable)]
pub fn main() {
    let init_val = c!(1e-16, 0, 0, 1.28e-3);
    let results = bdf1(init_val, tov_derivative_dual, 1e-3, 1e-15, 16_000);
    results.write("data/tov_dual.csv");
}

// TOV
//fn main() {
//    let mut tov_solver = ODE::new(
//        1e-16,
//        c!(0, 0, 1.28*1e-3),
//        1e-3,
//        16 * 1000,
//    );
//
//    tov_solver.integrate(|r, rs| tov_derivative(r, rs));
//
//    tov_solver.param.print();
//    tov_solver.values.print();
//
//    tov_solver.records.write("tov.csv");
//}

//fn tov_derivative(r: f64, rs: Vec<f64>) -> Vec<f64> {
//    let r2 = r.powi(2);
//    let r3 = r.powi(3);
//
//    let m_old = rs[0];
//    let phi_old = rs[1];
//    let rho_old = rs[2];
//
//    let rho2 = rho_old.powi(2);
//    let rho_gamma = rho_old.powi(GAMMA);
//
//    let m_new = 4f64*PI*r2*rho_old;
//    let phi_new = phi_old
//        + (4f64*PI*K*r3 * rho_gamma + m_old)
//        / (r2 - 2f64*r*m_old);
//    let rho_new = -(4f64*PI*K*K*r3 * rho_old.powi(GAMMA+1)
//        + m_old
//        + (4f64*PI*K*r3*rho2 + K*m_old*rho_old))
//        / ((GAMMAF*K*r2 - 2f64*GAMMAF*K*r*m_old));
//
//    vec![m_new, phi_new, rho_new]
//}

fn tov_derivative_dual(rs: Vec<Dual>) -> Vec<Dual> {
    let r = rs[0];
    let r2 = r.pow(2);
    let r3 = r.pow(3);

    let m_old = rs[1];
    let phi_old = rs[2];
    let rho_old = rs[3];

    let rho2 = rho_old.pow(2);
    let rho_gamma = rho_old.pow(GAMMA);

    let m_new = 4f64 * PI * r2 * rho_old;
    let phi_new = phi_old + (4f64 * PI * K * r3 * rho_gamma + m_old) / (r2 - 2f64 * r * m_old);
    let rho_new = -(4f64 * PI * K * K * r3 * rho_old.pow(GAMMA + 1)
        + m_old
        + (4f64 * PI * K * r3 * rho2 + K * m_old * rho_old))
        / (GAMMAF * K * r2 - 2f64 * GAMMAF * K * r * m_old);

    vec![m_new, phi_new, rho_new]
}
