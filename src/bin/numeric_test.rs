extern crate peroxide;
use peroxide::*;
use std::f64::consts::PI;

fn main() {
    let K = 100f64;
    let GAMMA = 2f64;
    let r = 1e-3f64;

    let rho0_c = 1.28e-3f64;

    let p_c = K * rho0_c.powf(GAMMA);

    let rho_0 = (p_c / K).powf(1. / GAMMA);
    let rho_old = rho_0 + p_c / (GAMMA - 1.0);
    let m_old = 0f64;

    let dm = 4f64 * PI * r.powi(2) * rho_old;
    let dphi = (m_old + 4f64 * PI * r.powi(3) * p_c) / (r.powi(2) * (1f64 - 2f64 * m_old / r));
    let dp = - (rho_old + p_c) * (m_old + 4f64 * PI * r.powi(3) * p_c) / (r.powi(2) * (1f64 - 2f64 * m_old / r));

    dm.print();
    dphi.print();
    dp.print();

    (1f64 - 2f64 * m_old / r).print();
}