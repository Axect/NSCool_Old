extern crate peroxide;
use peroxide::*;

use std::f64::consts::PI;

const K: f64 = 100.;
const GAMMA: usize = 2;
const GAMMAF: f64 = 2.;
const RHO0C: f64 = 1.28e-3;

/// Tolman-Oppenheimer-Volkoff Equations for Spherically Symmetric Eqaulibrium Stars
///
/// # Equation
/// ```latex
/// dm/dr = 4πr^2ρ
/// dP/dr = -ρm/r^2 (1 + P/ρ) (1 + 4πPr^3/m) (1 - 2m/r)^{-1}
/// dΦ/dr = -1/ρ dP/dr (1 + P/ρ)^{-1}
/// P = Kρ_0^Γ
/// Γ = 1 + 1/n
/// ρ = ρ_0(1 + ε)
/// m = 0, P = P_c at r = 0
/// ```
pub fn tov_rhs(rs: Vec<Dual>) -> Vec<Dual> {
    let r = rs[0];
    let r2 = r.pow(2);

    let m_old = rs[1];
    let _phi_old = rs[2];
    let p_old = rs[3];

    let rho_old = rho_dual(p_old.clone());

    let m_r: f64;
    let m_rr: f64;
    let m_rrr: f64;

    // Using Taylor Series m = 1/6 m'''(0) r^3
    if r.value() < 1e-3 {
        m_rrr = 4.0/3.0 * PI * (RHO0C.powi(2)/(GAMMAF - 1.)
                                + 4f64.powf(-1./GAMMAF) * (RHO0C).powf(2./GAMMAF));
        m_rr = m_rrr * r.value();
        m_r = m_rr * r.value();
    } else {
        m_r = m_old.value() / r.value();
        m_rr = m_r / r.value();
        m_rrr = m_rr / r.value();
    }

    let m_new = 4f64 * PI * r2 * rho_old;
    let p_new = - rho_old * m_rr * (1.0 + p_old/rho_old) * (1.0 + 4.0*PI*p_old / m_rrr)
        * 1.0/(1.0 - 2.0*m_r);
    let phi_new = - 1.0/rho_old * p_new / (1.0 + p_old/rho_old);

    vec![m_new, phi_new, p_new]
}

pub fn rho_dual(p: Dual) -> Dual {
    let mut p_result = p.clone();
    
    if p.value() < 0.0 {
        p_result = Dual::new(p.value().abs(), p.slope().abs());
    }

    (p_result/K).powf(1.0/GAMMAF) + p_result / (GAMMAF - 1.0)
}
