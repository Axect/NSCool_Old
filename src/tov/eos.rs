extern crate peroxide;
use peroxide::*;

use std::f64::consts::PI;

pub const K: f64 = 100.;
pub const GAMMA: usize = 2;
pub const GAMMAF: f64 = 2.;
pub const RHO0C: f64 = 1.28e-3;

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
pub fn tov_rhs(rs: Vec<f64>) -> Vec<f64> {
    let r = rs[0];
    let r2 = r.powi(2);

    let m_old = rs[1];
    let _phi_old = rs[2];
    let p_old = rs[3];

    let rho_old = rho_of_p(p_old);

    let m_r: f64;
    let m_rr: f64;
    let m_rrr: f64;

    // Using Taylor Series m = 1/6 m'''(0) r^3
    if r < 1e-3 {
        m_rrr = 4.0/3.0 * PI * (RHO0C.powi(2)/(GAMMAF - 1.)
                                + 4f64.powf(-1./GAMMAF) * (RHO0C).powf(2./GAMMAF));
        m_rr = m_rrr * r;
        m_r = m_rr * r;
    } else {
        m_r = m_old / r;
        m_rr = m_r / r;
        m_rrr = m_rr / r;
    }

    let m_new = 4f64 * PI * r2 * rho_old;
    let p_new = - rho_old * m_rr * (1.0 + p_old/rho_old) 
        * (1.0 + 4.0*PI*p_old / m_rrr)
        * 1.0/(1.0 - 2.0*m_r);
    let phi_new = - 1.0/rho_old * p_new / (1.0 + p_old/rho_old);

    vec![m_new, phi_new, p_new]
}

pub fn rho_of_p(p: f64) -> f64 {
    let mut p_result = p;
    
    if p < 0.0 {
        p_result = p.abs();
    }

    (p_result/K).powf(1.0/GAMMAF) + p_result / (GAMMAF - 1.0)
}
