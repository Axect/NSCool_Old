extern crate peroxide;
use peroxide::*;

use std::f64::consts::PI;
use crate::consts::consts::*;
pub use self::EOSModel::*;

// ===========================================================
// Declare Constants
// ===========================================================
pub const K: f64 = 300.;
pub const GAMMA: usize = 2;
pub const GAMMAF: f64 = 2.;
pub const RHO0C: f64 = 0.4;

#[derive(Debug, Clone)]
pub enum EOSModel {
    BSk19,
    BSk20,
    BSk21
}

/// Pressure as a function of density
///
/// # Variables
///
/// * input: log(rho / g cm^{-3})
/// * output: log(P / dyn cm^{-2})
///
/// # Reference
///
/// [Potekhin et al.](https://arxiv.org/pdf/1310.0049.pdf)
pub fn log_p_of_log_rho(xi: f64, model: EOSModel) -> f64 {
    let a: [f64; 23] = match model {
        BSk19 => PARAM_BSk19,
        BSk20 => PARAM_BSk20,
        BSk21 => PARAM_BSk21
    };

    (a[0] + a[1]*xi + a[2]*xi.powi(3)) / (1f64 + a[3]*xi) / ((a[4] * (xi - a[5])).exp() + 1f64)
    + (a[6] + a[7]*xi) / ((a[8] * (a[5] - xi)).exp() + 1f64)
    + (a[9] + a[10]*xi) / ((a[11] * (a[12] - xi)).exp() + 1f64)
    + (a[13] + a[14]*xi) / ((a[15] * (a[16] - xi)).exp() + 1f64)
    + a[17] / (1f64 + (a[18] * (xi - a[19])).powi(2))
    + a[20] / (1f64 + (a[21] * (xi - a[22])).powi(2))
}



/// Tolman-Oppenheimer-Volkoff Equations for Spherically Symmetric Equilibrium Stars
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
pub fn tov_rhs(r: Dual, rs: Vec<Dual>) -> Vec<Dual> {
    let m_old = rs[0];
    let _phi_old = rs[1];
    let p_old = rs[2];

    let rho_0 = (p_old / K).powf(1. / GAMMAF);
    let rho_old = rho_0 + p_old / (GAMMAF - 1f64);

    let dm = 4f64 * PI * r.pow(2) * rho_old;

    let dphi = if r.value() != 0f64 { (m_old + 4f64 * PI * r.pow(3) * p_old) / (r.pow(2) * (1f64 - 2f64 * m_old / r)) } else { dual(0, 0) };
    let dp = if r.value() != 0f64 {- (rho_old + p_old) * (m_old + 4f64 * PI * r.pow(3) * p_old) / (r.pow(2) * (1f64 - 2f64 * m_old / r)) } else { dual(0, 0) };

    vec![dm, dphi, dp]
}

/// Import eos
///
/// # Column (Each units are CGS)
///
///   Rho        Press         nbar        Ye         Ymu        Yn
/// [g/cm3]   [dyne/cm2]     [#/fm3]     [A_cell]    [A_ion]     [Z]
pub fn import_eos(title: &str, unit: UnitSystem) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let m = Matrix::read(title, false, ' ').expect("Can't read file");
    let rho_eos = m.col(0);
    let p_eos = m.col(1);
    let nbar_eos = m.col(2);

    match unit {
        UnitSystem::CGS => {
            (rho_eos, p_eos, nbar_eos)
        },
        UnitSystem::Geometrized => {
            (
                rho_eos.fmap(|t| cgs_to_geom(t, Density)),
                p_eos.fmap(|t| cgs_to_geom(t, Pressure)),
                nbar_eos
            )
        },
        UnitSystem::Natural => {
            (
                rho_eos.fmap(|t| cgs_to_natural(t, Density)),
                p_eos.fmap(|t| cgs_to_natural(t, Pressure)),
                nbar_eos
            )
        }
    }

}

pub const PARAM_BSk19: [f64; 23] = [
    3.916,
    7.701,
    0.00858,
    0.22114,
    3.269,
    11.964,
    13.349,
    1.3683,
    3.254,
    -12.953,
    0.9237,
    6.20,
    14.383,
    16.693,
    -1.0514,
    2.486,
    15.362,
    0.085,
    6.23,
    11.68,
    -0.029,
    20.1,
    14.19
];

pub const PARAM_BSk20: [f64; 23] = [
    4.078,
    7.587,
    0.00839,
    0.21695,
    3.614,
    11.942,
    13.751,
    1.3373,
    3.606,
    -22.996,
    1.6229,
    4.88,
    14.274,
    23.560,
    -1.5564,
    2.095,
    15.294,
    0.084,
    6.36,
    11.67,
    -0.042,
    14.8,
    14.18
];

pub const PARAM_BSk21: [f64; 23] = [
    4.857,
    6.981,
    0.00706,
    0.19351,
    4.085,
    12.065,
    10.521,
    1.5905,
    4.104,
    -28.726,
    2.0845,
    4.89,
    14.302,
    22.881,
    -1.7690,
    0.989,
    15.313,
    0.091,
    4.68,
    11.65,
    -0.086,
    10.0,
    14.15
];