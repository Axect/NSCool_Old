extern crate peroxide;
extern crate natural_unit;
use natural_unit::*;

use std::f64::consts::PI;
pub use self::EOSModel::*;
use peroxide::{State, dual, Real};
use self::peroxide::operation::number::Number::D;
use self::peroxide::operation::number::Number;
use self::peroxide::{PowOps, ExpLogOps, seq};

#[derive(Debug, Clone, Copy)]
pub enum EOSModel {
    APR,
    BSk19,
    BSk20,
    BSk21,
    BSk22,
    BSk24,
    BSk25,
    BSk26,
}

/// Pressure as a function of density
///
/// # Variables
///
/// * input: log(rho / g cm^{-3})
/// * output: log(P / dyn cm^{-2})
///
/// # Caution
///
/// * Analytical representation of Potekhin is slightly different with Pearson
///     * Peroxide follows Potekhin's notation (p19 <-> p20 & p22 <-> p23 of Pearson et al.)
///
/// # Reference
///
/// [Potekhin et al.](https://arxiv.org/abs/1310.0049)
/// [Pearson et al.](https://arxiv.org/abs/1903.04981)
pub fn log_p_of_log_rho(xi: Number, model: EOSModel) -> Number {
    let a: [f64; 23] = match model {
        APR => PARAM_APR,
        BSk19 => PARAM_BSk19,
        BSk20 => PARAM_BSk20,
        BSk21 => PARAM_BSk21,
        BSk22 => PARAM_BSk22,
        BSk24 => PARAM_BSk24,
        BSk25 => PARAM_BSk25,
        BSk26 => PARAM_BSk26,
    };

    (a[0] + a[1]*xi + a[2]*xi.powi(3)) / ((1f64 + a[3]*xi) * ((a[4] * (xi - a[5])).exp() + 1f64))
    + (a[6] + a[7]*xi) / ((a[8] * (a[5] - xi)).exp() + 1f64)
    + (a[9] + a[10]*xi) / ((a[11] * (a[12] - xi)).exp() + 1f64)
    + (a[13] + a[14]*xi) / ((a[15] * (a[16] - xi)).exp() + 1f64)
    + a[17] / (1f64 + (a[18] * (xi - a[19])).powi(2))
    + a[20] / (1f64 + ((xi - a[22]) * a[21]).powi(2)) // MeV fm-3
}

pub fn polytropes(rho_0: Number, K: f64) -> Number {
    K * rho_0.powi(2)
}

pub const PARAM_APR: [f64; 23] = [
    6.22,       // 1
    6.121,      // 2
    0.006035,   // 3
    0.16354,    // 4
    4.73,       // 5
    11.5765,    // 6
    12.589,     // 7
    1.4365,     // 8
    4.75,       // 9
    -42.489,    // 10
    3.8175,     // 11
    2.3,        // 12
    14.81,      // 13
    29.80,      // 14
    -2.976,     // 15
    1.99,       // 16
    14.93,      // 17
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
];

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

pub const PARAM_BSk22: [f64; 23] = [
    6.682,
    5.651,
    0.00459,
    0.14359,
    2.681,
    11.972,
    13.993,
    1.2904,
    2.665,
    -27.787,
    2.014,
    4.09,
    14.135,
    28.03,
    -1.921,
    1.08,
    14.89,
    0.098,
    4.75,
    11.67,
    -0.037,
    11.9,
    14.1,
];

pub const PARAM_BSk24: [f64; 23] = [
    6.795,
    5.552,
    0.00435,
    0.13963,
    3.636,
    11.943,
    13.848,
    1.3031,
    3.644,
    -30.840,
    2.2322,
    4.65,
    14.29,
    30.08,
    -2.080,
    1.1,
    14.71,
    0.099,
    5.,
    11.66,
    -0.095,
    9.1,
    14.15,
];

pub const PARAM_BSk25: [f64; 23] = [
    7.21,
    5.196,
    0.00328,
    0.12516,
    4.624,
    12.16,
    9.348,
    1.6624,
    4.66,
    -28.232,
    2.0638,
    5.27,
    14.365,
    29.1,
    -2.130,
    0.865,
    14.66,
    0.069,
    6.3,
    11.65,
    -0.172,
    8.6,
    14.18,
];

pub const PARAM_BSk26: [f64; 23] = [
    3.672,
    7.844,
    0.00876,
    0.22604,
    3.129,
    11.939,
    13.738,
    1.3389,
    3.112,
    -23.031,
    1.6264,
    4.83,
    14.272,
    23.28,
    -1.542,
    2.1,
    15.31,
    0.083,
    6.16,
    11.66,
    -0.042,
    14.8,
    14.18,
];