extern crate peroxide;

#[allow(unused_imports)]
use peroxide::*;
use crate::structure::eos::PiecewisePolytrope;

#[derive(Debug, Clone)]
pub struct NSStructure {
    pub r: Vec<f64>,
    pub rho: Vec<f64>,
    pub p: Vec<f64>,
    pub n: Vec<f64>,
    pub m: f64,
}

pub fn tov_piecewise_polytrope(st: &mut State<f64>, p: &PiecewisePolytrope) {
    let r = st.param;
    let xs = &st.value;             // m, rho, p
    let mut dx = &mut st.deriv;

    let m = xs[0];
    let rho = xs[1];

    let kgs = p.extract_k_gamma();
    let rhos = p.get_interval();

    let mut j = 0usize;

    unimplemented!()
}
