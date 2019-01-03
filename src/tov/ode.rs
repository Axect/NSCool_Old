extern crate peroxide;
use peroxide::*;

#[derive(Debug)]
pub struct ODE {
    pub param: f64,
    pub values: Vec<f64>,
    pub step: f64,
    pub records: Matrix,
    pub stage: usize,
    pub num: usize,
}
