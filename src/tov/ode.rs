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

#[allow(non_snake_case)]
pub fn jacobian<F>(x: Vec<f64>, f: F) -> Matrix
    where F: Fn(Vec<Dual>) -> Vec<Dual>
{
    let x_var = x.clone().into_iter()
        .map(|t| Dual::new(t, 1.))
        .collect::<Vec<Dual>>();

    let x_const = x.clone().into_iter()
        .map(|t| Dual::new(t, 0.))
        .collect::<Vec<Dual>>();

    let l = x.len();

    let mut J: Matrix = zeros!(l, l);

    let mut x_temp = x_const.clone();

    for i in 0 .. l {
        x_temp[i] = x_var[i];
        let dual_temp = f(x_temp.clone());
        let slope_temp = dual_temp.into_iter()
            .map(|dx| dx.slope())
            .collect::<Vec<f64>>();
        for j in 0 .. l {
            J[(j, i)] = slope_temp[j];
        }
        x_temp = x_const.clone();
    }
    J
}