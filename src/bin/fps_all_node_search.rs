extern crate natural_unit;
extern crate peroxide;
extern crate rayon;
use natural_unit::*;
use peroxide::*;
use rayon::prelude::*;
use std::f64::consts::PI;

#[allow(non_snake_case)]
fn main() {
    // Define new unit (c=G=M=1)
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;
    let cgs_to_geom = ConversionFactor::new(1f64 / M, c.powi(2) / (G * M), c.powi(3) / (G * M));

    let fps_data = Matrix::read("data/fps.dat", false, ' ').expect("Can't read APR");

    let log_rho = fps_data
        .col(0)
        .fmap(|x| convert(x, Density, cgs_to_geom.clone()).log10());
    let log_p = fps_data
        .col(1)
        .fmap(|x| convert(x, Pressure, cgs_to_geom.clone()).log10());

    log_rho.len().print();
    log_p.len().print();

    let kgs = vec![100f64, 2f64, 2f64, 2f64, 2f64, 2f64];

    let data = hstack!(log_rho.clone(), log_p.clone());

    // let mut opt = Optimizer::new(data, piecewise_polytrope);
    // let param = opt.set_init_param(kgs.clone())
    //     .set_method(LevenbergMarquardt)
    //     .set_max_iter(100)
    //     .optimize();
    // param.print();
}

#[allow(non_snake_case)]
fn piecewise_polytrope(
    log_rho: &Vec<f64>,
    kr: Vec<Number>,
    RHO_interval: &Vec<f64>,
) -> Vec<Number> {
    let k0 = kr[0];
    let g0 = kr[1];
    let g1 = kr[2];
    let g2 = kr[3];
    let g3 = kr[4];
    let g4 = kr[5];
    let rho0 = RHO_interval[0];
    let rho1 = RHO_interval[1];
    let rho2 = RHO_interval[2];
    let rho3 = RHO_interval[3];
    let k1 = (k0 * Number::F(10f64.powf(rho0)).powf(g0)) / Number::F(10f64.powf(rho0)).powf(g1);
    let k2 = (k1 * Number::F(10f64.powf(rho1)).powf(g1)) / Number::F(10f64.powf(rho1)).powf(g2);
    let k3 = (k2 * Number::F(10f64.powf(rho2)).powf(g2)) / Number::F(10f64.powf(rho2)).powf(g3);
    let k4 = (k3 * Number::F(10f64.powf(rho3)).powf(g3)) / Number::F(10f64.powf(rho3)).powf(g4);
    log_rho
        .clone()
        .into_iter()
        .map(|x| {
            if x <= rho0 {
                k0.log10() + g0 * x
            } else if x <= rho1 {
                k1.log10() + g1 * x
            } else if x <= rho2 {
                k2.log10() + g2 * x
            } else if x <= rho3 {
                k3.log10() + g3 * x
            } else {
                k4.log10() + g4 * x
            }
        })
        .collect()
}
