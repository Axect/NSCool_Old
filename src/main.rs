#![feature(proc_macro_hygiene)]
extern crate peroxide;
//extern crate inline_python;
extern crate natural_unit;
use peroxide::*;
use natural_unit::*;
use std::f64::consts::PI;
use NSCool::tov::eos::*;
//use inline_python::python;

pub const K: f64 = 100f64;
pub const GAMMA: f64 = 2f64;

#[allow(non_snake_case)]
fn main() {
    // Define new unit (c=G=M=1)
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;
    let cgs_to_geom = ConversionFactor::new(
        1f64 / M,
        c.powi(2) / (G*M),
        c.powi(3) / (G*M),
    );

    // Set Initial Conditions
    let rho_c = convert(3.48E+15f64, Density, cgs_to_geom);
    let r_step = convert(10_000_00f64, Length, cgs_to_geom);
    let m_c = 0f64;
    let init_state = State::<f64>::new(0f64, vec![m_c, rho_c], vec![0f64; 2]);

    // Insert ODE
    let mut tov_solver = ExplicitODE::new(tov_polytrope);
    tov_solver
        .set_step_size(1e-5*r_step)
        .set_times(200_000)
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state)
        .set_stop_condition(stop_by_p);

    // Integration
    let results = tov_solver.integrate();
    println!("Integrate finish");

    // Prepare vectors to input Python
    let result_r = results.col(0);
    let result_m = results.col(1);
    let result_rho = results.col(2);
    let result_p = result_rho.fmap(|x| K * x.powf(GAMMA));

//    // Python plot code
//    python! {
//        import pylab as plt
//        import numpy as np
//
//        # Use latex
//        plt.rc("text", usetex=True)
//        plt.rc("font", family="serif")
//
//        # Plot
//        plt.figure(figsize=(10,6), dpi=300)
//        plt.title(r"$R$ vs $\rho$", fontsize=16)
//        plt.xlabel(r"$r$", fontsize=14)
//        plt.ylabel(r"$\rho$", fontsize=14)
//        plt.plot('result_r, 'result_rho)
//        plt.grid()
//        plt.savefig("data/r_vs_rho.png")
//
//        # Plot2
//        plt.figure(figsize=(10,6), dpi=300)
//        plt.title(r"$R$ vs $m$", fontsize=16)
//        plt.xlabel(r"$r$", fontsize=14)
//        plt.ylabel(r"$m$", fontsize=14)
//        plt.plot('result_r, 'result_m)
//        plt.grid()
//        plt.savefig("data/r_vs_m.png")
//
//        # Plot3
//        plt.figure(figsize=(10,6), dpi=300)
//        plt.title(r"$R$ vs $P$", fontsize=16)
//        plt.xlabel(r"$r$", fontsize=14)
//        plt.ylabel(r"$p$", fontsize=14)
//        plt.plot('result_r, 'result_p)
//        plt.grid()
//        plt.savefig("data/r_vs_p.png")
//    }
}

pub fn tov_polytrope(st: &mut State<f64>) {
    let r = st.param;
    let xs = &st.value; // m, rho, p
    let mut dx = &mut st.deriv;

    let m = xs[0];
    let rho = xs[1];
    let p = K * rho.powf(GAMMA);

    let eps = rho + 1f64 / (GAMMA - 1f64) * K * rho.powf(GAMMA);
    dx[0] = 4f64*PI*r.powi(2)*eps;
    if r < 1e-4 {
        let m_r3 = 4f64 * PI / 3f64 * eps;
        dx[1] = -(eps + p) * (m_r3 * r + 4f64*PI*r*p) / (1f64 - 2f64*m_r3*r.powi(2)) * rho / (p * GAMMA);
    } else {
        dx[1] = -(eps + p) * (m + 4f64*PI*r.powi(3)*p) / (r * (r - 2f64*m)) * rho / (p * GAMMA);
    }
}

pub fn stop_by_p(st: &ExplicitODE) -> bool {
    let rho = st.get_state().value[1];
    (K*rho.powf(GAMMA) - 0f64).abs() < 1e-6
}