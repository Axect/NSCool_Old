#![feature(proc_macro_hygiene)]
extern crate peroxide;
extern crate inline_python;
extern crate natural_unit;
use peroxide::*;
use natural_unit::*;
use std::f64::consts::PI;
use inline_python::python;

pub const K0: f64 = 0.3715;
pub const GAMMA0: f64 = 1.3399;
pub const GAMMA1: f64 = 0.7224;
pub const GAMMA2: f64 = 1.2221;
pub const GAMMA3: f64 = 2.7666;
pub const GAMMA4: f64 = 1.1997;
pub const RHO0: f64 = -6.4f64;
pub const RHO1: f64 = -5.55f64;
pub const RHO2: f64 = -3.6f64;
pub const RHO3: f64 = -2.4f64;

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
    let rho_c = convert(0.1057380000E+18, Density, cgs_to_geom);
    let r_step = convert(10_000_00f64, Length, cgs_to_geom);
    let m_c = 0f64;
    let init_state = State::<f64>::new(0f64, vec![m_c, rho_c], vec![0f64; 2]);

    // Insert ODE
    let mut tov_solver = ExplicitODE::new(tov_piecewise_polytrope);
    tov_solver
        .set_step_size(1e-6*r_step)
        .set_times(2_000_000)
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state)
        .set_stop_condition(stop_by_p);

    let results1 = tov_solver.integrate();
    println!("Integrate finish!");

    let value = results1.row(results1.row-1);
    println!("Radius: {} km", invert(value[0], Length, cgs_to_geom) / 100f64 / 1000f64);
    println!("Mass: {} solar mass", value[1]);

    // Prepare vectors to input Python
    let result_r = results1.col(0);
    let result_m = results1.col(1);
    let result_rho = results1.col(2);
    let result_p = results1.col(2).fmap(|rho| {
        let k1 = K0 * 10f64.powf(RHO0).powf(GAMMA0) / 10f64.powf(RHO0).powf(GAMMA1);
        let k2 = k1 * 10f64.powf(RHO1).powf(GAMMA1) / 10f64.powf(RHO1).powf(GAMMA2);
        let k3 = k2 * 10f64.powf(RHO2).powf(GAMMA2) / 10f64.powf(RHO2).powf(GAMMA3);
        let k4 = k3 * 10f64.powf(RHO3).powf(GAMMA3) / 10f64.powf(RHO3).powf(GAMMA4);
        let k: f64;
        let g: f64;

        if rho > 10f64.powf(RHO3) {
            k = k4;
            g = GAMMA4;
        } else if rho > 10f64.powf(RHO2) {
            k = k3;
            g = GAMMA3;
        } else if rho > 10f64.powf(RHO1) {
            k = k2;
            g = GAMMA2;
        } else if rho > 10f64.powf(RHO0) {
            k = k1;
            g = GAMMA1;
        } else {
            k = K0;
            g = GAMMA0;
        }
        k * rho.powf(g)
    });

    // Python plot code
    python! {
        import pylab as plt
        import numpy as np

        # Use latex
        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")

        # Plot
        plt.figure(figsize=(10,6), dpi=300)
        plt.title(r"$R$ vs $\rho$", fontsize=16)
        plt.xlabel(r"$r$", fontsize=14)
        plt.ylabel(r"$\rho$", fontsize=14)
        plt.plot('result_r, 'result_rho)
        plt.grid()
        plt.savefig("figure/FPS/r_vs_rho.png")

        # Plot2
        plt.figure(figsize=(10,6), dpi=300)
        plt.title(r"$R$ vs $m$", fontsize=16)
        plt.xlabel(r"$r$", fontsize=14)
        plt.ylabel(r"$m$", fontsize=14)
        plt.plot('result_r, 'result_m)
        plt.grid()
        plt.savefig("figure/FPS/r_vs_m.png")

        # Plot3
        plt.figure(figsize=(10,6), dpi=300)
        plt.title(r"$R$ vs $P$", fontsize=16)
        plt.xlabel(r"$r$", fontsize=14)
        plt.ylabel(r"$p$", fontsize=14)
        plt.plot('result_r, 'result_p)
        plt.grid()
        plt.savefig("figure/FPS/r_vs_p.png")
    }
}

pub fn tov_piecewise_polytrope(st: &mut State<f64>) {
    let r = st.param;
    let xs = &st.value; // m, rho, p
    let mut dx = &mut st.deriv;

    let m = xs[0];
    let rho = xs[1];


    let k1 = K0 * 10f64.powf(RHO0).powf(GAMMA0) / 10f64.powf(RHO0).powf(GAMMA1);
    let k2 = k1 * 10f64.powf(RHO1).powf(GAMMA1) / 10f64.powf(RHO1).powf(GAMMA2);
    let k3 = k2 * 10f64.powf(RHO2).powf(GAMMA2) / 10f64.powf(RHO2).powf(GAMMA3);
    let k4 = k3 * 10f64.powf(RHO3).powf(GAMMA3) / 10f64.powf(RHO3).powf(GAMMA4);
    let k: f64;
    let g: f64;

    if rho > 10f64.powf(RHO3) {
        k = k4;
        g = GAMMA4;
    } else if rho > 10f64.powf(RHO2) {
        k = k3;
        g = GAMMA3;
    } else if rho > 10f64.powf(RHO1) {
        k = k2;
        g = GAMMA2;
    } else if rho > 10f64.powf(RHO0) {
        k = k1;
        g = GAMMA1;
    } else {
        k = K0;
        g = GAMMA0;
    }

    let p = k * rho.powf(g);

    let eps = rho + 1f64 / (g - 1f64) * k * rho.powf(g);
    dx[0] = 4f64*PI*r.powi(2)*eps;
    if r < 1e-5 {
        let m_r3 = 4f64 * PI / 3f64 * eps;
        dx[1] = -(eps + p) * (m_r3 * r + 4f64*PI*r*p) / (1f64 - 2f64*m_r3*r.powi(2)) * rho / (p * g);
    } else {
        dx[1] = -(eps + p) * (m + 4f64*PI*r.powi(3)*p) / (r * (r - 2f64*m)) * rho / (p * g);
    }
}

pub fn stop_by_p(st: &ExplicitODE) -> bool {
    let rho = st.get_state().value[1];
    K0 * rho.powf(GAMMA0) < 1e-11
}