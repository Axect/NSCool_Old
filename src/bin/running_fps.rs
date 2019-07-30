#![feature(proc_macro_hygiene)]
extern crate peroxide;
extern crate inline_python;
extern crate natural_unit;
use peroxide::*;
use natural_unit::*;
use std::f64::consts::PI;
use inline_python::python;

pub const K: f64 = 8.2185f64;
pub const GAMMA: f64 = 1.7708f64;
pub const K2: f64 = 2.6709f64;
pub const GAMMA2: f64 = 1.4586f64;

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
        .set_step_size(1e-6*r_step)
        .set_times(2_000_000)
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state)
        .set_stop_condition(stop_by_crust);

    let results1 = tov_solver.integrate();

    // Integration
    tov_solver
        .set_stop_condition(stop_by_p);
    let results2 = tov_solver.integrate();
    println!("Integrate finish!");

    let value = results1.row(results1.row-1);
    println!("To crust: {} km", invert(value[0], Length, cgs_to_geom) / 100f64 / 1000f64);
    println!("To crust: {} solar mass", value[1]);

    let value2 = results2.row(results2.row-1);
    println!("To surface: {} km", invert(value2[0], Length, cgs_to_geom) / 100f64 / 1000f64);
    println!("To surface: {} solar mass", value2[1]);

    // Prepare vectors to input Python
    let result_r = concat(results1.col(0), results2.col(0));
    let result_m = concat(results1.col(1), results2.col(1));
    let result_rho = concat(results1.col(2), results2.col(2));
    let result_p = concat(
        results1.col(2).fmap(|x| K * x.powf(GAMMA)),
        results2.col(2).fmap(|x| K2 * x.powf(GAMMA2))
    );

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

pub fn tov_polytrope_crust(st: &mut State<f64>) {
    let r = st.param;
    let xs = &st.value; // m, rho, p
    let mut dx = &mut st.deriv;

    let m = xs[0];
    let rho = xs[1];
    let p = K2 * rho.powf(GAMMA2);

    let eps = rho + 1f64 / (GAMMA2 - 1f64) * K2 * rho.powf(GAMMA2);
    dx[0] = 4f64*PI*r.powi(2)*eps;
    dx[1] = -(eps + p) * (m + 4f64*PI*r.powi(3)*p) / (r * (r - 2f64*m)) * rho / (p * GAMMA2);
}

pub fn stop_by_p(st: &ExplicitODE) -> bool {
    let rho = st.get_state().value[1];
    K2*rho.powf(GAMMA2) < 1e-10
}

pub fn stop_by_crust(st: &ExplicitODE) -> bool {
    let rho = st.get_state().value[1];
    rho < 1e-6
}