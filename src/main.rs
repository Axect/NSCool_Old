extern crate peroxide;
extern crate natural_unit;
use peroxide::*;
use natural_unit::*;

#[allow(non_snake_case)]
use std::f64::consts::PI;
use NSCool::tov::eos::*;

pub const K: f64 = 30000f64;
pub const Gamma: f64 = 2f64;

fn main() {
    // Define new unit
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;

    // c=G=M=1
    let cgs_to_geom = ConversionFactor::new(
        1f64 / M,
        c.powi(2) / (G*M),
        c.powi(3) / (G*M),
    );

    let rho_c = convert(3.48E+15f64, Density, cgs_to_geom);
    let r_step = convert(10_000_00f64, Length, cgs_to_geom);
    let m_c = 0f64;

    let solar_rad = convert(CONSTANT_CGS.r_solar, Length, cgs_to_geom);

    let init_state = State::<f64>::new(0f64, vec![m_c, rho_c], vec![0f64; 2]);

    let mut tov_solver = ExplicitODE::new(tov_polytrope);
    tov_solver
        .set_step_size(1e-5*r_step)
        .set_times(100_000)
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state);

    solar_rad.print();
    r_step.print();

    let results = tov_solver.integrate();

    let l = results.row;
    let final_m = results.row(l-1)[1];
    final_m.print();

    let mut simple_writer = SimpleWriter::new();
    simple_writer
        .set_path("data/tov_BSk19_RK4.pickle")
        .insert_matrix(results)
        .write_pickle();
}

pub fn tov_polytrope(st: &mut State<f64>) {
    let r = st.param;
    let xs = &st.value; // m, rho, p
    let mut dx = &mut st.deriv;

    let m = xs[0];
    let rho = xs[1];
    let p = K * rho.powf(Gamma);

    let eps = rho + 1f64 / (Gamma - 1f64) * K * rho.powf(Gamma);
    dx[0] = 4f64*PI*r.powi(2)*eps;
    if r < 1e-4 {
        let m_r3 = 4f64 * PI / 3f64 * eps;
        dx[1] = -(eps + p) * (m_r3 * r + 4f64*PI*r*p) / (1f64 - 2f64*m_r3*r.powi(2)) * rho / (p * Gamma);
    } else {
        dx[1] = -(eps + p) * (m + 4f64*PI*r.powi(3)*p) / (r * (r - 2f64*m)) * rho / (p * Gamma);
    }
}