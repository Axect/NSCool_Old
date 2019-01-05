extern crate peroxide;

use peroxide::*;

#[allow(non_snake_case)]
use NSCool::tov::eos::*;

use std::f64::consts::PI;

pub fn main() {
    let p_c = K*RHO0C.powf(GAMMAF);
    let init_val = c!(0, 0, 0, p_c);
    let mut records = zeros(16_000 + 1, 4);
    records.subs_row(0, init_val.clone());

    let mut next_val = one_step_rk(init_val, tov_rhs, 1e-3);
    let mut idx: usize = 16_000;
    
    for i in 1 .. 16001 {
        records.subs_row(i, next_val.clone());
        next_val = one_step_rk(next_val.clone(), tov_rhs, 1e-3);
        if next_val[3] <= 0f64 {
            idx = i;
            break;
        }
    }

    let new_data = records.data.into_iter().take(idx * 4).collect::<Vec<f64>>();
    matrix(new_data, idx, 4, Row).write("data/tov_rk4.csv").expect("Something Wrong");
}

// TOV
//fn main() {
//    let mut tov_solver = ODE::new(
//        1e-16,
//        c!(0, 0, 1.28*1e-3),
//        1e-3,
//        16 * 1000,
//    );
//
//    tov_solver.integrate(|r, rs| tov_derivative(r, rs));
//
//    tov_solver.param.print();
//    tov_solver.values.print();
//
//    tov_solver.records.write("tov.csv");
//}

//fn tov_derivative(r: f64, rs: Vec<f64>) -> Vec<f64> {
//    let r2 = r.powi(2);
//    let r3 = r.powi(3);
//
//    let m_old = rs[0];
//    let phi_old = rs[1];
//    let rho_old = rs[2];
//
//    let rho2 = rho_old.powi(2);
//    let rho_gamma = rho_old.powi(GAMMA);
//
//    let m_new = 4f64*PI*r2*rho_old;
//    let phi_new = phi_old
//        + (4f64*PI*K*r3 * rho_gamma + m_old)
//        / (r2 - 2f64*r*m_old);
//    let rho_new = -(4f64*PI*K*K*r3 * rho_old.powi(GAMMA+1)
//        + m_old
//        + (4f64*PI*K*r3*rho2 + K*m_old*rho_old))
//        / ((GAMMAF*K*r2 - 2f64*GAMMAF*K*r*m_old));
//
//    vec![m_new, phi_new, rho_new]
//}

fn tov_derivative_dual(rs: Vec<Dual>) -> Vec<Dual> {
    let r = rs[0];
    let r2 = r.pow(2);
    let r3 = r.pow(3);

    let m_old = rs[1];
    let phi_old = rs[2];
    let rho_old = rs[3];

    let rho2 = rho_old.pow(2);
    let rho_gamma = rho_old.pow(GAMMA);

    let m_new = 4f64 * PI * r2 * rho_old;
    let phi_new = phi_old + (4f64 * PI * K * r3 * rho_gamma + m_old) / (r2 - 2f64 * r * m_old);
    let rho_new = -(4f64 * PI * K * K * r3 * rho_old.pow(GAMMA + 1)
        + m_old
        + (4f64 * PI * K * r3 * rho2 + K * m_old * rho_old))
        / (GAMMAF * K * r2 - 2f64 * GAMMAF * K * r * m_old);

    vec![m_new, phi_new, rho_new]
}
