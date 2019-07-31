extern crate peroxide;
extern crate natural_unit;
extern crate rayon;
use peroxide::*;
use natural_unit::*;
use rayon::prelude::*;
use std::f64::consts::PI;

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

    let fps_data = Matrix::read("data/fps.dat", false, ' ').expect("Can't read APR");

    let rho = fps_data.col(0).fmap(|x| convert(x, Density, cgs_to_geom.clone()).log10());
    let p = fps_data.col(1).fmap(|x| convert(x, Pressure, cgs_to_geom.clone()).log10());

    rho.len().print();
    p.len().print();

    let kgs = vec![100f64, 2f64, 2f64, 2f64, 2f64, 2f64];

    let data = hstack!(rho.clone(), p.clone());

    let mut opt = Optimizer::new(data, piecewise_polytrope);
    let param = opt.set_init_param(kgs.clone())
        .set_method(LevenbergMarquardt)
        .set_max_iter(100)
        .optimize();
    param.print();

//    let K0 = param[0];
//    let GAMMA0 = param[1];
//    let GAMMA2 = param[2];
//    let GAMMA3 = param[3];
//    let GAMMA4 = param[4];

    let results = (1..101).into_par_iter().map(|idx| {
        // Set Initial Conditions
        let rho_c = convert(0.1057380000E+18 * (idx as f64), Density, cgs_to_geom);
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

        (invert(value[0], Length, cgs_to_geom) / 100f64 / 1000f64, value[1])
    }).collect::<Vec<(f64, f64)>>();

    let results_r: Vec<f64> = results.clone().into_iter().map(|(x,_)| x).collect();
    let results_m: Vec<f64> = results.clone().into_iter().map(|(_,x)| x).collect();

    let mut plt = Plot2D::new();
    plt.set_path("figure/FPS/MR.png")
        .set_domain(results_r)
        .insert_image(results_m)
        .set_legends(vec!["$M-R$"])
        .set_title("$M-R$ relation of FPS")
        .set_xlabel("$R$")
        .set_ylabel("$M$")
        .savefig().expect("Can't draw a figure");
}

fn piecewise_polytrope(rho: &Vec<f64>, kr: Vec<Number>) -> Vec<Number> {
    let k0 = kr[0];
    let g0 = kr[1];
    let g1 = kr[2];
    let g2 = kr[3];
    let g3 = kr[4];
    let g4 = kr[5];
    let rho0 = RHO0;
    let rho1 = RHO1;
    let rho2 = RHO2;
    let rho3 = RHO3;
    let k1 = (k0 * Number::F(10f64.powf(rho0)).powf(g0)) / Number::F(10f64.powf(rho0)).powf(g1);
    let k2 = (k1 * Number::F(10f64.powf(rho1)).powf(g1)) / Number::F(10f64.powf(rho1)).powf(g2);
    let k3 = (k2 * Number::F(10f64.powf(rho2)).powf(g2)) / Number::F(10f64.powf(rho2)).powf(g3);
    let k4 = (k3 * Number::F(10f64.powf(rho3)).powf(g3)) / Number::F(10f64.powf(rho3)).powf(g4);
    rho.clone().into_iter()
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