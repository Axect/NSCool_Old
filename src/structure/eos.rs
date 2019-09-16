use peroxide::{Matrix, Number, ExpLogOps, PowOps, zeros_shape, MutMatrix};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum EOSModel {
    APR,
    FPS,
    SLy4,
}

pub fn load_table(eos: EOSModel) -> Matrix {
    match eos {
        EOSModel::APR => {
            let m = Matrix::read("data/Manufactured/APR.csv", true, ',').expect("Can't read APR");
            let mut n = zeros_shape(m.row, m.col, m.shape);
            for i in 0 .. m.row {
                unsafe {
                    let mut r = n.row_mut(i);
                    let ref_r = &m.row(m.row - 1 - i);
                    for j in 0 .. m.col {
                        *r[j] = ref_r[j];
                    }
                }
            }
            n
        }
        EOSModel::FPS => {
            Matrix::read("data/Manufactured/FPS.csv", true, ',').expect("Can't read FPS")
        }
        EOSModel::SLy4 => {
            Matrix::read("data/Manufactured/SLY4.csv", true, ',').expect("Can't read SLy4")
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Polytrope {
    K: f64,
    Gamma: f64,
}

#[derive(Debug, Clone)]
pub struct PiecewisePolytrope {
    interval: Vec<f64>,
    polytropes: Vec<Polytrope>,
}

pub fn piecewise_poly_fit(data: &Matrix, pieces: usize) -> PiecewisePolytrope {

}

fn piecewise_polytrope(rho: &Vec<f64>, kr: Vec<Number>, ics: Vec<usize>) -> Vec<Number> {
    let k0 = kr[0];
    let gs = kr.into_iter().skip(1).collect::<Vec<Number>>();
    let rhos = ics.into_iter().map(|i| rho[i]).collect::<Vec<f64>>();
    let k1 = (k0 * Number::F(10f64.powf(rho0)).powf(g0)) / Number::F(10f64.powf(rho0)).powf(g1);
    let k2 = (k1 * Number::F(10f64.powf(rho1)).powf(g1)) / Number::F(10f64.powf(rho1)).powf(g2);
    let k3 = (k2 * Number::F(10f64.powf(rho2)).powf(g2)) / Number::F(10f64.powf(rho2)).powf(g3);
    rho.clone()
        .into_iter()
        .map(|x| {
            if x <= rho0 {
                k0.log10() + g0 * x
            } else if x <= rho1 {
                k1.log10() + g1 * x
            } else if x <= rho2 {
                k2.log10() + g2 * x
            } else {
                k3.log10() + g3 * x
            }
        })
        .collect()
}
