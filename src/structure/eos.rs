#[allow(unused_imports)]
use peroxide::*;
use rayon::prelude::*;
use natural_unit::{CONSTANT_CGS, ConversionFactor, convert};
use natural_unit::Dimension::{Density, Pressure};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum EOSModel {
    APR,
    FPS,
    SLy4,
}

pub fn load_table(eos: EOSModel) -> Matrix {
    match eos {
        EOSModel::APR => {
            let m = Matrix::read("./data/Manufactured/APR.csv", true, ',').expect("Can't read APR");
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
            Matrix::read("./data/Manufactured/FPS.csv", true, ',').expect("Can't read FPS")
        }
        EOSModel::SLy4 => {
            Matrix::read("./data/Manufactured/SLY4.csv", true, ' ').expect("Can't read SLy4")
        }
    }
}


#[derive(Debug, Clone)]
pub struct PiecewisePolytrope {
    pub log_interval: Vec<f64>,
    pub interval: Vec<f64>,
    pub param: Vec<f64>,
}

impl PiecewisePolytrope {
    pub fn extract_k_gamma(&self) -> Vec<(f64, f64)> {
        let mut ks = vec![0f64; self.log_interval.len() + 1];
        let mut gs = vec![0f64; ks.len()];
        
        ks[0] = self.param[0];
        for (dst, src) in gs.iter_mut().zip(&self.param[1..]) {
            *dst = *src;
        }

        for i in 1 .. ks.len() {
            ks[i] = ks[i-1] + (gs[i] - gs[i-1]) * self.log_interval[i-1];
        }

        ks.into_iter().zip(gs).collect()
    }

    pub fn get_log_interval<'a>(&'a self) -> &'a Vec<f64> {
        &self.log_interval
    }

    pub fn get_interval<'a>(&'a self) -> &'a Vec<f64> {
        &self.interval
    } 
}

#[allow(non_snake_case)]
pub fn piecewise_poly_fit(data: &Matrix, pieces: usize) -> PiecewisePolytrope {
    // Define new unit (c=G=M=1)
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;
    let cgs_to_geom = ConversionFactor::new(1f64 / M, c.powi(2) / (G * M), c.powi(3) / (G * M));

    // Convert unit
    let log_rho = data.col(1).fmap(|x| convert(x, Density, cgs_to_geom.clone()).log10());
    let log_p = data.col(2).fmap(|x| convert(x, Pressure, cgs_to_geom.clone()).log10());

    let kgs = vec![2f64; pieces + 1];
    let new_data = hstack!(log_rho.clone(), log_p.clone());

    // Generate whole candidates
    let ics_vec = iter_candidate(log_rho.len(), pieces - 1);


    let results = ics_vec.into_par_iter()
        .map(|ics| {
            let mut opt = Optimizer::new(new_data.clone(), |rho, kr| piecewise_polytrope(rho, kr, ics.clone()));
            let param = opt
                .set_init_param(kgs.clone())
                .set_method(LevenbergMarquardt)
                .set_max_iter(100)
                .optimize();
            let err = opt.get_error();
            err.print();
            (err, ics, param)
            //if err < rss_error {
            //    rss_error = err;
            //    index = ics.clone();
            //    params = param;
            //}
        })
        .reduce_with(|a, b| {
            if a.0 <= b.0 {
                a
            } else {
                b
            }
        }).unwrap();
    
    let (rss_error, index, params) = results;
    
    let log_interval: Vec<f64> = index.into_iter().map(|i| log_rho[i]).collect();
    let interval = log_interval.fmap(|x| 10f64.powf(x));

    rss_error.print();

    PiecewisePolytrope {
        log_interval,
        interval,
        param: params,
    }
}

fn piecewise_polytrope(rho: &Vec<f64>, kr: Vec<Number>, ics: Vec<usize>) -> Vec<Number> {
    let k0 = kr[0];
    let gs = kr.into_iter().skip(1).collect::<Vec<Number>>();
    let rhos = ics.into_iter().map(|i| rho[i]).collect::<Vec<f64>>();
    let mut ks: Vec<Number> = vec![k0];
    for i in 0 .. gs.len() - 1 {
        let k = ks[i];
        let g = gs[i];
        let g_next = gs[i+1];
        let exp_rho = Number::F(10f64.powf(rhos[i]));

        let next_k = k * exp_rho.powf(g) / exp_rho.powf(g_next);
        ks.push(next_k);
    }

    let mut result: Vec<Number> = vec![Number::F(0f64); rho.len()];
    for i in 0 .. rho.len() {
        let x = rho[i];
        let mut j = 0usize;
        loop {
            let rho_j = rhos[j];
            let g_j = gs[j];
            let k_j = ks[j];
            if x <= rho_j {
                result[i] = k_j.log10() + g_j * x;
                break;
            } else {
                j += 1;
                if j >= rhos.len() {
                    result[i] = ks[j].log10() + gs[j] * x;
                    break;
                }
            }
        }
    }
    result
}

pub fn iter_candidate(n: usize, r: usize) -> Vec<Vec<usize>> {
    let mut result = vec![vec![0usize; r]; C(n, r)];
    let mut curr = (0 .. r).into_iter().collect::<Vec<usize>>();
    result[0].copy_from_slice(&curr[..]);
    for i in 1 .. result.len() {
        let mut j = curr.len() - 1;
        if curr[j] < n-1 {
            curr[j] += 1;
        } else {
            while j > 0 && curr[j] >= n - r + j {
                j -= 1;
            }
            // Escape
            if j == 0 && curr[0] >= n - r {
                break;
            }
            curr[j] += 1;
            for k in (j+1) .. curr.len() {
                curr[k] = curr[k-1] + 1;
            }
        }
        result[i].copy_from_slice(&curr[..]);
    }
    result
}
