extern crate peroxide;
use peroxide::*;
extern crate natural_unit;
use natural_unit::CONSTANT_CGS;

pub const k_b_cgs: f64 = CONSTANT_CGS.k_b;
pub const c_cgs: f64 = CONSTANT_CGS.c;
pub const m_e_cgs: f64 = 9.1093837015e-28;

#[derive(Debug, Clone)]
pub struct NSMech {
    pub t: f64,
    pub rho: f64,
    pub a: f64,
    pub z: f64,
}

pub trait Emissivity {
    fn nplasma(&self) -> f64;
}

impl Emissivity for NSMech {
    fn nplasma(&self) -> f64 {
        assert!(self.t > 1e+7 && self.t < 1e+11, "Non-valid plasma temperature");

        if self.z == 0f64 {
            return 0f64;
        }

        let t = self.t;
        let mu_e = self.a / self.z;
        let rho = self.rho;

        let lambda = t * k_b_cgs / (m_e_cgs * c_cgs.powi(2));
        let gamma2 = 1.10985e+11 * rho / mu_e;
        let gamma = gamma2.sqrt();
        let gamma_half = gamma.sqrt();

        let f_t = 2.4 + 0.6 * gamma_half + 0.51 * gamma + 1.25 * gamma_half.powi(3);
        let f_l = (8.6 * gamma2 + 1.35 * gamma_half.powi(7)) / (225f64 - 17f64 * gamma + gamma2);

        let x = 1f64 / 6f64 * (17.5 + (2f64 * rho / mu_e).log10() - 3f64 * t.log10());
        let y = 1f64 / 6f64 * (-24.5 + (2f64 * rho / mu_e).log10() + 3f64 * t.log10());

        let f_xy = if x.abs() > 0.7 || y < 0f64 {
            1f64
        } else {
            1.05 + (0.39
                - 1.25 * x
                - 0.35 * (4.5 * x).sin()
                - 0.3 * (-(4.5 * x + 0.9).powi(2)).exp())
                * (-(0f64.min(y - 1.6 + 1.25 * x) / (0.57 - 0.25 * x)).powi(2)).exp()
        };

        3e+21 * lambda.powi(9) * gamma2.powi(3) * (-gamma).exp() * (f_t + f_l) * f_xy
    }
}
