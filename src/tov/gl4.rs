extern crate peroxide;
use crate::tov::ode::ODE;
use peroxide::*;

pub trait GL4 {
    fn gl_new(init_param: f64, init_val: Vec<f64>, step: f64, num: usize) -> ODE;
    fn gl_update<F>(&mut self, f: &F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>;
    fn gl_integrate<F>(&mut self, f: F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>;
}

impl GL4 for ODE {
    fn gl_new(init_param: f64, init_val: Vec<f64>, step: f64, num: usize) -> ODE {
        let l = init_val.len();

        let v1_temp1 = vec![init_param]; // len = 1
        let v1_temp2 = init_val.clone(); // len = l
        let v1 = c!(v1_temp1; v1_temp2);
        let v2 = vec![0f64; (l + 1) * (num)];
        ODE {
            param: init_param,
            values: init_val,
            step: step,
            records: matrix(
                c!(v1; v2),
                num + 1, // Original Value
                l + 1,   // param + ys
                Row,
            ),
            stage: 0,
            num: num,
        }
    }

    fn gl_update<F>(&mut self, f: &F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>,
    {
        let t = self.param;
        let ys = self.values.clone();
        let h = self.step;

        let mnum_11: f64 = 0.5 - 3f64.sqrt() / 6f64;
        let mnum_12: f64 = 0.25 - 3f64.sqrt() / 6f64;
        let mnum_21: f64 = 0.5 + 3f64.sqrt() / 6f64;
        let mnum_22: f64 = 0.25 + 3f64.sqrt() / 6f64;

        let t_1 = t + mnum_11 * h;
        let t_2 = t + mnum_21 * h;

        unimplemented!()
    }

    fn gl_integrate<F>(&mut self, f: F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>,
    {
        unimplemented!()
    }
}
