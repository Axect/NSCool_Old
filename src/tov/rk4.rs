extern crate peroxide;
use peroxide::*;

use crate::tov::ode::ODE;

pub trait ERK4 {
    fn new(init_param: f64, init_val: Vec<f64>, step: f64, num: usize) -> ODE;
    fn update<F>(&mut self, f: &F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>;
    fn integrate<F>(&mut self, f: F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>;
}

impl ERK4 for ODE {
    fn new(init_param: f64, init_val: Vec<f64>, step: f64, num: usize) -> ODE {
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

    fn update<F>(&mut self, f: &F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>,
    {
        let t = self.param;
        let ys = self.values.clone();
        let h = self.step;
        let h2 = h / 2f64;
        let h3 = h / 3f64;
        let h6 = h / 6f64;

        let k1 = f(t, ys.clone());

        let t_1 = t + h2;
        let t_2 = t + h;

        let k2_add = k1.fmap(|x| x * h);
        let k2 = f(t_1, ys.add(&k2_add));

        let k3_add = k2.fmap(|x| x * h2);
        let k3 = f(t_1, ys.add(&k3_add));

        let k4_add = k3.fmap(|x| x * h2);
        let k4 = f(t_2, ys.add(&k4_add));

        let total_add_part1 = k1.zip_with(|x, y| h6 * x + h3 * y, &k2);
        let total_add_part2 = k3.zip_with(|x, y| h3 * x + h6 * y, &k4);

        let total_add = total_add_part1.zip_with(|x, y| x + y, &total_add_part2);

        self.param = t + h;
        self.values = ys.add(&total_add);

        self.stage += 1;
        self.records[(self.stage, 0)] = self.param;
        for i in 1..self.values.len() + 1 {
            self.records[(self.stage, i)] = self.values[i - 1];
        }
    }

    fn integrate<F>(&mut self, f: F)
    where
        F: Fn(f64, Vec<f64>) -> Vec<f64>,
    {
        for _i in 0..self.num {
            self.update(&f);
        }
    }
}
