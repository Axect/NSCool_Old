extern crate peroxide;
use peroxide::*;

use NSCool::tov::ode::jacobian;

pub fn main() {
    let mut x = c!(1,2,3);
    let j = jacobian(x.clone(), f);
    j.print();

    for i in 0 .. 10 {
        x = update(x, f);
        x.print();
    }
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let x1 = xs[0];
    let x2 = xs[1];
    let x3 = xs[2];

    vec![
        x1.pow(2) - 2.*x1 + x2.pow(2) - x3 + 1.,
        x1*x2.pow(2) - x1 - 3.*x2 + x2*x3 + 2.,
        x1*x3.pow(2) - 3.*x3 + x2*x3.pow(2) + x1*x2,
    ]
}

fn update<F>(xs: Vec<f64>, f: F) -> Vec<f64>
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    let j = jacobian(xs.clone(), f);
    let pinv_j = j.pseudo_inv().unwrap();
    let xs_dual = xs.dualize(Kind::Const);
    let fx = f(xs_dual).extract();

    xs.sub(&(pinv_j % fx).col(0))
}

enum Kind {
    Var,
    Const
}

trait Vec2Dual {
    fn dualize(&self, kind: Kind) -> Vec<Dual>;
}

impl Vec2Dual for Vec<f64> {
    fn dualize(&self, kind: Kind) -> Vec<Dual> {
        match kind {
            Kind::Const => {
               self.clone().into_iter()
                   .map(|t| Dual::new(t, 0.))
                   .collect::<Vec<Dual>>()
            },
            Kind::Var => {
                self.clone().into_iter()
                    .map(|t| Dual::new(t, 1.))
                    .collect::<Vec<Dual>>()
            }
        }
    }
}

trait Dual2Vec {
    fn extract(&self) -> Vec<f64>;
}

impl Dual2Vec for Vec<Dual> {
    fn extract(&self) -> Vec<f64> {
        self.clone().into_iter()
            .map(|t| t.value())
            .collect::<Vec<f64>>()
    }
}