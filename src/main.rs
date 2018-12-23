extern crate peroxide;
use peroxide::*;

#[allow(non_snake_case)]
use NSCool::tov::rk4::*;

fn main() {
    let mut m = RK4::new(0f64, vec![1f64], 1e-6, 100);
    m.integrate(|t, xs| vec![t * xs[0].sqrt()]);
    let actuals = m.records.col(0).fmap(|x| actual(x));
    let errors = m.records.col(1).zip_with(|x, y| (x - y) / x, &actuals);
    for e in errors {
        println!("{:.e}", e);
    }
}

fn actual(x: f64) -> f64 {
    (1.0 / 16.0) * (x * x + 4.0).powi(2)
}
