extern crate peroxide;
use peroxide::*;

pub fn main() {
    let x1 = 1f64;
    let x2 = 1f64;
    let x3 = 1f64;
    
    println!("Hello");
}

fn f(x1: Dual, x2: Dual, x3: Dual) -> Dual {
    3f64*x1 - (x2 * x3).cos() - 3f64/2f64
}

fn g(x1: Dual, x2: Dual, x3: Dual) -> Dual {
    4f64*x1.pow(2) -625.*x2.pow(2) + 2f64*x3 - 1f64
}

fn h(x1: Dual, x2: Dual, x3: Dual) -> Dual {
    20.*x3 + (-x1 * x2).exp() + 9f64
}

fn jacobian<F>(x1: f64, x2: f64, x3: f64, f: F, g: F, h: F) -> Matrix 
    where F: Fn(Dual, Dual, Dual) -> Dual
{
    let (x1_var, x1_const) = var_const(x1);
    let (x2_var, x2_const) = var_const(x2);
    let (x3_var, x3_const) = var_const(x3);

    let x_var = vec![x1_var, x2_var, x3_var];
    let x_const = vec![x1_const, x2_const, x3_const];
    let f_list = vec![f, g, h];

    let mut j = matrix(vec![0f64; 9], 3, 3, Row);

    for i in 0 .. 3 {

        for j in 0 .. 3 {
            
        }
    }
    unimplemented!()
}

fn extract(x: Dual) -> (f64, f64) {
    (x.value(), x.slope())
}

fn var_const(x: f64) -> (Dual, Dual) {
    (Dual::new(x, 1.), Dual::new(x, 0.))
}
