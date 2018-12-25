extern crate peroxide;
use peroxide::*;
use NSCool::tov::ode::jacobian;

pub fn main() {
    let x = c!(1,1,1);

    let j = jacobian(x.clone(), fs);
    println!("Jacobian:");
    j.print();
    println!("");

    j.det().print();

    println!("Pseudo inverse of Jacobian:");
    j.pseudo_inv().unwrap().print();
    println!("");

    println!("Function apply x:");
    println!("{:?}", fs(dualize(x.clone())));
    println!("");


    let mut new_x = update(x, &fs);
    new_x.print();

    for _i in 0 .. 9 {
        new_x = update(new_x, &fs);
        new_x.print();
    }
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

fn fs(xs: Vec<Dual>) -> Vec<Dual> {
    let x1 = xs[0];
    let x2 = xs[1];
    let x3 = xs[2];

    vec![f(x1, x2, x3), g(x1, x2, x3), h(x1, x2, x3)]
}

//fn jacobian<F>(x: Vec<f64>, f: F) -> Matrix
//    where F: Fn(Vec<Dual>) -> Vec<Dual>
//{
//    let x_var = x.clone().into_iter()
//        .map(|t| Dual::new(t, 1.))
//        .collect::<Vec<Dual>>();
//    let x_const = x.clone().into_iter()
//        .map(|t| Dual::new(t, 0.))
//        .collect::<Vec<Dual>>();
//
//    let mut j = matrix(vec![0f64; 9], 3, 3, Row);
//
//    let mut vec_temp = x_const.clone();
//
//    for i in 0 .. 3 {
//        vec_temp[i] = x_var[i];
//        let dual_temp = f(vec_temp.clone());
//        let slope_temp = dual_temp.into_iter()
//            .map(|dx| dx.slope())
//            .collect::<Vec<f64>>();
//        for k in 0 .. 3 {
//            j[(k, i)] = slope_temp[k];
//        }
//        vec_temp = x_const.clone();
//    }
//    j
//}

fn extract(x: Dual) -> (f64, f64) {
    (x.value(), x.slope())
}

fn var_const(x: f64) -> (Dual, Dual) {
    (Dual::new(x, 1.), Dual::new(x, 0.))
}

fn update(x: Vec<f64>, f: &Fn(Vec<Dual>) -> Vec<Dual>) -> Vec<f64> {
    let j = jacobian(x.clone(), f);
    let x_dual = dualize(x.clone());
    let fx = f(x_dual).into_iter().map(|t| t.value()).collect::<Vec<f64>>();

    let target = (j.pseudo_inv().unwrap() % fx);

    x.sub(&target.col(0))
}

fn dualize(xs: Vec<f64>) -> Vec<Dual> {
    xs.into_iter().map(|x| Dual::new(x, 0.)).collect::<Vec<Dual>>()
}