extern crate peroxide;
use peroxide::*;

fn main() {
    let mut k1_var = Dual::new(0.9, 1.);
    let mut k2_var = Dual::new(0.9, 1.);

    let mut k1_const = Dual::new(0.9, 0.);
    let mut k2_const = Dual::new(0.9, 0.);

    let mut t11 = t1(k1_var, k2_const);
    let mut t12 = t1(k1_const, k2_var);
    let mut t21 = t2(k1_var, k2_const);
    let mut t22 = t2(k1_const, k2_var);

    let jacobian = matrix(
        vec![
            t11.slope(), t12.slope(), 
            t21.slope(), t22.slope()
            ], 
        2, 2, Row);

    let mut target = matrix(vec![t11.value(), t22.value()], 2, 1, Col);
    
    let mut deltas = - jacobian.inv().unwrap() % target;

    deltas.print();

    k1_var = k1_var + deltas[(0, 0)];
    k2_var = k2_var + deltas[(1, 0)];

    k1_const = Dual::new(k1_var.value(), 0.);
    k2_const = Dual::new(k2_var.value(), 0.);

    k1_const.print();
    k2_const.print();

    for _i in 0 .. 10 {
        t11 = t1(k1_var, k2_const);
        t12 = t1(k1_const, k2_var);
        t21 = t2(k1_var, k2_const);
        t22 = t2(k1_const, k2_var);

        let jacobian = matrix(
        vec![
            t11.slope(), t12.slope(), 
            t21.slope(), t22.slope()
            ], 
        2, 2, Row);

        target = matrix(c!(t11.value(), t21.value()), 2, 1, Col);
        deltas = - jacobian.inv().unwrap() % target;

        k1_var = k1_var + deltas[(0, 0)];
        k2_var = k2_var + deltas[(1, 0)];

        k1_const = Dual::new(k1_var.value(), 0.);
        k2_const = Dual::new(k2_var.value(), 0.);

        k1_const.print();
        k2_const.print();
    }
}

fn t1(k1: Dual, k2: Dual) -> Dual {
    k1.pow(3) + k2 - 1.
}

fn t2(k1: Dual, k2: Dual) -> Dual {
    k2.pow(3) - k1 + 1.
}