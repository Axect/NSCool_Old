extern crate peroxide;
use peroxide::*;

fn main() {
    let k1_var = Dual::new(0.9, 1.);
    let k2_var = Dual::new(0.9, 1.);

    let k1_const = Dual::new(0.9, 0.);
    let k2_const = Dual::new(0.9, 0.);

    let t11 = t1(k1_var, k2_const);
    let t12 = t1(k1_const, k2_var);
    let t21 = t2(k1_const, k2_var);
    let t22 = t2(k1_var, k2_const);

    t11.print();
    t12.print();
    t21.print();
    t22.print();
}

fn t1(k1: Dual, k2: Dual) -> Dual {
    k1.pow(3) + k2 - Dual::new(1, 0)
}

fn t2(k1: Dual, k2: Dual) -> Dual {
    k2.pow(3) - k1 + Dual::new(1, 0)
}