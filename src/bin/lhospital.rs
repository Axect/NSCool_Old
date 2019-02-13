extern crate peroxide;
use peroxide::*;

fn main() {
    let a = dual(0, 1);
    let b = a.sin();
    let c = a.slope() / b.slope();
    c.print();
}