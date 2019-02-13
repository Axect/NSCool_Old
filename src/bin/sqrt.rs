extern crate peroxide;
use peroxide::*;

fn main() {
    let a = dual(4, 1);
    a.powf(0.5).print();
}