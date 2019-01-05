extern crate peroxide;
use peroxide::*;

fn main() {
    let a = matrix(c!(1,2,2,3,1,4,5,2,3),3,3,Row);
    a.inv().unwrap().print();
    a.pseudo_inv().unwrap().print();
}