extern crate peroxide;
use peroxide::*;

fn main() {
    let a = DataFrame::read_csv("data/Manufactured/SLY4.csv", ' ').unwrap();
    a.print();
}
