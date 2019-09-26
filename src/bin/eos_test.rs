extern crate peroxide;
extern crate NSCool;

use peroxide::Printable;
use NSCool::structure::eos::*;
use NSCool::structure::eos::EOSModel::*;

use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let piece: usize = args[1].parse::<usize>().unwrap();

    let eos_data = load_table(SLy4);
    println!("Load table success");
    let p = piecewise_poly_fit(&eos_data, piece);
    let kgs = p.extract_k_gamma();

    println!("{:?}", kgs);
    p.get_log_interval().print();
}
