extern crate peroxide;
extern crate NSCool;

#[allow(unused_imports)]
use peroxide::*;
use std::error::Error;
use NSCool::structure::eos::*;
use NSCool::structure::eos::EOSModel::*;

use std::env;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let piece: usize = args[1].parse::<usize>().unwrap();

    let eos_data = load_table(SLy4);
    println!("Load table success");
    let p = piecewise_poly_fit(&eos_data, piece);
    let kgs = p.extract_k_gamma();
    let k = kgs.clone().into_iter().map(|x| x.0).collect::<Vec<f64>>();
    let g = kgs.clone().into_iter().map(|x| x.1).collect::<Vec<f64>>();
    let rho = p.get_log_interval();
    let rss = p.get_rss();

    let mut df = DataFrame::with_header(vec!["K", "Gamma", "Rho", "RSS"]);
    df["K"] = k;
    df["Gamma"] = g;
    df["Rho"] = rho.clone();
    df["RSS"] = vec![rss];
    df.print();

    df.write_nc(&format!("data/NC/sly4_piece_{}.nc", piece))?;
    println!("Everything is Okay");
    Ok(())
}
