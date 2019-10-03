extern crate peroxide;
use peroxide::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let a = DataFrame::read_nc("data/NC/sly4_piece_4.nc", vec!["K","Gamma","Rho","RSS"])?;
    a.print();
    Ok(())
}
