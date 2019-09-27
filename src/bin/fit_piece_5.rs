extern crate peroxide;
extern crate natural_unit;
extern crate NSCool;
use peroxide::*;
use natural_unit::*;
use NSCool::structure::eos::{load_table, EOSModel};

const RHO: [f64; 4] = [-11.7720, -6.9881, -3.9332, -2.5193];
const KG: [(f64, f64); 5] = [(355327.00867647934, 1.8470827227990216), (0.3640468463692981, 1.3382915899844015), (0.002789299645770813, 1.035539299065058), (4975.6233524624695, 2.6249154018728196), (17.056449830011807, 1.646480570265483)];

#[allow(non_snake_case)]
fn main() {
    // Define new unit (c=G=M=1)
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;
    let cgs_to_geom = ConversionFactor::new(1f64 / M, c.powi(2) / (G * M), c.powi(3) / (G * M));

    let sly4_data = load_table(EOSModel::SLy4);

    let rho = sly4_data
        .col(1)
        .fmap(|x| convert(x, Density, cgs_to_geom.clone()).log10());
    let p = sly4_data
        .col(2)
        .fmap(|x| convert(x, Pressure, cgs_to_geom.clone()).log10());

    rho.len().print();
    p.len().print();

    let data = hstack!(rho.clone(), p.clone());

    let fit = piecewise_polytrope(&rho);

    let mut plot = Plot2D::new();
    plot.set_domain(rho)
        .insert_image(p)
        .insert_image(fit)
        .set_title("$\\log \\rho$ vs $\\log P$ (Total)")
        .set_xlabel("$\\log\\rho$")
        .set_ylabel("$\\log P$")
        .set_path("figure/SLy/sly4_piece5_bruteforce.png")
        .set_legend(vec!["SLy4", "fit"])
        .set_marker(vec![Point, Line])
        .savefig()
        .expect("Can't draw a plot");
}

fn piecewise_polytrope(log_rho: &Vec<f64>) -> Vec<f64> {
    log_rho.clone()
        .into_iter()
        .map(|x| {
            if x <= RHO[0] {
                KG[0].0.log10() + KG[0].1 * x
            } else if x <= RHO[1] {
                KG[1].0.log10() + KG[1].1 * x
            } else if x <= RHO[2] {
                KG[2].0.log10() + KG[2].1 * x
            } else if x <= RHO[3] {
                KG[3].0.log10() + KG[3].1 * x
            } else {
                KG[4].0.log10() + KG[4].1 * x
            }
        })
        .collect()
}
