extern crate natural_unit;
extern crate peroxide;
extern crate NSCool;
use natural_unit::*;
use peroxide::*;
use NSCool::structure::eos::{load_table, EOSModel};

const RHO: [f64; 3] = [-11.7720, -6.8641, -4.2205];
const KG: [(f64, f64); 4] = [(373981.6209505603, 1.8485919406272393), (0.34056576252855564, 1.3354533321504343), (0.0017649970672523633, 1.002496127574712), (340.0320488034904, 2.2546726199881117)];

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
        .set_path("figure/SLy/sly4_fit_bruteforce.png")
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
            } else {
                KG[3].0.log10() + KG[3].1 * x
            }
        })
        .collect()
}
