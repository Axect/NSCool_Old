extern crate natural_unit;
extern crate peroxide;
use natural_unit::*;
use peroxide::*;

#[allow(non_snake_case)]
fn main() {
    // Define new unit (c=G=M=1)
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;
    let cgs_to_geom = ConversionFactor::new(1f64 / M, c.powi(2) / (G * M), c.powi(3) / (G * M));

    let sly4_data = Matrix::read("data/sly4.dat", false, ' ').expect("Can't read APR");

    let rho = sly4_data
        .col(0)
        .fmap(|x| convert(x, Density, cgs_to_geom.clone()).log10());
    let p = sly4_data
        .col(1)
        .fmap(|x| convert(x, Pressure, cgs_to_geom.clone()).log10());

    rho.len().print();
    p.len().print();

    let kgs = vec![100f64, 2f64, 2f64, 2f64, 2f64];

    let data = hstack!(rho.clone(), p.clone());

    let mut opt = Optimizer::new(data, piecewise_polytrope);
    let param = opt
        .set_init_param(kgs.clone())
        .set_method(LevenbergMarquardt)
        .set_max_iter(100)
        .optimize();
    param.print();

    let fit = piecewise_polytrope(&rho, NumberVector::from_f64_vec(param)).to_f64_vec();

    let mut plot = Plot2D::new();
    plot.set_domain(rho)
        .insert_image(p)
        .insert_image(fit)
        .set_title("$\\log \\rho$ vs $\\log P$ (Total)")
        .set_xlabel("$\\log\\rho$")
        .set_ylabel("$\\log P$")
        .set_path("figure/SLy/sly4_fit.png")
        .set_legend(vec!["SLy4", "fit"])
        .set_marker(vec![Point, Line])
        .savefig()
        .expect("Can't draw a plot");
}

fn piecewise_polytrope(rho: &Vec<f64>, kr: Vec<Number>) -> Vec<Number> {
    let k0 = kr[0];
    let g0 = kr[1];
    let g1 = kr[2];
    let g2 = kr[3];
    let g3 = kr[4];
    let rho0 = -6.4f64;
    let rho1 = -5.6f64;
    let rho2 = -3.6f64;
    let k1 = (k0 * Number::F(10f64.powf(rho0)).powf(g0)) / Number::F(10f64.powf(rho0)).powf(g1);
    let k2 = (k1 * Number::F(10f64.powf(rho1)).powf(g1)) / Number::F(10f64.powf(rho1)).powf(g2);
    let k3 = (k2 * Number::F(10f64.powf(rho2)).powf(g2)) / Number::F(10f64.powf(rho2)).powf(g3);
    rho.clone()
        .into_iter()
        .map(|x| {
            if x <= rho0 {
                k0.log10() + g0 * x
            } else if x <= rho1 {
                k1.log10() + g1 * x
            } else if x <= rho2 {
                k2.log10() + g2 * x
            } else {
                k3.log10() + g3 * x
            }
        })
        .collect()
}
