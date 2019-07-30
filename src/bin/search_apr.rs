extern crate peroxide;
extern crate natural_unit;
use peroxide::*;
use natural_unit::*;

#[allow(non_snake_case)]
fn main() {
    // Define new unit (c=G=M=1)
    let c = CONSTANT_CGS.c;
    let G = CONSTANT_CGS.G;
    let M = CONSTANT_CGS.m_solar;
    let cgs_to_geom = ConversionFactor::new(
        1f64 / M,
        c.powi(2) / (G*M),
        c.powi(3) / (G*M),
    );

    let apr_data = Matrix::read("data/APR_EOS_Acc_Fe.dat", false, ' ').expect("Can't read APR");
    let apr_crust = Matrix::read("data/APR_EOS_Acc_Fe_crust.dat", false, ' ').expect("Can't read crust");

    let rho = apr_data.col(0).fmap(|x| convert(x, Density, cgs_to_geom.clone()).log10());
    let p = apr_data.col(1).fmap(|x| convert(x, Pressure, cgs_to_geom.clone()).log10());
    let rho_crust = apr_crust.col(0).fmap(|x| convert(x, Density, cgs_to_geom).log10());
    let p_crust = apr_crust.col(1).fmap(|x| convert(x, Pressure, cgs_to_geom).log10());

    rho.len().print();
    p.len().print();

    let kg = vec![100f64, 2f64];

    let data = hstack!(rho.clone(), p.clone());
    let data_crust = hstack!(rho_crust.clone(), p_crust.clone());

    let mut opt = Optimizer::new(data, polytrope);
    let param = opt.set_init_param(kg.clone())
        .set_method(LevenbergMarquardt)
        .set_max_iter(100)
        .optimize();
    param.print();

    let mut opt_crust = Optimizer::new(data_crust, polytrope);
    let param_crust = opt_crust.set_init_param(kg.clone())
        .set_method(LevenbergMarquardt)
        .set_max_iter(100)
        .optimize();
    param_crust.print();

    let fit = polytrope(&rho, NumberVector::from_f64_vec(param)).to_f64_vec();
    let fit_crust = polytrope(&rho_crust, NumberVector::from_f64_vec(param_crust)).to_f64_vec();

    let mut plot = Plot2D::new();
    plot.set_domain(rho)
        .insert_image(p)
        .insert_image(fit)
        .set_title("$\\log \\rho$ vs $\\log P$ (Inner)")
        .set_xlabel("$\\log\\rho$")
        .set_ylabel("$\\log P$")
        .set_path("figure/APR/apr_fit.png")
        .set_legends(vec!["APR", "fit"])
        .set_marker(vec![Point, Line])
        .savefig().expect("Can't draw a plot");

    let mut plot = Plot2D::new();
    plot.set_domain(rho_crust)
        .insert_image(p_crust)
        .insert_image(fit_crust)
        .set_title("$\\log \\rho$ vs $\\log P$ (Crust)")
        .set_xlabel("$\\log\\rho$")
        .set_ylabel("$\\log P$")
        .set_path("figure/APR/apr_crust_fit.png")
        .set_legends(vec!["APR", "fit"])
        .set_marker(vec![Point, Line])
        .savefig().expect("Can't draw a plot");
}

fn polytrope(rho: &Vec<f64>, kr: Vec<Number>) -> Vec<Number> {
    let k = kr[0];
    let g = kr[1];
    rho.clone().into_iter()
        .map(|x| k.log10() + g * x)
        .collect()
}