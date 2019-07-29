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

    let apr_data = Matrix::read("data/APR_EOS_Acc_Fe.dat", false, ' ').expect("Can't read matrix");
    let sly_data = Matrix::read("data/sly4.dat", false, ' ').expect("Can't read matrix");
    let fps_data = Matrix::read("data/fps.dat", false, ' ').expect("Can't read matrix");

    let apr_skipped = apr_data.skip(196, Row);

    let rho = apr_data.col(0).fmap(|x| x.ln());
    let p = apr_data.col(1).fmap(|x| x.ln());

    let rho_skipped = apr_skipped.col(0).fmap(|x| x.ln());
    let p_skipped = apr_skipped.col(1).fmap(|x| x.ln());

    rho_skipped[0].print();
    p_skipped[0].print();

    rho.len().print();
    p.len().print();

    let kg = vec![2f64];
    let data = hstack!(rho_skipped.clone(), p_skipped.clone());
    let mut opt = Optimizer::new(data, polytrope);
    let param = opt.set_init_param(kg)
        .set_method(LevenbergMarquardt)
        .set_max_iter(30)
        .optimize();
    param.print();

    let mut plot = Plot2D::new();
    plot.set_domain(rho)
        .insert_image(p)
        .set_title("$\\log \\rho$ vs $\\log P$")
        .set_xlabel("$\\log\\rho$")
        .set_ylabel("$\\log P$")
        .set_path("figure/test_rho_p.png")
        .set_legends(vec!["$P(\\rho)$"])
        .savefig().expect("Can't draw a plot");
}

fn polytrope(rho: &Vec<f64>, kr: Vec<Number>) -> Vec<Number> {
    let g = kr[0];
    rho.clone().into_iter()
        .map(|x| g * Number::from_f64(x))
        .collect()
}