extern crate peroxide;
use peroxide::*;

fn main() {
    let normal = Normal(0f64, 0.1f64);
    let normal2 = Normal(0f64, 100f64);

    let mut x = seq(0., 99., 1f64);
    x = zip_with(|a, b| a + b, &x, &normal.sample(x.len()));

    let mut y = x.fmap(|t| t.powi(2));
    y = zip_with(|a, b| a + b, &y, &normal2.sample(y.len()));

    let n_init = vec![1f64];
    let data = hstack!(x.clone(), y.clone());

    let mut opt = Optimizer::new(data, quad);
    opt.set_init_param(n_init)
        .set_max_iter(30)
        .set_method(LevenbergMarquardt)
        .optimize().print();

    let mut plt = Plot2D::new();
    plt.set_domain(x)
        .insert_image(y)
        .set_legend(vec!["$y=x^2$"])
        .set_title("Test")
        .set_path("figure/lm_test.png")
        .savefig().expect("Can't draw a plot");
}

fn quad(x: &Vec<f64>, n: Vec<Number>) -> Vec<Number> {
    let exp = n[0];
    x.clone().into_iter()
        .map(|t| Number::from_f64(t).powf(exp))
        .collect()
}