use ndarray::{Array2, array};
use ordered_float::NotNan;
use signature_rs::log_sig::LogSignatureBuilder;

#[allow(clippy::too_many_lines)]
pub fn main() {
    let builder = LogSignatureBuilder::<u8>::new()
        .with_max_degree(5)
        .with_num_dimensions(4);
    let path: Array2<f64> = array![
        [0.000, 0.000, 0.000, 0.000],
        [1.000, 2.000, 3.000, 4.000],
        [6.000, 5.000, 4.000, 3.000],
        [7.000, 8.000, 9.000, 8.000],
        [12.000, 11.000, 10.000, 9.000],
        [-0.077, 0.042, -0.067, 1.230],
        [-0.154, 0.675, 0.006, -2.340],
        [0.916, 1.177, -0.139, 3.450],
        [1.095, 0.823, -0.261, -5.670],
    ];
    let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));
    let log_sig = builder.build_from_path(&path);
    for (i, c) in log_sig.series.coefficients.iter().enumerate() {
        println!("[{i}]: {c}");
    }
}
