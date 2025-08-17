use ndarray::{Array2, array};
use ordered_float::NotNan;
use signature_rs::log_sig::LogSignatureBuilder;

#[allow(clippy::too_many_lines)]
pub fn main() {
    let builder = LogSignatureBuilder::<u8>::new()
        .with_max_degree(3)
        .with_num_dimensions(20);
    let path: Array2<f64> = array![
        [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
            20.
        ],
        [
            1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.10, 11.11, 12.12, 13.13, 14.14, 15.15,
            16.16, 17.17, 18.18, 19.19, 20.20
        ],
        [
            1.1, 2.2, 3.3, 4.4, 5.5, 6.5, 7.7, 8.8, 9.9, 10., 11., 12., 13., 14., 15., 16., 17.,
            18., 19., 20.
        ],
        [
            20., 19., 18., 17., 16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2.,
            1.
        ],
        [
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16,
            0.17, 0.18, 0.19, 0.20
        ],
    ];
    let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));
    let log_sig = builder.build_from_path(&path);

    for (i, c) in log_sig.series.coefficients.iter().enumerate() {
        println!("[{i}]: {c}");
    }
}
