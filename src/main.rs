use num_rational::Ratio;
use signature_rs::BCHCoefficientGenerator;
use signature_rs::bch_series_generator::BchSeriesGenerator;
use signature_rs::lyndon::{LyndonBasis, Sort};
pub fn main() {
    // let mut generator = BCHCoefficientGenerator::<u8>::new(5);
    // let series = generator.generate_coefficients();
    // println!("{series:#?}");

    let basis = LyndonBasis::<char>::new(2, Sort::Lexicographical);
    let bch_series_generator = BchSeriesGenerator::<char>::new(basis, 20);
    let bch_coefficients = bch_series_generator.generate_bch_coefficients::<Ratio<i128>>();

    for (i, c) in bch_coefficients.iter().enumerate() {
        println!("{i} \t {c}",);
    }
}
