use num_bigint::BigInt;
use signature_rs::BCHCoefficientGenerator;
use signature_rs::lie_series::LieSeries;
use signature_rs::lyndon::{Lexicographical, LyndonBasis};
pub fn main() {
    // let mut generator = BCHCoefficientGenerator::<u8>::new(5);
    // let series = generator.generate_coefficients();
    // println!("{series:#?}");

    let basis = LyndonBasis::<2, char, Lexicographical>::generate_basis(20);
    let lie_series = LieSeries::<2, char>::new(basis);
    let bch_coefficients = lie_series.generate_bch_coefficients::<BigInt>();

    for i in 0..bch_coefficients.len() {
        println!("{i} \t {}", bch_coefficients[i]);
    }
}
