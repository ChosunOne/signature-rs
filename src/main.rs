use num_bigint::BigInt;
use signature_rs::BCHCoefficientGenerator;
use signature_rs::lie_series::LieSeries;
pub fn main() {
    // let mut generator = BCHCoefficientGenerator::<u8>::new(5);
    // let series = generator.generate_coefficients();
    // println!("{series:#?}");

    let bch_series = LieSeries::<2, char, BigInt>::new(5);

    for i in 0..bch_series.bch_coefficients.len() {
        println!("{i} \t {}", bch_series.bch_coefficients[i]);
    }
}
