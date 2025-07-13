use signature_rs::BCHCoefficientGenerator;
use signature_rs::lie_series::LieSeries;
pub fn main() {
    // let mut generator = BCHCoefficientGenerator::<u8>::new(20);
    // let series = generator.generate_coefficients();
    // println!("{series:#?}");
    // println!("{}", series.len());

    let bch_series = LieSeries::<5, char>::new(5);

    for i in 0..bch_series.bch_coefficients.len() {
        println!("{i} \t {}", bch_series.bch_coefficients[i]);
    }
}
