use signature_rs::BCHCoefficientGenerator;
pub fn main() {
    let mut generator = BCHCoefficientGenerator::<u8>::new(20);
    let series = generator.generate_coefficients();
    println!("{series:#?}");
    println!("{}", series.len());
}
