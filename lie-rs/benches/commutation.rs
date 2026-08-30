use commutator_rs::Commutator;
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use lie_rs::LieSeries;
use lyndon_rs::LyndonBasis;
use ordered_float::NotNan;

type S = NotNan<f64>;

fn grid() -> Vec<(usize, usize)> {
    std::env::var("BENCH_GRID")
        .unwrap_or_else(|_| "2x8,2x12,3x8,4x8".to_owned())
        .split(',')
        .map(|tok| {
            let (d, m) = tok.split_once('x').expect("grid entries are `dxm`");
            (d.parse().unwrap(), m.parse().unwrap())
        })
        .collect()
}

fn dense_coefficients(len: usize) -> Vec<S> {
    (0..len)
        .map(|i| NotNan::new((i % 17 + 1) as f64 / 19.0 - 0.4).unwrap())
        .collect()
}

fn make_series(d: usize, m: usize) -> (usize, LieSeries<u8, S>, LieSeries<u8, S>, LyndonBasis<u8>) {
    let basis = LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical);
    let words = basis.generate_basis(m);
    let len = words.len();
    let a = LieSeries::new(words.clone(), dense_coefficients(len));
    let b = LieSeries::new(words, dense_coefficients(len));
    (len, a, b, basis)
}

fn bench_commutation_kernel(c: &mut Criterion) {
    let mut group = c.benchmark_group("commutator_kernel");
    for (d, m) in grid() {
        let (len, a, b, _) = make_series(d, m);
        group.bench_with_input(
            BenchmarkId::new("dxm", format!("{d}x{m}")),
            &(d, m),
            |bencher, _| {
                bencher.iter(|| {
                    let mut out = vec![S::default(); len];
                    LieSeries::commutator_coefficients(
                        &a,
                        &a.coefficients,
                        &b.coefficients,
                        &mut out,
                    );
                    std::hint::black_box(out)
                });
            },
        );
    }
    group.finish();
}

fn bench_commutator(c: &mut Criterion) {
    let mut group = c.benchmark_group("commutator");
    for (d, m) in grid() {
        let (_, a, b, _) = make_series(d, m);
        group.bench_with_input(
            BenchmarkId::new("dxm", format!("{d}x{m}")),
            &(d, m),
            |bencher, _| {
                bencher.iter(|| std::hint::black_box(a.commutator(&b)));
            },
        );
    }
    group.finish();
}

fn bench_lie_series_new(c: &mut Criterion) {
    let mut group = c.benchmark_group("lie_series_new");
    for (d, m) in grid() {
        let (len, _, _, basis) = make_series(d, m);
        group.bench_with_input(
            BenchmarkId::new("dxm", format!("{d}x{m}")),
            &(d, m),
            |bencher, _| {
                bencher.iter(|| LieSeries::new(basis.generate_basis(m), vec![S::default(); len]));
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_commutation_kernel,
    bench_commutator,
    bench_lie_series_new
);
criterion_main!(benches);
