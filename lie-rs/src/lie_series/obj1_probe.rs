use super::*;
use lyndon_rs::lyndon::{LyndonBasis, Sort};
use ordered_float::NotNan;

#[test]
#[ignore = "probe"]
fn probe_entry_contribution_ratio() {
    for (d, m) in [(2usize, 12usize), (3usize, 8usize), (4usize, 8usize), (2usize, 8usize), (3usize, 12usize)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let len = words.len();
        let a: LieSeries<u8, NotNan<f64>> =
            LieSeries::new(words, vec![NotNan::new(1.0).unwrap(); len]);
        let table = &a.feasible_decompositions;
        let entries = table.entries();
        let e = entries.len() - 1;
        let c: usize = (0..e)
            .map(|ei| {
                entries[ei + 1].decomp_start as usize - entries[ei].decomp_start as usize
            })
            .sum();
        let rows_gt1 = (0..e)
            .filter(|&ei| {
                entries[ei + 1].decomp_start - entries[ei].decomp_start > 1
            })
            .count();
        println!(
            "d={d} m={m}: basis={len} entries={e} contribs={c} rho={:.3} rows_gt1={rows_gt1} ({:.1}%)",
            c as f64 / e as f64,
            100 * rows_gt1 / e.max(1)
        );
    }
}
