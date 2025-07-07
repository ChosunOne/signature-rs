use std::collections::HashMap;

use crate::{lyndon::LyndonBasis, rooted_tree::RootedTree};
pub mod bch;
pub mod lyndon;
pub mod rooted_tree;

pub const FACTORIALS: [u64; 21] = [
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600,
    6227020800,
    87178291200,
    1307674368000,
    20922789888000,
    355687428096000,
    6402373705728000,
    121645100408832000,
    2432902008176640000,
];
pub const PRIMES: [usize; 25] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 6767, 71, 73, 79, 83, 89,
    97,
];

/// Generates the truncated Baker-Campbell-Hausdorff series of a Lyndon basis
/// with the given truncation depth `n` and alphabet size `2`.
pub fn generate_bch_series_2(n: usize) -> Vec<f64> {
    let mut t_n = LyndonBasis::<2, u8>::generate_basis(n)
        .into_iter()
        .map(|w| RootedTree::from(w))
        .collect::<Vec<_>>();
    let m_n = t_n.len();
    let bernoulli = bernoulli_sequence(n);
    let mut degree = Vec::with_capacity(m_n);
    let mut tree_t_n_map = HashMap::<RootedTree<u8>, usize>::new();
    for (i, tree) in t_n.iter().enumerate() {
        degree.push(tree.degree());
        tree_t_n_map.insert(tree.clone(), i);
    }
    let s = build_s(&mut t_n, &mut tree_t_n_map, &mut degree);
    let tm_n = t_n.len();
    let mut z_ui = vec![0.; tm_n];
    let mut is_computed_z_ui = vec![false; tm_n];
    let mut x_minus_y = vec![0.; tm_n];
    x_minus_y[0] = 1.;
    x_minus_y[1] = -1.;
    let x_minus_y = x_minus_y;
    let mut x_plus_y = vec![0.; tm_n];
    x_plus_y[0] = 1.;
    x_plus_y[1] = 1.;
    let x_plus_y = x_plus_y;
    let mut prime = vec![0; m_n];
    let mut d_prime = vec![0; m_n];
    let mut sigma = vec![0_usize; m_n];
    sigma[0] = 1;
    sigma[1] = 1;
    let mut kappa = vec![0; m_n];
    kappa[0] = 1;
    kappa[1] = 1;
    for i in 0..m_n {
        z_ui[i] = compute_z_ui(
            i,
            &s,
            &x_minus_y,
            &x_plus_y,
            &bernoulli,
            &degree,
            &mut is_computed_z_ui,
            &mut z_ui,
        );

        if i > 1 {
            prime[i] = s[i][0].0;
            d_prime[i] = s[i][0].1;
            kappa[i] = {
                if d_prime[prime[i]] != d_prime[i] {
                    1
                } else {
                    kappa[prime[i]] + 1
                }
            };
            sigma[i] = kappa[i] * sigma[prime[i]] * sigma[d_prime[i]];
        }
    }
    z_ui[..m_n]
        .iter()
        .zip(sigma.iter())
        .map(|(&z, &sig)| z / sig as f64)
        .collect()
}

fn build_s(
    t_n: &mut Vec<RootedTree<u8>>,
    tree_t_n_map: &mut HashMap<RootedTree<u8>, usize>,
    degree: &mut Vec<usize>,
) -> Vec<Vec<(usize, usize)>> {
    let mut s = vec![vec![]; t_n.len()];
    let mut i = 0;
    while i < t_n.len() {
        let tree = &t_n[i];
        let Some((v, w)) = tree.factorize() else {
            i += 1;
            continue;
        };
        let v_idx = tree_t_n_map[&v];
        let w_idx = tree_t_n_map[&w];
        s[i].push((v_idx, w_idx));

        for p in 0..v.degree() - 1 {
            let s_v = &s[v_idx];
            let mut v_root_w = t_n[s_v[p].0].clone();
            v_root_w.graft(w.clone());
            let v_comp = t_n[s_v[p].1].clone();
            if !tree_t_n_map.contains_key(&v_root_w) {
                t_n.push(v_root_w.clone());
                degree.push(v_root_w.degree());
                tree_t_n_map.insert(v_root_w.clone(), t_n.len() - 1);
                s.push(vec![]);
            }
            s[i].push((tree_t_n_map[&v_root_w], tree_t_n_map[&v_comp]));
        }
        for q in 0..w.degree() - 1 {
            let s_w = &s[w_idx];
            let mut v_w_root = v.clone();
            let w_root = t_n[s_w[q].0].clone();
            v_w_root.graft(w_root);
            let w_comp = t_n[s_w[q].1].clone();
            if !tree_t_n_map.contains_key(&v_w_root) {
                t_n.push(v_w_root.clone());
                degree.push(v_w_root.degree());
                tree_t_n_map.insert(v_w_root.clone(), t_n.len() - 1);
                s.push(vec![]);
            }
            s[i].push((tree_t_n_map[&v_w_root], tree_t_n_map[&w_comp]));
        }
        i += 1;
    }

    s
}

fn compute_z_ui(
    i: usize,
    s: &[Vec<(usize, usize)>],
    x_minus_y: &[f64],
    x_plus_y: &[f64],
    bernoulli: &[f64],
    degree: &[usize],
    is_computed_z_ui: &mut [bool],
    z_ui: &mut [f64],
) -> f64 {
    if i == 0 || i == 1 {
        is_computed_z_ui[i] = true;
        return 1.;
    }
    // 1/2 [X-Y, Z] (u_i)
    let s_ui = &s[i];
    for &(j, k) in s_ui {
        // Check if there are any Z_u terms that we need to complete first
        if !is_computed_z_ui[j] {
            z_ui[j] = compute_z_ui(
                j,
                s,
                x_minus_y,
                x_plus_y,
                bernoulli,
                degree,
                is_computed_z_ui,
                z_ui,
            );
        }
        if !is_computed_z_ui[k] {
            z_ui[k] = compute_z_ui(
                k,
                s,
                x_minus_y,
                x_plus_y,
                bernoulli,
                degree,
                is_computed_z_ui,
                z_ui,
            );
        }
    }
    let left_term = 0.5 * lie_bracket(&x_minus_y, &z_ui, s_ui);

    // Bernoulli term
    let mut bernoulli_term = 0.0;
    for p in 1..=((degree[i] - 1) / 2) {
        let bernoulli_coef = bernoulli[2 * p] / FACTORIALS[2 * p] as f64;
        let adjoint_term = adjoint_operator(&z_ui, &x_plus_y, s, i, 2 * p);
        bernoulli_term += bernoulli_coef * adjoint_term;
    }
    is_computed_z_ui[i] = true;
    let result = (left_term + bernoulli_term) / degree[i] as f64;
    result
}

/// Computes the Lie bracket [alpha, beta] (u_i)
fn lie_bracket(alpha: &[f64], beta: &[f64], s_ui: &[(usize, usize)]) -> f64 {
    let mut sum = 0.0;
    for &(j, k) in s_ui.iter() {
        sum += alpha[j] * beta[k] - alpha[k] * beta[j];
    }

    sum
}

fn adjoint_operator(
    alpha: &[f64],
    beta: &[f64],
    s: &[Vec<(usize, usize)>],
    i: usize,
    power: usize,
) -> f64 {
    if power == 0 {
        return beta[i];
    }
    if power == 1 {
        return lie_bracket(alpha, beta, &s[i]);
    }
    let mut sum = 0.0;
    for &(j, k) in s[i].iter() {
        sum += alpha[j] * adjoint_operator(alpha, beta, s, k, power - 1)
            - alpha[k] * adjoint_operator(alpha, beta, s, j, power - 1);
    }

    sum
}

fn binomial(n: usize, k: usize) -> f64 {
    if k == 0 {
        return 1.0;
    }

    let k = k.min(n - k);
    let mut result = 1.0;

    for i in 0..k {
        result *= (n - i) as f64 / (i + 1) as f64;
    }

    result
}

pub fn bernoulli_sequence(max_n: usize) -> Vec<f64> {
    let mut b = vec![0f64; max_n + 1];

    b[0] = 1.0;
    if max_n > 0 {
        b[1] = -0.5;
    }

    for n in 2..=max_n {
        if n % 2 == 1 {
            continue;
        }
        let mut sum = 0.0;
        for k in 0..n {
            sum += binomial(n + 1, k) * b[k];
        }

        b[n] = -sum / (n + 1) as f64;
    }

    // Use positive sign convention
    if max_n > 0 {
        b[1] = 0.5;
    }

    b
}

#[cfg(test)]
mod test {
    use crate::lyndon::LyndonWord;

    use super::*;

    #[test]
    fn test_lie_bracket() {
        let x = vec![1., 2., 3.];
        let y = vec![4., 5., 6.];
        let s_ui = vec![(0, 1), (1, 2)];
        let z = lie_bracket(&x, &y, &s_ui);
        assert_eq!(z, x[0] * y[1] - x[1] * y[0] + x[1] * y[2] - x[2] * y[1]);
    }

    // #[test]
    // fn test_adjoint_operator() {
    //     let mut t_n = LyndonWord::<u8>::generate_basis(2, 3)
    //         .into_iter()
    //         .map(|w| RootedTree::from(w))
    //         .collect::<Vec<_>>();
    //     let mut degree = vec![];
    //     let mut tree_t_n_map = HashMap::<RootedTree<u8>, usize>::new();
    //     for (i, tree) in t_n.iter().enumerate() {
    //         degree.push(tree.degree());
    //         tree_t_n_map.insert(tree.clone(), i);
    //     }
    //     let z_ui = vec![1., 1., 0.5, 0., 0., 0.];
    //     let x_plus_y = vec![1., 1., 0., 0., 0.];
    //     let s = build_s(&mut t_n, &mut tree_t_n_map, &mut degree);
    //     let z = adjoint_operator(&z_ui, &x_plus_y, &s, 3, 3);
    //     assert_eq!(z, 0.0);
    //
    //     let z = adjoint_operator(&z_ui, &x_plus_y, &s, 0, 1);
    //     let expected_z = lie_bracket(&z_ui, &x_plus_y, &s[0]);
    //     assert_eq!(z, expected_z);
    // }

    #[test]
    fn test_bernoulli_sequence() {
        let seq = bernoulli_sequence(10);
        assert_eq!(
            seq,
            [
                1.0,
                0.5,
                0.16666666666666666,
                0.0,
                -0.033333333333333305,
                0.0,
                0.023809523809523662,
                0.0,
                -0.03333333333333233,
                0.0,
                0.0757575757575662,
            ]
        );
    }

    #[test]
    fn test_build_s() {
        let mut t_n = vec![
            vec![0],
            vec![1],
            vec![0, 1],
            vec![0, 1, 1],
            vec![0, 0, 1],
            vec![0, 1, 1, 1],
            vec![0, 0, 1, 1],
            vec![0, 0, 0, 1],
            vec![0, 1, 1, 1, 1],
            vec![0, 0, 1, 0, 1],
            vec![0, 1, 0, 1, 1],
            vec![0, 0, 1, 1, 1],
            vec![0, 0, 0, 1, 1],
            vec![0, 0, 0, 0, 1],
        ]
        .into_iter()
        .map(|v| LyndonWord::<2, u8>::try_from(v).expect("To make a lyndon word"))
        .map(|w| RootedTree::from(w))
        .collect::<Vec<_>>();
        let m_n = t_n.len();
        let mut degree = vec![];
        let mut tree_t_n_map = HashMap::<RootedTree<u8>, usize>::new();
        for (i, tree) in t_n.iter().enumerate() {
            degree.push(tree.degree());
            tree_t_n_map.insert(tree.clone(), i);
        }

        let s = build_s(&mut t_n, &mut tree_t_n_map, &mut degree);
        let expected_s = vec![
            vec![],                                             // X
            vec![],                                             // Y
            vec![(0, 1)],                                       // XY
            vec![(2, 1), (2, 1)],                               // XYY
            vec![(0, 2), (m_n, 1)],                             // XXY
            vec![(3, 1), (3, 1), (3, 1)],                       // XYYY
            vec![(0, 3), (4, 1), (4, 1)],                       // XXYY
            vec![(0, 4), (m_n, 2), (m_n + 1, 1)],               // XXXY
            vec![(5, 1), (5, 1), (5, 1), (5, 1)],               // XYYYY
            vec![(4, 2), (4, 2), (m_n + 2, 1), (m_n + 2, 1)],   // XXYXY
            vec![(2, 3), (6, 1), (m_n + 3, 1), (m_n + 3, 1)],   // XYXYY
            vec![(0, 5), (6, 1), (6, 1), (6, 1)],               // XXYYY
            vec![(0, 6), (m_n, 3), (7, 1), (7, 1)],             // XXXYY
            vec![(0, 7), (m_n, 4), (m_n + 1, 2), (m_n + 4, 1)], // XXXXY
            // Auxiliary S values
            vec![(0, 0)],                                 // (X) -> (X)
            vec![(0, m_n), (m_n, 0)],                     // (X) -> (X) -> (X)
            vec![(m_n, 2), (4, 0), (m_n + 5, 1)],         // (X) <- (X) -> (X) -> (Y)
            vec![(2, 2), (4, 1), (m_n + 6, 1)],           // (Y) <- (X) -> (X) -> (Y)
            vec![(0, m_n + 1), (m_n, m_n), (m_n + 1, 0)], // (X) -> (X) -> (X) -> (X)
            vec![(m_n, 0), (m_n, 0)],                     // (X) <- (X) -> (X)
            vec![(m_n, 1), (2, 0)],                       // (Y) <- (X) -> (X)
        ];
        for (i, (s_ui, expected_s_ui)) in s.iter().zip(expected_s.iter()).enumerate() {
            dbg!(i);
            assert_eq!(s_ui, expected_s_ui);
        }
    }

    #[test]
    fn test_bch_series() {
        let series = generate_bch_series_2(5);
        let expected_z_ui_series = vec![
            1.,
            1.,
            0.5,
            1. / 12.,
            1. / 12.,
            0.,
            1. / 24.,
            0.,
            -1. / 720.,
            1. / 360.,
            1. / 120.,
            1. / 180.,
            1. / 180.,
            -1. / 720.,
        ];
        dbg!(&series);
        dbg!(&expected_z_ui_series);
        // assert_eq!(series, expected_z_ui_series);
    }
}
