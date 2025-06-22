use std::collections::HashMap;

const FACTORIALS: [u64; 21] = [
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

/// Generates the truncated Baker-Campbell-Hausdorff series of a Lyndon basis
/// with the given truncation depth `n` and alphabet size `2`.
pub fn generate_bch_series_2(n: usize) -> Vec<f64> {
    let mut t_n = generate_lyndon_basis(2, n);
    let m_n = t_n.len();
    let mut s = Vec::with_capacity(m_n);
    let bernoulli = bernoulli_sequence(n);
    let mut prime = Vec::with_capacity(m_n);
    let mut dprime = Vec::with_capacity(m_n);
    let mut degree = Vec::with_capacity(m_n);
    let mut word_t_n_map = HashMap::<Vec<usize>, usize>::new();
    for (i, word) in t_n.iter().enumerate() {
        word_t_n_map.insert(word.to_vec(), i);
    }
    let mut i = 0;
    while i < t_n.len() {
        let word = &t_n[i].clone();
        degree.push(word.len());
        let mut s_ui = Vec::new();
        if word.len() == 1 {
            s.push(s_ui);
            i += 1;
            continue;
        }
        let (w, v) = standard_factorization(word);
        if !word_t_n_map.contains_key(v) {
            t_n.push(v.to_vec());
            word_t_n_map.insert(v.to_vec(), t_n.len() - 1);
        }
        if !word_t_n_map.contains_key(w) {
            t_n.push(w.to_vec());
            word_t_n_map.insert(w.to_vec(), t_n.len() - 1);
        }
        prime.push(word_t_n_map[v]);
        dprime.push(word_t_n_map[w]);

        s_ui.push((word_t_n_map[w], word_t_n_map[v]));
        for p in 1..v.len() {
            let v_root = &v[..p];
            let v_comp_w = [(&v[p..]).to_vec(), w.to_vec()].concat();
            if !word_t_n_map.contains_key(v_root) {
                t_n.push(v_root.to_vec());
                word_t_n_map.insert(v_root.to_vec(), t_n.len() - 1);
            }
            if !word_t_n_map.contains_key(&v_comp_w) {
                t_n.push(v_comp_w.clone());
                word_t_n_map.insert(v_comp_w.clone(), t_n.len() - 1);
            }
            s_ui.push((word_t_n_map[v_root], word_t_n_map[&v_comp_w]));
        }
        for q in 1..w.len() {
            let w_root = &w[..q];
            let w_comp_v = [v.to_vec(), (&w[q..]).to_vec()].concat();
            if !word_t_n_map.contains_key(w_root) {
                t_n.push(w_root.to_vec());
                word_t_n_map.insert(w_root.to_vec(), t_n.len() - 1);
            }
            if !word_t_n_map.contains_key(&w_comp_v) {
                t_n.push(w_comp_v.clone());
                word_t_n_map.insert(w_comp_v.clone(), t_n.len() - 1);
            }
            s_ui.push((word_t_n_map[w_root], word_t_n_map[&w_comp_v]));
        }
        s.push(s_ui);
        i += 1;
    }
    let tm_n = t_n.len();
    let mut z_ui = vec![0.; tm_n];
    let mut x_minus_y = vec![0.; tm_n];
    x_minus_y[0] = 1.;
    x_minus_y[1] = -1.;
    let x_minus_y = x_minus_y;
    let mut x_plus_y = vec![0.; tm_n];
    x_plus_y[0] = 1.;
    x_plus_y[1] = 1.;
    let x_plus_y = x_plus_y;
    for i in 0..m_n {
        if i == 0 || i == 1 {
            z_ui[i] = 1.;
            continue;
        }
        // 1/2 [X-Y, Z] (u_i)
        let s_ui = &s[i];
        dbg!(i);
        dbg!(&s_ui);
        let left_term = 0.5 * lie_bracket(&x_minus_y, &z_ui, s_ui);

        // Bernoulli term
        let mut bernoulli_term = 0.0;
        let mut beta = x_plus_y.clone();
        for p in 1..(n / 2) {
            let bernoulli_coef = bernoulli[2 * p] / FACTORIALS[2 * p] as f64;
            for _ in 1..2 * p {
                let mut beta_new = vec![0.; tm_n];
                for x in 0..tm_n {
                    let s_ux = &s[x];
                    beta_new[x] = lie_bracket(&z_ui, &beta, &s_ux);
                }
                beta = beta_new;
            }
            bernoulli_term += bernoulli_coef * lie_bracket(&z_ui, &beta, &s_ui);
        }
        z_ui[i] = (left_term + bernoulli_term) / degree[i] as f64;
    }
    dbg!(&z_ui);

    todo!()
}

/// Computes the Lie bracket [alpha, beta] (u_i)
fn lie_bracket(alpha: &[f64], beta: &[f64], s_ui: &[(usize, usize)]) -> f64 {
    let mut sum = 0.0;
    for &(j, k) in s_ui.iter() {
        sum += alpha[j] * beta[k] - alpha[k] * beta[j];
    }

    sum
}

/// Implements the Fredricksen-Kessler-Maiorana algorithm to generate
/// Lyndon words.
/// `n` is the alphabet size, and `k` is the maximum word length.
fn generate_lyndon_basis(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut basis = Vec::new();
    if k == 0 {
        return basis;
    }
    let mut w = vec![];

    loop {
        if w.is_empty() {
            w = vec![0];
        } else {
            *w.last_mut().unwrap() += 1;
        }

        if !w.is_empty() && w.len() <= k && *w.last().unwrap() < n {
            basis.push(w.clone());

            let m = w.len();
            while w.len() < k {
                w.push(w[w.len() % m]);
            }
        }

        while !w.is_empty() && *w.last().unwrap() >= n {
            w.pop();
        }

        if w.is_empty() {
            break;
        }
    }

    basis.sort_by_key(|word| (word.len(), word.clone()));
    basis
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

pub fn calculate_log_signature<const N: usize>(path: &[[f64; N]], depth: usize) -> Vec<f64> {
    let lyndon_basis = generate_lyndon_basis(N, depth);
    let mut log_sig = vec![0f64; lyndon_basis.len()];
    let mut displacement = vec![0f64; lyndon_basis.len()];

    for (i, window) in path.windows(2).enumerate() {
        displacement.fill(0f64);
        let path_i = window[0];
        let path_i_plus_1 = window[1];

        for j in 0..N {
            displacement[j] = path_i_plus_1[j] - path_i[j];
        }

        if i == 0 {
            log_sig.copy_from_slice(&displacement);
        }
    }

    log_sig
}

fn is_lyndon(word: &[usize]) -> bool {
    let n = word.len();
    if n == 0 {
        return false;
    }

    let mut i = 0;
    let mut j = 1;
    while j < n {
        match word[i].cmp(&word[j]) {
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                i = 0;
                j += 1;
            }
            std::cmp::Ordering::Greater => {
                return false;
            }
        }
    }

    let period = j - i;
    period == n
}

fn standard_factorization(word: &[usize]) -> (&[usize], &[usize]) {
    let n = word.len();
    assert!(n > 1, "Word length must be greater than 1.");

    for split in 1..n {
        if is_lyndon(&word[split..]) {
            return (&word[..split], &word[split..]);
        }
    }
    unreachable!()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_lyndon_basis_generation() {
        let basis = generate_lyndon_basis(5, 5);
        assert_eq!(basis.len(), 829);
        for word in basis {
            assert!(is_lyndon(&word), "not a Lyndon word: {:?}", word);
        }
    }

    #[test]
    fn test_lyndon_standard_factorization() {
        let test_cases = vec![
            // (input, (expected_left, expected_right))
            (vec![0, 1], (vec![0], vec![1])),
            (vec![0, 0, 1], (vec![0], vec![0, 1])),
            (vec![0, 1, 1], (vec![0, 1], vec![1])),
            (vec![0, 1, 1, 1], (vec![0, 1, 1], vec![1])),
            (vec![0, 0, 1, 1], (vec![0], vec![0, 1, 1])),
            (vec![0, 0, 0, 1], (vec![0], vec![0, 0, 1])),
            (vec![0, 0, 1, 0, 1], (vec![0, 0, 1], vec![0, 1])),
            (vec![0, 1, 0, 1, 1], (vec![0, 1], vec![0, 1, 1])),
            (vec![0, 0, 1, 1, 1], (vec![0], vec![0, 1, 1, 1])),
            (vec![0, 0, 0, 1, 1], (vec![0], vec![0, 0, 1, 1])),
            (vec![0, 0, 0, 0, 1], (vec![0], vec![0, 0, 0, 1])),
            // Three-letter alphabet
            (vec![0, 1, 2], (vec![0], vec![1, 2])),
            (vec![0, 0, 1, 2], (vec![0], vec![0, 1, 2])),
            (vec![0, 1, 1, 2], (vec![0], vec![1, 1, 2])),
            (vec![0, 1, 2, 2], (vec![0], vec![1, 2, 2])),
            (vec![0, 1, 2, 3], (vec![0], vec![1, 2, 3])),
        ];

        for (word, (expected_left, expected_right)) in test_cases {
            let (left, right) = standard_factorization(&word);
            assert_eq!(left, &expected_left[..], "Failed for {:?}", word);
            assert_eq!(right, &expected_right[..], "Failed for {:?}", word);
        }
    }

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
    fn test_bch_series() {
        let series = generate_bch_series_2(5);
    }
}
