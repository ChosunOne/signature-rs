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

/// Implements the Fredricksen-Kessler-Maiorana algorithm to generate
/// Lyndon words.
fn generate_lyndon_basis(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut basis = Vec::new();
    if k == 0 {
        return basis;
    }
    let mut w = vec![0];

    while !w.is_empty() {
        *w.last_mut().unwrap() += 1;
        let m = w.len();
        if m <= n && w.iter().all(|&x| x >= 1) {
            basis.push(w.clone());
        }

        while w.len() < k {
            w.push(w[w.len() % m]);
        }

        while !w.is_empty() && *w.last().unwrap() == n {
            w.pop();
        }
    }
    basis.sort_by_key(|word| (word.len(), word.clone()));
    basis
}

fn bch_series(max_n: usize) -> Vec<f64> {
    let mut series = vec![1., 1.];
    let mut degree = vec![1, 1];
    let mut prime = vec![1, 1];
    let mut dprime = vec![0, 0];

    let mut i = 3;
    for n in 2..=max_n {
        for j in 1..i {
            for k in j + 1..i {
                if degree[j - 1] + degree[k - 1] != n || j < dprime[k - 1] {
                    continue;
                }
                dprime.push(j);
                prime.push(k);
                degree.push(n);
                i += 1;
            }
        }
    }
    series.truncate(max_n);
    series
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
    assert!(is_lyndon(&word));

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
}
