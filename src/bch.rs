use num_rational::Ratio;

use crate::{
    Int,
    constants::{FACTORIALS, PRIMES},
};

pub fn p_adic_expansion(n: usize, p: i128) -> Vec<i128> {
    let mut num = n as i128;
    let mut max_power = 0;
    while p.saturating_pow(max_power) <= num {
        max_power += 1;
    }
    let mut alphas = vec![0i128; max_power as usize];
    for (i, power) in (0..max_power).enumerate().rev() {
        alphas[i] = num / p.pow(power);
        num -= alphas[i] * p.pow(power);
    }

    alphas
}
pub fn s_p(n: usize, p: i128) -> i128 {
    p_adic_expansion(n, p).iter().sum()
}

pub fn bch_denominator<U: Int>(n: usize) -> U {
    let primes = PRIMES.iter().filter(|&&x| x < n as i128);
    let mut prod = U::from(1);
    for &p in primes {
        let s_p_n = s_p(n, p);
        let mut t = 0;
        while p.pow(t) <= s_p_n {
            t += 1;
        }
        prod *= U::from(p.pow(t - 1));
    }
    prod
}

/// Calculates the numerator of the Goldberg coefficient.
pub fn goldberg_coeff_numerator<U: Int>(q: Vec<usize>, a_first: bool) -> U {
    let n = q.iter().sum();
    let d = U::from(FACTORIALS[n]) * bch_denominator(n);
    let m = q.len();
    let mut c = vec![U::from(0); n * n];
    let mut a_current = a_first;
    if m % 2 == 0 {
        a_current = !a_first;
    }
    let mut l = 0;
    for i in (0..m).rev() {
        for r in 1..=q[i] {
            l += 1;
            let mut h = U::from(0);
            if i == m - 1 {
                h = d.clone() / FACTORIALS[l].into();
            } else if a_current && i == m - 2 {
                h = d.clone() / (FACTORIALS[r] * FACTORIALS[q[i + 1]]).into();
            }
            c[(l - 1) * n] = h;
            for k in 2..l {
                h = U::from(0);
                for j in 1..=r {
                    if l > j && c[k - 2 + n * (l - j - 1)] != U::from(0) {
                        h += c[k - 2 + n * (l - j - 1)].clone() / U::from(FACTORIALS[j]);
                    }
                }
                if a_current && i <= m - 2 {
                    for j in 1..=q[i + 1] {
                        if l > r + j && c[k - 2 + n * (l - r - j - 1)] != U::from(0) {
                            h += c[k - 2 + n * (l - r - j - 1)].clone()
                                / U::from(FACTORIALS[r] * FACTORIALS[j]);
                        }
                    }
                }
                c[k - 1 + n * (l - 1)] = h;
            }
            c[l - 1 + n * (l - 1)] = d.clone();
        }
        a_current = !a_current;
    }
    let mut sum = U::from(0);
    for k in 1..=n {
        let denom: U = if k % 2 != 0 {
            (k as u64).into()
        } else {
            (-(k as i64)).into()
        };
        sum += c[k - 1 + n * (n - 1)].clone() / denom;
    }
    sum
}

#[cfg(test)]
mod test {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(100, 11, vec![1, 9])]
    #[case(11, 2, vec![1, 1, 0, 1])]
    #[case(123, 2, vec![1, 1, 0, 1, 1, 1, 1])]
    fn test_p_adic_expansion(
        #[case] n: usize,
        #[case] p: i128,
        #[case] expected_expansion: Vec<i128>,
    ) {
        let expansion = p_adic_expansion(n, p);
        assert_eq!(expansion, expected_expansion);
    }

    #[test]
    fn test_bch_denominators() {
        let n = 25;
        let denominators = (1..n + 1)
            .map(|x| bch_denominator::<i128>(x))
            .collect::<Vec<_>>();
        let expected_denominators = [
            1, 1, 2, 1, 6, 2, 6, 3, 10, 2, 6, 2, 210, 30, 12, 3, 30, 10, 210, 42, 330, 30, 60, 30,
            546,
        ];
        assert_eq!(denominators, expected_denominators);
    }

    #[rstest]
    #[case(vec![1], true, Ratio::new(1.into(), 1.into()))]
    #[case(vec![1], false, Ratio::new(1.into(), 1.into()))]
    #[case(vec![1, 1], true, Ratio::new(1.into(), 2.into()))]
    #[case(vec![2, 1], true, Ratio::new(1.into(), 12.into()))]
    #[case(vec![2, 2], true, Ratio::new(1.into(), 24.into()))]
    #[case(vec![3, 1], true, Ratio::default())]
    #[case(vec![3, 2], true, Ratio::new(1.into(), 180.into()))]
    #[case(vec![4, 1], true, Ratio::new((-1).into(), 720.into()))]
    #[case(vec![2, 1, 1, 1], true, Ratio::new((-1).into(), 120.into()))]
    fn test_goldberg_coeff(
        #[case] q_m: Vec<usize>,
        #[case] a_first: bool,
        #[case] expected_coeff: Ratio<i128>,
    ) {
        dbg!(&q_m);
        dbg!(a_first);
        let n: usize = q_m.iter().sum();
        let coeff = Ratio::new(
            goldberg_coeff_numerator::<i128>(q_m, a_first),
            i128::from(FACTORIALS[n] * bch_denominator::<i128>(n)),
        );
        assert_eq!(coeff, expected_coeff);
    }
}
