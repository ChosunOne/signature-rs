use num_rational::Rational64;

use crate::{FACTORIALS, PRIMES};

pub fn p_adic_expansion(n: usize, p: usize) -> Vec<usize> {
    let mut num = n;
    let mut max_power = 0;
    while p.saturating_pow(max_power) <= num {
        max_power += 1;
    }
    let mut alphas = vec![0; max_power as usize];
    for (i, power) in (0..max_power).enumerate().rev() {
        alphas[i] = num / p.pow(power);
        num -= alphas[i] * p.pow(power);
    }

    alphas
}
pub fn s_p(n: usize, p: usize) -> usize {
    p_adic_expansion(n, p).iter().sum()
}

pub fn bch_denominator(n: usize) -> usize {
    let primes = PRIMES.iter().filter(|&&x| x < n);
    let mut prod = 1;
    for &p in primes {
        let s_p_n = s_p(n, p);
        let mut t = 0;
        while p.pow(t) <= s_p_n {
            t += 1;
        }
        prod *= p.pow(t - 1);
    }
    prod
}

pub fn goldberg_coeff(q: Vec<usize>, a_first: bool) -> Rational64 {
    let n = q.iter().sum();
    let d = FACTORIALS[n] * bch_denominator(n) as u64;
    let m = q.len();
    let mut c = vec![0isize; n * n];
    let mut a_current = a_first;
    if m % 2 == 0 {
        a_current = !a_first;
    }
    let mut l = 0;
    for i in (0..m).rev() {
        for r in 1..=q[i] {
            l += 1;
            let mut h = 0;
            if i == m - 1 {
                h = (d / FACTORIALS[l]) as isize;
            } else if a_current && i == m - 2 {
                h = (d / (FACTORIALS[r] * FACTORIALS[q[i + 1]])) as isize;
            }
            c[(l - 1) * n] = h;
            for k in 2..l {
                h = 0;
                for j in 1..=r {
                    if l > j && c[k - 2 + n * (l - j - 1)] != 0 {
                        h += c[k - 2 + n * (l - j - 1)] / FACTORIALS[j] as isize;
                    }
                }
                if a_current && i <= m - 2 {
                    for j in 1..=q[i + 1] {
                        if l > r + j && c[k - 2 + n * (l - r - j - 1)] != 0 {
                            h += c[k - 2 + n * (l - r - j - 1)]
                                / (FACTORIALS[r] * FACTORIALS[j]) as isize;
                        }
                    }
                }
                c[k - 1 + n * (l - 1)] = h;
            }
            c[l - 1 + n * (l - 1)] = d as isize;
        }
        a_current = !a_current;
    }
    let mut sum = Rational64::default();
    for k in 1..=n {
        let denom = if k % 2 != 0 {
            k as isize
        } else {
            -(k as isize)
        };
        sum += Rational64::new(c[k - 1 + n * (n - 1)] as i64, denom as i64);
    }
    sum * Rational64::new(1, d as i64)
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
        #[case] p: usize,
        #[case] expected_expansion: Vec<usize>,
    ) {
        let expansion = p_adic_expansion(n, p);
        assert_eq!(expansion, expected_expansion);
    }

    #[test]
    fn test_bch_denominators() {
        let n = 25;
        let denominators = (1..n + 1).map(|x| bch_denominator(x)).collect::<Vec<_>>();
        let expected_denominators = [
            1, 1, 2, 1, 6, 2, 6, 3, 10, 2, 6, 2, 210, 30, 12, 3, 30, 10, 210, 42, 330, 30, 60, 30,
            546,
        ];
        assert_eq!(denominators, expected_denominators);
    }

    #[rstest]
    #[case(vec![1], true, Rational64::new(1, 1))]
    #[case(vec![1], false, Rational64::new(1, 1))]
    #[case(vec![1, 1], true, Rational64::new(1, 2))]
    #[case(vec![2, 1], true, Rational64::new(1, 12))]
    #[case(vec![2, 2], true, Rational64::new(1,24))]
    #[case(vec![3, 1], true, Rational64::default())]
    #[case(vec![3, 2], true, Rational64::new(1, 180))]
    #[case(vec![4, 1], true, Rational64::new(-1, 720))]
    #[case(vec![2, 1, 1, 1], true, Rational64::new(-1, 120))]
    fn test_goldberg_coeff(
        #[case] q_m: Vec<usize>,
        #[case] a_first: bool,
        #[case] expected_coeff: Rational64,
    ) {
        dbg!(&q_m);
        dbg!(a_first);
        let coeff = goldberg_coeff(q_m, a_first);
        assert_eq!(coeff, expected_coeff);
    }
}
