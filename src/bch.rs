use num_rational::Ratio;

use crate::{FACTORIALS, Int, PRIMES};

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

pub fn bch_denominator(n: usize) -> i128 {
    let primes = PRIMES.iter().filter(|&&x| x < n as i128);
    let mut prod = 1i128;
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

pub fn goldberg_coeff<U: Int>(q: Vec<usize>, a_first: bool) -> Ratio<U> {
    let n = q.iter().sum();
    let d = FACTORIALS[n] * bch_denominator(n);
    let m = q.len();
    let mut c = vec![Ratio::<U>::default(); n * n];
    let mut a_current = a_first;
    if m % 2 == 0 {
        a_current = !a_first;
    }
    let mut l = 0;
    for i in (0..m).rev() {
        for r in 1..=q[i] {
            l += 1;
            let mut h = Ratio::<U>::default();
            if i == m - 1 {
                h = Ratio::new(d.into(), FACTORIALS[l].into());
            } else if a_current && i == m - 2 {
                h = Ratio::new(d.into(), (FACTORIALS[r] * FACTORIALS[q[i + 1]]).into());
            }
            c[(l - 1) * n] = h;
            for k in 2..l {
                h = Ratio::from_integer(0.into());
                for j in 1..=r {
                    if l > j && c[k - 2 + n * (l - j - 1)] != Ratio::from_integer(0.into()) {
                        h += &c[k - 2 + n * (l - j - 1)]
                            * Ratio::<U>::new(1.into(), FACTORIALS[j].into());
                    }
                }
                if a_current && i <= m - 2 {
                    for j in 1..=q[i + 1] {
                        if l > r + j
                            && c[k - 2 + n * (l - r - j - 1)] != Ratio::from_integer(0.into())
                        {
                            h += &c[k - 2 + n * (l - r - j - 1)]
                                * Ratio::<U>::new(1.into(), (FACTORIALS[r] * FACTORIALS[j]).into());
                        }
                    }
                }
                c[k - 1 + n * (l - 1)] = h;
            }
            c[l - 1 + n * (l - 1)] = Ratio::<U>::new(d.into(), 1.into());
        }
        a_current = !a_current;
    }
    let mut sum = Ratio::<U>::default();
    for k in 1..=n {
        let denom: U = if k % 2 != 0 {
            (k as u64).into()
        } else {
            (-(k as i64)).into()
        };
        sum += &c[k - 1 + n * (n - 1)] * Ratio::new(1.into(), denom);
    }
    sum * Ratio::<U>::new(1.into(), d.into())
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
        let denominators = (1..n + 1).map(|x| bch_denominator(x)).collect::<Vec<_>>();
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
        let coeff = goldberg_coeff(q_m, a_first);
        assert_eq!(coeff, expected_coeff);
    }
}
