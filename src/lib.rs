use std::{
    collections::HashMap,
    ops::{Index, IndexMut},
};

use bitvec::prelude::*;
use num_rational::Rational64;

use crate::{
    lyndon::{Generator, LyndonBasis},
    rooted_tree::RootedTree,
};
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

pub struct BCHCoefficientGenerator<T: Generator<Letter = T>> {
    graph_partition_table: GraphPartitionTable<T>,
    is_computed_z: BitVec,
    z: Vec<Rational64>,
    bernoulli: Vec<Rational64>,
    x_minus_y: Vec<Rational64>,
    x_plus_y: Vec<Rational64>,
    prime: Vec<usize>,
    d_prime: Vec<usize>,
    sigma: Vec<usize>,
    kappa: Vec<usize>,
    m_n: usize,
    adjoint_cache: HashMap<(usize, usize), Rational64>,
}

impl<T: Generator<Letter = T>> BCHCoefficientGenerator<T> {
    pub fn new(n: usize) -> Self {
        let t_n = LyndonBasis::<2, T>::generate_basis(n)
            .into_iter()
            .map(|w| RootedTree::from(w))
            .collect::<Vec<_>>();
        let m_n = t_n.len();
        let graph_partition_table = GraphPartitionTable::new(t_n);
        let is_computed_z = bitvec![0; graph_partition_table.tm_n()];
        let bernoulli = bernoulli_sequence(n);
        let z = vec![Rational64::default(); graph_partition_table.tm_n()];
        let mut x_minus_y = vec![Rational64::default(); graph_partition_table.tm_n()];
        x_minus_y[0] = Rational64::new(1, 1);
        x_minus_y[1] = Rational64::new(-1, 1);
        let mut x_plus_y = vec![Rational64::default(); graph_partition_table.tm_n()];
        x_plus_y[0] = Rational64::new(1, 1);
        x_plus_y[1] = Rational64::new(1, 1);
        let prime = vec![0; m_n];
        let d_prime = vec![0; m_n];
        let mut sigma = vec![0_usize; m_n];
        sigma[0] = 1;
        sigma[1] = 1;
        let mut kappa = vec![0; m_n];
        kappa[0] = 1;
        kappa[1] = 1;
        let adjoint_cache = HashMap::new();

        Self {
            graph_partition_table,
            is_computed_z,
            bernoulli,
            z,
            x_minus_y,
            x_plus_y,
            prime,
            d_prime,
            sigma,
            kappa,
            m_n,
            adjoint_cache,
        }
    }

    fn lie_bracket(&self, alpha: &[Rational64], beta: &[Rational64], i: usize) -> Rational64 {
        let mut sum = Rational64::default();
        for &(j, k) in self.graph_partition_table.partitions(i).partitions.iter() {
            sum += alpha[j] * beta[k] - alpha[k] * beta[j];
        }

        sum
    }

    fn adjoint_operator(&mut self, i: usize, power: usize) -> Rational64 {
        if power == 0 {
            return self.x_plus_y[i];
        }
        if let Some(&v) = self.adjoint_cache.get(&(i, power)) {
            return v;
        }
        if power == 1 {
            let result = self.lie_bracket(&self.z, &self.x_plus_y, i);
            self.adjoint_cache.insert((i, power), result);
            return result;
        }
        let mut sum = Rational64::default();
        let partitions = self
            .graph_partition_table
            .partitions(i)
            .partitions
            .iter()
            .copied()
            .collect::<Vec<_>>();
        for (j, k) in partitions {
            sum += self.z[j] * self.adjoint_operator(k, power - 1)
                - self.z[k] * self.adjoint_operator(j, power - 1);
        }
        self.adjoint_cache.insert((i, power), sum);

        sum
    }

    fn compute_z(&mut self, i: usize) -> Rational64 {
        if i == 0 || i == 1 {
            self.is_computed_z.set(i, true);
            return Rational64::new(1, 1);
        }
        // 1/2 [X-Y, Z] (u_i)
        let s_ui = self.graph_partition_table.partitions(i).clone();
        for &(j, k) in &s_ui.partitions {
            // Check if there are any Z_u terms that we need to complete first
            if !self.is_computed_z[j] {
                self.z[j] = self.compute_z(j);
            }
            if !self.is_computed_z[k] {
                self.z[k] = self.compute_z(k);
            }
        }
        let left_term = Rational64::new(1, 2) * self.lie_bracket(&self.x_minus_y, &self.z, i);

        // Bernoulli term
        let mut bernoulli_term = Rational64::default();
        let degree = self.graph_partition_table.degree(i);
        for p in 1..=((degree - 1) / 2) {
            let bernoulli_coef =
                self.bernoulli[2 * p] * Rational64::new(1, FACTORIALS[2 * p] as i64);
            let adjoint_term = self.adjoint_operator(i, 2 * p);
            bernoulli_term += bernoulli_coef * adjoint_term;
        }
        self.is_computed_z.set(i, true);
        let result = (left_term + bernoulli_term)
            * Rational64::new(1, self.graph_partition_table.degree(i) as i64);
        result
    }

    /// Generates the truncated Baker-Campbell-Hausdorff series of a Lyndon basis
    pub fn generate_coefficients(&mut self) -> Vec<Rational64> {
        for i in 0..self.m_n {
            if !self.is_computed_z[i] {
                self.z[i] = self.compute_z(i);
                if i > 1 {
                    self.prime[i] = self.graph_partition_table.partitions(i)[0].0;
                    self.d_prime[i] = self.graph_partition_table.partitions(i)[0].1;
                    self.kappa[i] = {
                        if self.d_prime[self.prime[i]] != self.d_prime[i] {
                            1
                        } else {
                            self.kappa[self.prime[i]] + 1
                        }
                    };
                    self.sigma[i] =
                        self.kappa[i] * self.sigma[self.prime[i]] * self.sigma[self.d_prime[i]];
                }
            }
        }

        self.z[..self.m_n]
            .iter()
            .zip(self.sigma.iter())
            .map(|(&z, &sig)| z * Rational64::new(1, sig as i64))
            .collect()
    }
}

#[derive(Clone, Debug, Default)]
pub struct EdgePartitions {
    pub partitions: Vec<(usize, usize)>,
}

impl Index<usize> for EdgePartitions {
    type Output = (usize, usize);

    fn index(&self, index: usize) -> &Self::Output {
        &self.partitions[index]
    }
}

impl IndexMut<usize> for EdgePartitions {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.partitions[index]
    }
}

#[derive(Debug, Clone)]
pub struct GraphPartitionTable<T: Generator<Letter = T>> {
    t_n: Vec<RootedTree<T>>,
    degree: Vec<usize>,
    s: Vec<EdgePartitions>,
}

impl<T: Generator<Letter = T>> GraphPartitionTable<T> {
    pub fn new(mut t_n: Vec<RootedTree<T>>) -> Self {
        let mut degree = Vec::with_capacity(t_n.len());
        let mut tree_t_n_map = HashMap::<RootedTree<T>, usize>::new();
        for (i, tree) in t_n.iter().enumerate() {
            degree.push(tree.degree());
            tree_t_n_map.insert(tree.clone(), i);
        }
        let mut s = vec![EdgePartitions::default(); t_n.len()];
        let mut i = 0;
        while i < t_n.len() {
            let tree = &t_n[i];
            let Some((v, w)) = tree.factorize() else {
                i += 1;
                continue;
            };
            let v_idx = tree_t_n_map[&v];
            let w_idx = tree_t_n_map[&w];
            s[i].partitions.push((v_idx, w_idx));

            for p in 0..v.degree() - 1 {
                let s_v = &s[v_idx];
                let mut v_root_w = t_n[s_v[p].0].clone();
                v_root_w.graft(w.clone());
                let v_comp = t_n[s_v[p].1].clone();
                if !tree_t_n_map.contains_key(&v_root_w) {
                    t_n.push(v_root_w.clone());
                    degree.push(v_root_w.degree());
                    tree_t_n_map.insert(v_root_w.clone(), t_n.len() - 1);
                    s.push(EdgePartitions::default());
                }
                s[i].partitions
                    .push((tree_t_n_map[&v_root_w], tree_t_n_map[&v_comp]));
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
                    s.push(EdgePartitions::default());
                }
                s[i].partitions
                    .push((tree_t_n_map[&v_w_root], tree_t_n_map[&w_comp]));
            }
            i += 1;
        }

        Self { degree, t_n, s }
    }

    pub fn partitions(&self, i: usize) -> &EdgePartitions {
        &self.s[i]
    }

    pub fn degree(&self, i: usize) -> usize {
        self.degree[i]
    }

    pub fn tree(&self, i: usize) -> &RootedTree<T> {
        &self.t_n[i]
    }

    pub fn tm_n(&self) -> usize {
        self.t_n.len()
    }
}

fn binomial(n: usize, k: usize) -> Rational64 {
    if k == 0 {
        return Rational64::new(1, 1);
    }

    let k = k.min(n - k);
    let mut result = Rational64::new(1, 1);

    for i in 0..k {
        result *= Rational64::new((n - i) as i64, (i + 1) as i64);
    }

    result
}

pub fn bernoulli_sequence(max_n: usize) -> Vec<Rational64> {
    let mut b = vec![Rational64::default(); max_n + 1];

    b[0] = Rational64::new(1, 1);
    if max_n > 0 {
        b[1] = Rational64::new(-1, 2);
    }

    for n in 2..=max_n {
        if n % 2 == 1 {
            continue;
        }
        let mut sum = Rational64::default();
        for k in 0..n {
            sum += binomial(n + 1, k) * b[k];
        }

        b[n] = -sum * Rational64::new(1, (n + 1) as i64);
    }

    // Use positive sign convention
    if max_n > 0 {
        b[1] = Rational64::new(1, 2);
    }

    b
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bernoulli_sequence() {
        let seq = bernoulli_sequence(10);
        let expected_seq = vec![
            Rational64::new(1, 1),
            Rational64::new(1, 2),
            Rational64::new(1, 6),
            Rational64::default(),
            Rational64::new(-1, 30),
            Rational64::default(),
            Rational64::new(1, 42),
            Rational64::default(),
            Rational64::new(-1, 30),
            Rational64::default(),
            Rational64::new(5, 66),
        ];
        assert_eq!(seq.len(), expected_seq.len());
        for (&term, &expected_term) in seq.iter().zip(&expected_seq) {
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_graph_partition_table() {
        let t_n = LyndonBasis::<2, char>::generate_basis(5)
            .into_iter()
            .map(|w| RootedTree::from(w))
            .collect::<Vec<_>>();
        let m_n = t_n.len();
        let graph_partition_table = GraphPartitionTable::new(t_n);
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
        for (i, (s_ui, expected_s_ui)) in graph_partition_table
            .s
            .iter()
            .zip(expected_s.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(&s_ui.partitions, expected_s_ui);
        }
    }

    #[test]
    fn test_bch_series() {
        let mut generator = BCHCoefficientGenerator::<char>::new(5);
        let series = generator.generate_coefficients();
        let expected_z_ui_series = vec![
            Rational64::new(1, 1),
            Rational64::new(1, 1),
            Rational64::new(1, 2),
            Rational64::new(1, 12),
            Rational64::new(1, 12),
            Rational64::default(),
            Rational64::new(1, 24),
            Rational64::default(),
            Rational64::new(-1, 720),
            Rational64::new(1, 360),
            Rational64::new(1, 120),
            Rational64::new(1, 180),
            Rational64::new(1, 180),
            Rational64::new(-1, 720),
        ];
        for (&term, &expected_term) in series.iter().zip(&expected_z_ui_series) {
            assert_eq!(term, expected_term);
        }
    }
}
