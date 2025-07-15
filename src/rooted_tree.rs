use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Index, IndexMut},
};

use crate::lyndon::{Generator, LyndonWord, LyndonWordError};

#[derive(Default, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RootedTree<T: Generator<Letter = T>> {
    color: T,
    children: Vec<RootedTree<T>>,
    degree: usize,
}

impl<T: Generator<Letter = T>> RootedTree<T> {
    pub fn new(color: T) -> Self {
        Self {
            color,
            children: Vec::new(),
            degree: 1,
        }
    }

    pub fn graft(&mut self, tree: RootedTree<T>) {
        self.degree += tree.degree;
        self.children.push(tree);
        self.canonicalize();
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn factorize(&self) -> Option<(RootedTree<T>, RootedTree<T>)> {
        if self.degree == 1 {
            return None;
        }

        let mut children = self.children.clone();

        let w = children
            .pop()
            .expect("There to be a child for a tree with degree > 1");
        let mut v = Self::new(self.color);
        for child in children {
            v.graft(child);
        }

        Some((v, w))
    }

    fn get_node(&self, path: &[usize]) -> &Self {
        let mut current = self;
        for &index in path {
            current = &current.children[index];
        }
        current
    }

    fn get_node_mut(&mut self, path: &[usize]) -> &mut Self {
        let mut current = self;
        for &index in path {
            current = &mut current.children[index];
        }
        current
    }

    fn to_letters(&self) -> Vec<T> {
        let mut letters = vec![self.color];

        for child in &self.children {
            letters.extend(child.to_letters());
        }

        letters
    }

    /// Sort the tree in canonical order, which is by color then by children
    fn canonicalize(&mut self) {
        for child in self.children.iter_mut() {
            child.canonicalize();
        }

        self.children.sort_by_key(|x| (x.degree, x.color));
    }

    fn count_degrees(&mut self) -> usize {
        let mut degree = 1;

        for child in &mut self.children {
            degree += child.count_degrees();
        }

        self.degree = degree;
        degree
    }
}

impl<const N: usize, T: Generator<Letter = T>> From<LyndonWord<N, T>> for RootedTree<T> {
    fn from(value: LyndonWord<N, T>) -> Self {
        if value.len() == 1 {
            return Self::new(value.letters[0]);
        }

        let (v, w) = value.factorize();
        let mut root = RootedTree::from(v);
        let comp_tree = RootedTree::from(w);
        root.graft(comp_tree);

        root.count_degrees();
        root.canonicalize();

        root
    }
}

impl<const N: usize, T: Generator<Letter = T>> TryFrom<&RootedTree<T>> for LyndonWord<N, T> {
    type Error = LyndonWordError;
    fn try_from(value: &RootedTree<T>) -> Result<Self, Self::Error> {
        let letters = value.to_letters();
        let candidate_word = LyndonWord::try_from(letters)?;
        let candidate_tree = RootedTree::from(candidate_word.clone());
        if candidate_tree == *value {
            return Ok(candidate_word);
        }
        Err(LyndonWordError::InvalidWord)
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
        #[cfg(feature = "progress")]
        let mut pb = ProgressBar::new(t_n.len() as u64);
        #[cfg(feature = "progress")]
        {
            pb.set_message("Generating Partition Table ");
            pb.tick();
        }

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
                #[cfg(feature = "progress")]
                pb.inc(1);

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
                    #[cfg(feature = "progress")]
                    {
                        pb.inc_length(1);
                        pb.inc(1);
                    }
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
                    #[cfg(feature = "progress")]
                    {
                        pb.inc_length(1);
                        pb.inc(1);
                    }
                    t_n.push(v_w_root.clone());
                    degree.push(v_w_root.degree());
                    tree_t_n_map.insert(v_w_root.clone(), t_n.len() - 1);
                    s.push(EdgePartitions::default());
                }
                s[i].partitions
                    .push((tree_t_n_map[&v_w_root], tree_t_n_map[&w_comp]));
            }
            i += 1;
            #[cfg(feature = "progress")]
            pb.inc(1);
        }

        #[cfg(feature = "progress")]
        pb.finish();

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

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use rstest::rstest;

    use crate::lyndon::LyndonBasis;

    use super::*;

    #[rstest]
    #[case("X", vec![])]
    #[case("Y", vec![])]
    #[case("XY", vec![vec![0]])]
    #[case("XYY", vec![vec![0], vec![1]])]
    #[case("XXY", vec![vec![0], vec![0, 0]])]
    #[case("XYYY", vec![vec![0], vec![1], vec![2]])]
    #[case("XXYY", vec![vec![0], vec![0, 0], vec![0, 1]])]
    #[case("XXXY", vec![vec![0], vec![0, 0], vec![0, 0, 0]])]
    #[case("XYYYY", vec![vec![0], vec![1], vec![2], vec![3]])]
    #[case("XXYXY", vec![vec![0], vec![0, 0], vec![1], vec![1, 0]])]
    #[case("XYXYY", vec![vec![0], vec![1], vec![1, 0], vec![1, 1]])]
    #[case("XXYYY", vec![vec![0], vec![0, 0], vec![0, 1], vec![0, 2]])]
    #[case("XXXYY", vec![vec![0], vec![0, 0], vec![0, 0, 0], vec![0, 0, 1]])]
    #[case("XXXXY", vec![vec![0], vec![0, 0], vec![0, 0, 0], vec![0, 0, 0, 0]])]
    fn it_creates_a_rooted_tree_from_a_lyndon_word(
        #[case] letters: &str,
        #[case] paths: Vec<Vec<usize>>,
    ) -> Result<(), LyndonWordError> {
        let word = letters.parse::<LyndonWord<2, char>>()?;
        dbg!(&word);
        let tree = RootedTree::from(word);
        dbg!(&tree);
        assert_eq!(tree.degree(), letters.len());
        assert_eq!(tree.color, letters.chars().collect::<Vec<_>>()[0]);
        for (path, expected_color) in paths.iter().zip(letters.chars().skip(1)) {
            assert_eq!(tree.get_node(path).color, expected_color);
        }

        Ok(())
    }

    #[rstest]
    #[case("X")]
    #[case("Y")]
    #[case("XY")]
    #[case("XYY")]
    #[case("XXY")]
    #[case("XYYY")]
    #[case("XXYY")]
    #[case("XXXY")]
    #[case("XYYYY")]
    #[case("XXYXY")]
    #[case("XYXYY")]
    #[case("XXYYY")]
    #[case("XXXYY")]
    #[case("XXXXY")]
    fn it_converts_a_rooted_tree_back_to_a_lyndon_word(
        #[case] letters: &str,
    ) -> Result<(), LyndonWordError> {
        let word = letters.parse::<LyndonWord<2, char>>()?;
        let tree = RootedTree::from(word.clone());
        let reconstructed_word = LyndonWord::try_from(&tree)?;
        assert_eq!(word, reconstructed_word);
        Ok(())
    }

    #[rstest]
    #[case("X", "Y")]
    #[case("XY", "Y")]
    #[case("X", "XY")]
    #[case("XYY", "Y")]
    #[case("X", "XYY")]
    #[case("X", "XXY")]
    #[case("XYYY", "Y")]
    #[case("XXY", "XY")]
    #[case("XY", "XYY")]
    #[case("X", "XXYY")]
    #[case("X", "XXXY")]
    fn it_grafts_a_rooted_tree(
        #[case] word_1_letters: &str,
        #[case] word_2_letters: &str,
    ) -> Result<(), LyndonWordError> {
        let word_1 = word_1_letters.parse::<LyndonWord<2, char>>()?;
        dbg!(&word_1);
        let word_2 = word_2_letters.parse::<LyndonWord<2, char>>()?;
        dbg!(&word_2);
        let mut tree_1 = RootedTree::from(word_1);
        dbg!(&tree_1);
        let tree_2 = RootedTree::from(word_2);
        dbg!(&tree_2);
        tree_1.graft(tree_2);
        dbg!(&tree_1);
        let reconstructed_word = LyndonWord::try_from(&tree_1)?;
        let expected_word = (word_1_letters.to_owned() + word_2_letters)
            .as_str()
            .parse::<LyndonWord<2, char>>()?;
        assert_eq!(reconstructed_word, expected_word);
        Ok(())
    }

    #[rstest]
    #[case("XY", "X", "Y")]
    #[case("XYY", "XY", "Y")]
    #[case("XXY", "X", "XY")]
    #[case("XYYY", "XYY", "Y")]
    #[case("XXYY", "X", "XYY")]
    #[case("XXXY", "X", "XXY")]
    #[case("XYYYY", "XYYY", "Y")]
    #[case("XXYXY", "XXY", "XY")]
    #[case("XYXYY", "XY", "XYY")]
    #[case("XXXYY", "X", "XXYY")]
    #[case("XXXXY", "X", "XXXY")]
    fn it_factorizes_a_rooted_tree_1(
        #[case] letters: &str,
        #[case] expected_v: &str,
        #[case] expected_w: &str,
    ) -> Result<(), LyndonWordError> {
        let word = LyndonWord::<2, char>::try_from(letters.chars().collect::<Vec<_>>())?;
        let tree = RootedTree::from(word);
        let (v, w) = tree.factorize().expect("To factorize a lyndon tree.");
        let expected_v_word =
            LyndonWord::<2, char>::try_from(expected_v.chars().collect::<Vec<_>>())?;
        let expected_w_word =
            LyndonWord::<2, char>::try_from(expected_w.chars().collect::<Vec<_>>())?;
        let v_word = LyndonWord::try_from(&v)?;
        let w_word = LyndonWord::try_from(&w)?;

        assert_eq!(v_word, expected_v_word);
        assert_eq!(w_word, expected_w_word);
        Ok(())
    }

    #[test]
    fn it_factorizes_a_non_lyndon_rooted_tree() {
        let letters = vec!['X', 'Y'];
        let word = LyndonWord::<2, char>::try_from(letters).expect("To make a lyndon word");
        let mut tree = RootedTree::from(word.clone());
        tree.graft(RootedTree::from(word));
        let Some((v, w)) = tree.factorize() else {
            panic!("Failed to factorize tree");
        };
        assert_eq!(v.degree(), 2);
        assert_eq!(v.color, 'X');
        assert_eq!(v.get_node(&[0]).color, 'Y');
        assert_eq!(w.degree(), 2);
        assert_eq!(w.color, 'X');
        assert_eq!(w.get_node(&[0]).color, 'Y');
    }

    #[test]
    fn it_correctly_identifies_isomorphisms() {
        let mut tree_set = HashSet::new();

        let mut tree_1 = RootedTree::new('X');
        tree_1.graft(RootedTree::new('X'));

        // (X) <- (X) -> (X) -> (Y)
        let mut grafted_tree_1 = tree_1.clone();
        grafted_tree_1.graft(RootedTree::from(
            LyndonWord::<2, char>::try_from(vec!['X', 'Y']).expect("To make a lyndon word"),
        ));
        tree_set.insert(grafted_tree_1);

        // (Y) <- (X) <- (X) -> (X)
        let mut grafted_tree_2 = RootedTree::new('X');
        grafted_tree_2.graft(RootedTree::from(
            LyndonWord::<2, char>::try_from(vec!['X', 'Y']).expect("To make a lyndon word"),
        ));
        grafted_tree_2.graft(RootedTree::new('X'));

        assert!(tree_set.contains(&grafted_tree_2));
    }

    #[test]
    fn test_graph_partition_table() {
        let t_n = LyndonBasis::<2, char>::generate_basis(5, true)
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
}
