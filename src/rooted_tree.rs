use std::{fmt::Display, hash::Hash};

use crate::lyndon_word::{LyndonWord, LyndonWordError};

#[derive(Default, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RootedTree<T: Ord + Copy + Clone + Display + Eq + Hash> {
    color: T,
    children: Vec<RootedTree<T>>,
    degree: usize,
}

impl<T: Ord + Copy + Clone + Display + Eq + Hash> RootedTree<T> {
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

impl<T: Ord + Copy + Clone + Display + Hash> From<LyndonWord<T>> for RootedTree<T> {
    fn from(value: LyndonWord<T>) -> Self {
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

impl<T: Ord + Copy + Clone + Display + Eq + Hash> TryFrom<&RootedTree<T>> for LyndonWord<T> {
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

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case("x", vec![])]
    #[case("y", vec![])]
    #[case("xy", vec![vec![0]])]
    #[case("xyy", vec![vec![0], vec![1]])]
    #[case("xxy", vec![vec![0], vec![0, 0]])]
    #[case("xyyy", vec![vec![0], vec![1], vec![2]])]
    #[case("xxyy", vec![vec![0], vec![0, 0], vec![0, 1]])]
    #[case("xxxy", vec![vec![0], vec![0, 0], vec![0, 0, 0]])]
    #[case("xyyyy", vec![vec![0], vec![1], vec![2], vec![3]])]
    #[case("xxyxy", vec![vec![0], vec![0, 0], vec![1], vec![1, 0]])]
    #[case("xyxyy", vec![vec![0], vec![1], vec![1, 0], vec![1, 1]])]
    #[case("xxyyy", vec![vec![0], vec![0, 0], vec![0, 1], vec![0, 2]])]
    #[case("xxxyy", vec![vec![0], vec![0, 0], vec![0, 0, 0], vec![0, 0, 1]])]
    #[case("xxxxy", vec![vec![0], vec![0, 0], vec![0, 0, 0], vec![0, 0, 0, 0]])]
    fn it_creates_a_rooted_tree_from_a_lyndon_word(
        #[case] letters: &str,
        #[case] paths: Vec<Vec<usize>>,
    ) -> Result<(), LyndonWordError> {
        let word = letters.parse::<LyndonWord<char>>()?;
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
    #[case("x")]
    #[case("y")]
    #[case("xy")]
    #[case("xyy")]
    #[case("xxy")]
    #[case("xyyy")]
    #[case("xxyy")]
    #[case("xxxy")]
    #[case("xyyyy")]
    #[case("xxyxy")]
    #[case("xyxyy")]
    #[case("xxyyy")]
    #[case("xxxyy")]
    #[case("xxxxy")]
    fn it_converts_a_rooted_tree_back_to_a_lyndon_word(
        #[case] letters: &str,
    ) -> Result<(), LyndonWordError> {
        let word = letters.parse::<LyndonWord<char>>()?;
        let tree = RootedTree::from(word.clone());
        let reconstructed_word = LyndonWord::try_from(&tree)?;
        assert_eq!(word, reconstructed_word);
        Ok(())
    }

    #[rstest]
    #[case("x", "y")]
    #[case("xy", "y")]
    #[case("x", "xy")]
    #[case("xyy", "y")]
    #[case("x", "xyy")]
    #[case("x", "xxy")]
    #[case("xyyy", "y")]
    #[case("xxy", "xy")]
    #[case("xy", "xyy")]
    #[case("x", "xxyy")]
    #[case("x", "xxxy")]
    fn it_grafts_a_rooted_tree(
        #[case] word_1_letters: &str,
        #[case] word_2_letters: &str,
    ) -> Result<(), LyndonWordError> {
        let word_1 = word_1_letters.parse::<LyndonWord<char>>()?;
        dbg!(&word_1);
        let word_2 = word_2_letters.parse::<LyndonWord<char>>()?;
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
            .parse::<LyndonWord<char>>()?;
        assert_eq!(reconstructed_word, expected_word);
        Ok(())
    }

    #[rstest]
    #[case("xy", "x", "y")]
    #[case("xyy", "xy", "y")]
    #[case("xxy", "x", "xy")]
    #[case("xyyy", "xyy", "y")]
    #[case("xxyy", "x", "xyy")]
    #[case("xxxy", "x", "xxy")]
    #[case("xyyyy", "xyyy", "y")]
    #[case("xxyxy", "xxy", "xy")]
    #[case("xyxyy", "xy", "xyy")]
    #[case("xxxyy", "x", "xxyy")]
    #[case("xxxxy", "x", "xxxy")]
    fn it_factorizes_a_rooted_tree_1(
        #[case] letters: &str,
        #[case] expected_v: &str,
        #[case] expected_w: &str,
    ) -> Result<(), LyndonWordError> {
        let word = LyndonWord::try_from(letters.chars().collect::<Vec<_>>())?;
        let tree = RootedTree::from(word);
        let (v, w) = tree.factorize().expect("To factorize a lyndon tree.");
        let expected_v_word = LyndonWord::try_from(expected_v.chars().collect::<Vec<_>>())?;
        let expected_w_word = LyndonWord::try_from(expected_w.chars().collect::<Vec<_>>())?;
        let v_word = LyndonWord::try_from(&v)?;
        let w_word = LyndonWord::try_from(&w)?;

        assert_eq!(v_word, expected_v_word);
        assert_eq!(w_word, expected_w_word);
        Ok(())
    }

    #[test]
    fn it_factorizes_a_non_lyndon_rooted_tree() {
        let letters = vec!['x', 'y'];
        let word = LyndonWord::try_from(letters).expect("To make a lyndon word");
        let mut tree = RootedTree::from(word.clone());
        tree.graft(RootedTree::from(word));
        let Some((v, w)) = tree.factorize() else {
            panic!("Failed to factorize tree");
        };
        assert_eq!(v.degree(), 2);
        assert_eq!(v.color, 'x');
        assert_eq!(v.get_node(&[0]).color, 'y');
        assert_eq!(w.degree(), 2);
        assert_eq!(w.color, 'x');
        assert_eq!(w.get_node(&[0]).color, 'y');
    }

    #[test]
    fn it_correctly_identifies_isomorphisms() {
        let mut tree_set = HashSet::new();

        let mut tree_1 = RootedTree::new('x');
        tree_1.graft(RootedTree::new('x'));

        // (X) <- (X) -> (X) -> (Y)
        let mut grafted_tree_1 = tree_1.clone();
        grafted_tree_1.graft(RootedTree::from(
            LyndonWord::try_from(vec!['x', 'y']).expect("To make a lyndon word"),
        ));
        tree_set.insert(grafted_tree_1);

        // (Y) <- (X) <- (X) -> (X)
        let mut grafted_tree_2 = RootedTree::new('x');
        grafted_tree_2.graft(RootedTree::from(
            LyndonWord::try_from(vec!['x', 'y']).expect("To make a lyndon word"),
        ));
        grafted_tree_2.graft(RootedTree::new('x'));

        assert!(tree_set.contains(&grafted_tree_2));
    }
}
