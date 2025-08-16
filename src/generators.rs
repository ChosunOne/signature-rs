use std::{
    fmt::{Debug, Display},
    hash::Hash,
};

pub trait Generator:
    Debug + Copy + PartialEq + PartialOrd + Ord + Clone + Display + Eq + Hash
{
    fn alphabet(size: usize) -> Vec<Self>;
}

impl Generator for u8 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as u8;
        }

        letters
    }
}

impl Generator for u16 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as u16;
        }

        letters
    }
}

impl Generator for u32 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as u32;
        }

        letters
    }
}

impl Generator for u64 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as u64;
        }

        letters
    }
}

impl Generator for usize {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i;
        }

        letters
    }
}

impl Generator for u128 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as u128;
        }

        letters
    }
}

impl Generator for i8 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as i8;
        }

        letters
    }
}

impl Generator for i16 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as i16;
        }

        letters
    }
}

impl Generator for i32 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as i32;
        }

        letters
    }
}

impl Generator for i64 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as i64;
        }

        letters
    }
}

impl Generator for isize {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as isize;
        }

        letters
    }
}

impl Generator for i128 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for (i, l) in letters.iter_mut().enumerate() {
            *l = i as i128;
        }

        letters
    }
}

impl Generator for char {
    fn alphabet(size: usize) -> Vec<Self> {
        assert!(
            (size <= 26),
            "Only up to 26 generators are supported for 'char' based generators."
        );

        let alphabet_letters = [
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
            'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        ];

        let mut letters = vec!['A'; size];
        letters.copy_from_slice(&alphabet_letters[..size]);

        letters
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ENotation {
    symbol: [char; 2],
}

impl Display for ENotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.symbol.iter().collect::<String>())
    }
}

impl Generator for ENotation {
    fn alphabet(size: usize) -> Vec<Self> {
        assert!(
            size <= 9,
            "Only up to 9 generators are supported for `ENotation` based generators."
        );
        let letters = vec![
            ENotation {
                symbol: ['e', '₁']
            },
            ENotation {
                symbol: ['e', '₂']
            },
            ENotation {
                symbol: ['e', '₃']
            },
            ENotation {
                symbol: ['e', '₄']
            },
            ENotation {
                symbol: ['e', '₅']
            },
            ENotation {
                symbol: ['e', '₆']
            },
            ENotation {
                symbol: ['e', '₇']
            },
            ENotation {
                symbol: ['e', '₈']
            },
            ENotation {
                symbol: ['e', '₉']
            },
        ];

        letters
    }
}
