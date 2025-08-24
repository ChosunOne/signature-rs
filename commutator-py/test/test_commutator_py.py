import pytest
from commutator_py import CommutatorTerm, FormalIndeterminate
from lyndon_py import LyndonBasis, LyndonWord, Sort


class TestCommutatorPyGroup:
    def test_constructor(self):
        term = CommutatorTerm(1, 0)
        assert term.coefficient == 1
        assert term.atom == 0
        assert term.left is None
        assert term.right is None

        term = CommutatorTerm(1, left=CommutatorTerm(1, 0), right=CommutatorTerm(1, 1))
        assert term.coefficient == 1
        assert term.left == CommutatorTerm(1, 0)
        assert term.right == CommutatorTerm(1, 1)
        assert term.atom is None

    @pytest.mark.parametrize(
        "a,b,expected_term",
        [
            (
                CommutatorTerm(1, 0),
                CommutatorTerm(1, 1),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
            ),
            (
                CommutatorTerm(1, 1),
                CommutatorTerm(1, 0),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 1),
                    right=CommutatorTerm(1, 0),
                ),
            ),
            (
                CommutatorTerm(
                    2,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
                CommutatorTerm(
                    3,
                    left=CommutatorTerm(1, 1),
                    right=CommutatorTerm(1, 0),
                ),
                CommutatorTerm(
                    6,
                    left=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 0),
                        right=CommutatorTerm(1, 1),
                    ),
                    right=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 1),
                        right=CommutatorTerm(1, 0),
                    ),
                ),
            ),
        ],
    )
    def test_commutator_operation(
        self, a: CommutatorTerm, b: CommutatorTerm, expected_term: CommutatorTerm
    ):
        term = a.commutator(b)
        assert term == expected_term

    @pytest.mark.parametrize(
        "a,b,expected_lt",
        [
            (
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 0),
                        right=CommutatorTerm(1, 1),
                    ),
                    right=CommutatorTerm(1, 1),
                ),
                True,
            ),
            (
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 0),
                        right=CommutatorTerm(1, 1),
                    ),
                    right=CommutatorTerm(1, 1),
                ),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
                False,
            ),
        ],
    )
    def test_expression_ordering(
        self, a: CommutatorTerm, b: CommutatorTerm, expected_lt: bool
    ):
        assert (a < b) == expected_lt

    @pytest.mark.parametrize(
        "term,expected_term",
        [
            (
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 1),
                    right=CommutatorTerm(1, 0),
                ),
                CommutatorTerm(
                    -1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
            ),
            (
                CommutatorTerm(
                    -1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
                CommutatorTerm(
                    -1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
            ),
        ],
    )
    def test_lyndon_sorting(self, term: CommutatorTerm, expected_term: CommutatorTerm):
        term.lyndon_sort()
        assert term == expected_term

    @pytest.mark.parametrize(
        "term,expected_left_term,expected_right_term",
        [
            (
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 0),
                        right=CommutatorTerm(1, 1),
                    ),
                    right=CommutatorTerm(1, 2),
                ),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 1),
                        right=CommutatorTerm(1, 2),
                    ),
                ),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 0),
                        right=CommutatorTerm(1, 2),
                    ),
                    right=CommutatorTerm(1, 1),
                ),
            ),
        ],
    )
    def test_jacobi_identity(
        self,
        term: CommutatorTerm,
        expected_left_term: CommutatorTerm,
        expected_right_term: CommutatorTerm,
    ):
        ji = term.jacobi_identity()
        assert ji is not None
        left, right = ji
        assert left == expected_left_term
        assert right == expected_right_term

    @pytest.mark.parametrize(
        "word,expected_term",
        [
            (
                LyndonWord([0]),
                CommutatorTerm(1, 0),
            ),
            (
                LyndonWord([0, 1]),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
            ),
            (
                LyndonWord([0, 1, 2]),
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 1),
                        right=CommutatorTerm(1, 2),
                    ),
                ),
            ),
        ],
    )
    def test_from_lyndon_word(self, word: LyndonWord, expected_term: CommutatorTerm):
        term = CommutatorTerm.from_lyndon_word(word)
        assert term == expected_term

    @pytest.mark.parametrize(
        "num_generators,max_degree,term,expected_basis_terms",
        [
            (
                3,
                4,
                14
                * CommutatorTerm(
                    1,
                    left=CommutatorTerm(
                        1,
                        left=CommutatorTerm(
                            1,
                            left=CommutatorTerm(1, 0),
                            right=CommutatorTerm(1, 1),
                        ),
                        right=CommutatorTerm(1, 1),
                    ),
                    right=CommutatorTerm(1, 2),
                ),
                [
                    14
                    * CommutatorTerm(
                        1,
                        left=CommutatorTerm(
                            1,
                            left=CommutatorTerm(
                                1,
                                left=CommutatorTerm(1, 0),
                                right=CommutatorTerm(1, 2),
                            ),
                            right=CommutatorTerm(1, 1),
                        ),
                        right=CommutatorTerm(1, 1),
                    ),
                    28
                    * CommutatorTerm(
                        1,
                        left=CommutatorTerm(
                            1,
                            left=CommutatorTerm(1, 0),
                            right=CommutatorTerm(
                                1,
                                left=CommutatorTerm(1, 1),
                                right=CommutatorTerm(1, 2),
                            ),
                        ),
                        right=CommutatorTerm(1, 1),
                    ),
                    14
                    * CommutatorTerm(
                        1,
                        left=CommutatorTerm(1, 0),
                        right=CommutatorTerm(
                            1,
                            left=CommutatorTerm(1, 1),
                            right=CommutatorTerm(
                                1,
                                left=CommutatorTerm(1, 1),
                                right=CommutatorTerm(1, 2),
                            ),
                        ),
                    ),
                ],
            ),
        ],
    )
    def test_lyndon_basis_decomposition(
        self,
        num_generators: int,
        max_degree: int,
        term: CommutatorTerm,
        expected_basis_terms: list[CommutatorTerm],
    ):
        basis = LyndonBasis(num_generators, Sort.Lexicographical)
        basis_set = set(
            [
                CommutatorTerm.from_lyndon_word(x)
                for x in basis.generate_basis(max_degree)
            ]
        )
        basis_terms = term.lyndon_basis_decomposition(basis_set)
        expected_basis_terms.sort()
        for x in basis_terms:
            print(x)
        assert len(basis_terms) == len(expected_basis_terms)
        for basis_term, expected_basis_term in zip(basis_terms, expected_basis_terms):
            assert basis_term == expected_basis_term

    @pytest.mark.parametrize(
        "a,b,expected_indeterminates",
        [
            (
                CommutatorTerm(1, 0),
                CommutatorTerm(1, 1),
                [
                    FormalIndeterminate(1, [0, 1]),
                    FormalIndeterminate(-1, [1, 0]),
                ],
            ),
            (
                CommutatorTerm(
                    1,
                    left=CommutatorTerm(1, 0),
                    right=CommutatorTerm(1, 1),
                ),
                CommutatorTerm(
                    1, left=CommutatorTerm(1, 2), right=CommutatorTerm(1, 3)
                ),
                [
                    FormalIndeterminate(
                        1,
                        [0, 1, 2, 3],
                    ),
                    FormalIndeterminate(
                        -1,
                        [0, 1, 3, 2],
                    ),
                    FormalIndeterminate(
                        -1,
                        [1, 0, 2, 3],
                    ),
                    FormalIndeterminate(
                        1,
                        [1, 0, 3, 2],
                    ),
                    FormalIndeterminate(
                        -1,
                        [2, 3, 0, 1],
                    ),
                    FormalIndeterminate(
                        1,
                        [2, 3, 1, 0],
                    ),
                    FormalIndeterminate(
                        1,
                        [3, 2, 0, 1],
                    ),
                    FormalIndeterminate(
                        -1,
                        [3, 2, 1, 0],
                    ),
                ],
            ),
        ],
    )
    def test_formal_indeterminate_expansion(
        self,
        a: CommutatorTerm,
        b: CommutatorTerm,
        expected_indeterminates: list[FormalIndeterminate],
    ):
        term = a.commutator(b)
        indeterminates = term.formal_indeterminates()
        assert len(indeterminates) == len(expected_indeterminates)
        for indeterminate, expected_indeterminate in zip(
            indeterminates, expected_indeterminates
        ):
            assert indeterminate == expected_indeterminate
