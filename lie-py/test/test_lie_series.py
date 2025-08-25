import pytest
from commutator_py import CommutatorTerm
from lie_py import LieSeries
from lyndon_py import LyndonBasis, Sort


class TestLieSeriesGroup:
    def test_constructor(self):
        lyndon_basis = LyndonBasis(2, Sort.Lexicographical)
        basis = lyndon_basis.generate_basis(3)
        lie_series = LieSeries(basis, [0.0 for _ in basis])
        lie_basis = lie_series.basis
        for lie_word, expected_word in zip(lie_basis, basis):
            assert lie_word == expected_word
        commutator_basis = [CommutatorTerm.from_lyndon_word(w) for w in basis]
        lie_commutator_basis = lie_series.commutator_basis
        for lie_term, expected_term in zip(lie_commutator_basis, commutator_basis):
            assert lie_term == expected_term

        lie_coefficients = lie_series.coefficients
        for c in lie_coefficients:
            assert c == 0.0

        assert lie_series.max_degree == 3

    @pytest.mark.parametrize(
        "num_generators,basis_depth,a_coefficients,b_coefficients,expected_coefficients",
        [
            (
                2,
                2,
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [0.0, 0.0, -3.0],
            ),
            (
                2,
                2,
                [3.0, 2.0, 1.0],
                [1.0, 2.0, 3.0],
                [0.0, 0.0, 4.0],
            ),
            (
                2,
                3,
                [1.0, 2.0, 3.0, 4.0, 5.0],
                [6.0, 7.0, 8.0, 9.0, 10.0],
                [0.0, 0.0, -5.0, -10.0, 5.0],
            ),
            (
                3,
                3,
                [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [5.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [
                    0.0,
                    0.0,
                    0.0,
                    -7.0,
                    -14.0,
                    -7.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                3,
                3,
                [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [
                    0.0,
                    0.0,
                    0.0,
                    -7.0,
                    -14.0,
                    -7.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -7.0,
                    -14.0,
                    14.0,
                    14.0,
                    49.0,
                    42.0,
                    -14.0,
                    21.0,
                ],
            ),
        ],
    )
    def test_commutator(
        self,
        num_generators: int,
        basis_depth: int,
        a_coefficients: list[float],
        b_coefficients: list[float],
        expected_coefficients: list[float],
    ):
        basis = LyndonBasis(num_generators, Sort.Lexicographical).generate_basis(
            basis_depth
        )
        a = LieSeries(basis, a_coefficients)
        b = LieSeries(basis, b_coefficients)
        series = a.commutator(b)
        coefficients = series.coefficients
        assert len(coefficients) == len(expected_coefficients)
        for c, e_c in zip(coefficients, expected_coefficients):
            assert (c - e_c).__abs__() < 0.001, f"{c} != {e_c}"
