from lie_py import BCHSeriesGenerator
from lyndon_py import LyndonBasis, LyndonWord, Sort


class TestBCHSeriesGeneratorGroup:
    def test_constructor(self):
        basis = LyndonBasis(2, Sort.Lexicographical)
        bch_series_generator = BCHSeriesGenerator(basis, 3)

        assert bch_series_generator.alphabet_size == 2

        basis = bch_series_generator.basis
        assert basis[0] == LyndonWord([0])
        assert basis[1] == LyndonWord([1])
        assert basis[2] == LyndonWord([0, 1])
        assert basis[3] == LyndonWord([0, 0, 1])
        assert basis[4] == LyndonWord([0, 1, 1])

        left_factor = bch_series_generator.left_factor
        assert left_factor[0] == 0
        assert left_factor[1] == 1
        assert left_factor[2] == 0
        assert left_factor[3] == 0
        assert left_factor[4] == 2

        right_factor = bch_series_generator.right_factor
        assert right_factor[0] == 0
        assert right_factor[1] == 0
        assert right_factor[2] == 1
        assert right_factor[3] == 2
        assert right_factor[4] == 1

        multi_degree = bch_series_generator.multi_degree
        assert multi_degree[0] == 1
        assert multi_degree[1] == 2
        assert multi_degree[2] == 4
        assert multi_degree[3] == 7
        assert multi_degree[4] == 8

    def test_goldberg_series(self):
        basis = LyndonBasis(2, Sort.Lexicographical)
        bch_series_generator = BCHSeriesGenerator(basis, 5)
        goldberg_coefficients = (
            bch_series_generator.generate_goldberg_coefficient_numerators()
        )
        assert len(goldberg_coefficients) == 14

    def test_bch_series(self):
        basis = LyndonBasis(2, Sort.Lexicographical)
        bch_series_generator = BCHSeriesGenerator(basis, 5)
        bch_series = bch_series_generator.generate_lie_series()
        assert (bch_series[0] - 1.0).__abs__() < 0.001
        assert (bch_series[1] - 1.0).__abs__() < 0.001
        assert (bch_series[2] - 1.0 / 2.0).__abs__() < 0.001
        assert (bch_series[3] - 1.0 / 12.0).__abs__() < 0.001
        assert (bch_series[4] - 1.0 / 12.0).__abs__() < 0.001
        assert (bch_series[5] - 0.0).__abs__() < 0.001
        assert (bch_series[6] - 1.0 / 24.0).__abs__() < 0.001
        assert (bch_series[7] - 0.0).__abs__() < 0.001
        assert (bch_series[8] - -1.0 / 720.0).__abs__() < 0.001
        assert (bch_series[9] - 1.0 / 180.0).__abs__() < 0.001
        assert (bch_series[10] - 1.0 / 360.0).__abs__() < 0.001
        assert (bch_series[11] - 1.0 / 180.0).__abs__() < 0.001
        assert (bch_series[12] - 1.0 / 120.0).__abs__() < 0.001
        assert (bch_series[13] - -1.0 / 720.0).__abs__() < 0.001
