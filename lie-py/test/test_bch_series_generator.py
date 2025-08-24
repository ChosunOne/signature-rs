from lie_py import BCHSeriesGenerator
from lyndon_py import LyndonBasis, LyndonWord, Sort


class TestBCHSeriesGeneratorGroup:
    def test_constructor(self):
        basis = LyndonBasis(2, Sort.Lexicographical)
        bch_series_generator = BCHSeriesGenerator(basis, 3)
        basis = bch_series_generator.basis
        assert basis[0] == LyndonWord([0])
