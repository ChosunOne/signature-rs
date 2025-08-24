from lie_py import BCHSeriesGenerator
from lyndon_py import LyndonBasis, Sort


class TestBCHSeriesGeneratorGroup:
    def test_constructor(self):
        basis = LyndonBasis(2, Sort.Lexicographical)
        bch_series_generator = BCHSeriesGenerator(basis, 3)
