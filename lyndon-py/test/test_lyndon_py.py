from lyndon_py import LyndonBasis, Sort


class TestLyndonPyGroup:
    def test_constructor(self):
        basis = LyndonBasis(3, Sort.Lexicographical)

    def test_generate_a_lyndon_basis(self):
        basis = LyndonBasis(2, Sort.Topological).generate_basis(5)
        expected_basis = []
        assert basis == expected_basis
