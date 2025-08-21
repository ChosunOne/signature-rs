import pytest
from lyndon_py import LyndonBasis, LyndonWord, Sort


class TestLyndonPyGroup:
    def test_constructor(self):
        basis = LyndonBasis(3, Sort.Lexicographical)
        assert basis.alphabet_size == 3
        assert basis.sort == Sort.Lexicographical

    def test_str(self):
        word = LyndonWord([0, 1, 2, 3, 4])
        assert word.__str__() == "01234"

    def test_generate_a_lyndon_basis(self):
        basis = LyndonBasis(2, Sort.Topological).generate_basis(5)
        expected_basis = [
            LyndonWord([0]),
            LyndonWord([1]),
            LyndonWord([0, 1]),
            LyndonWord([0, 1, 1]),
            LyndonWord([0, 0, 1]),
            LyndonWord([0, 1, 1, 1]),
            LyndonWord([0, 0, 1, 1]),
            LyndonWord([0, 0, 0, 1]),
            LyndonWord([0, 1, 1, 1, 1]),
            LyndonWord([0, 0, 1, 0, 1]),
            LyndonWord([0, 1, 0, 1, 1]),
            LyndonWord([0, 0, 1, 1, 1]),
            LyndonWord([0, 0, 0, 1, 1]),
            LyndonWord([0, 0, 0, 0, 1]),
        ]
        assert basis == expected_basis

    @pytest.mark.parametrize(
        "word,expected_partition",
        [
            ([0], [1]),
            ([1], [1]),
            ([0, 1], [1, 1]),
            ([0, 1, 1], [2, 1]),
            ([0, 0, 1], [2, 1]),
            ([0, 1, 1, 1], [3, 1]),
            ([0, 1, 0, 1, 1], [2, 1, 1, 1]),
            ([0, 0, 1, 0, 1], [2, 1, 1, 1]),
        ],
    )
    def test_generate_goldberg_partitions(
        self, word: list[int], expected_partition: list[int]
    ):
        lyndon_word = LyndonWord(word)
        partition = lyndon_word.goldberg()
        assert partition == expected_partition

    def test_words_per_degree(self):
        basis = LyndonBasis(2, Sort.Lexicographical)
        num_words_per_degree = basis.number_of_words_per_degree(5)
        expected_num_words_per_degree = [2, 1, 2, 3, 6]
        assert len(num_words_per_degree) == len(expected_num_words_per_degree)
        for num, expected_num in zip(
            num_words_per_degree, expected_num_words_per_degree
        ):
            assert num == expected_num

    def test_letters(self):
        word = LyndonWord([0, 1, 2, 3])
        assert word.letters == [0, 1, 2, 3]

    def test_frozen_letters(self):
        word = LyndonWord([0, 1, 2, 3])
        with pytest.raises(AttributeError) as e_info:
            word.letters = [3, 2, 1, 0]
        assert (
            e_info.value.args[0]
            == "attribute 'letters' of 'builtins.LyndonWord' objects is not writable"
        )
        assert word.letters == [0, 1, 2, 3]

    def test_mul(self):
        word1 = LyndonWord([0, 1])
        word2 = LyndonWord([2, 3])
        word = word1 * word2
        assert word == LyndonWord([0, 1, 2, 3])

    def test_mul_error(self):
        word1 = LyndonWord([2, 3])
        word2 = LyndonWord([0, 1])
        with pytest.raises(ValueError) as e_info:
            _ = word1 * word2
        assert e_info.value.args[0] == "Invalid lyndon word"
