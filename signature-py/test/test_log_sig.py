import pytest
from jax import Array, jit
from jax.numpy import array
from numpy import float32
from numpy.typing import NDArray
from signature_py import LogSignatureBuilder


class TestLogSignatureGroup:
    def test_builder(self):
        builder = LogSignatureBuilder(3, 4)
        assert builder.max_degree == 3
        assert builder.num_dimensions == 4
        log_signature = builder.build()
        series_coefficients = log_signature.series.coefficients
        for c in series_coefficients:
            assert c == 0.0

        bch_series = log_signature.bch_series
        assert len(bch_series.coefficients) == 5

    @pytest.mark.parametrize(
        "num_dimensions,max_degree,path,expected_coefficients",
        [
            (
                3,
                3,
                array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]).__array__(),
                [
                    1.000000,
                    2.000000,
                    3.000000,
                    0.0,
                    0.0,
                    0.0,
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
                array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 2.0, 3.0],
                        [6.0, 5.0, 4.0],
                    ]
                ).__array__(),
                [
                    6.0,
                    5.0,
                    4.0,
                    -3.5,
                    -7.0,
                    -3.5,
                    2.333333,
                    4.666667,
                    -0.583333,
                    3.5,
                    0.0,
                    2.333333,
                    0.583333,
                    1.166667,
                ],
            ),
            (
                3,
                3,
                array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 2.0, 3.0],
                        [6.0, 5.0, 4.0],
                        [7.0, 8.0, 9.0],
                        [12.0, 11.0, 10.0],
                    ]
                ).__array__(),
                [
                    12.000000,
                    11.000000,
                    10.000000,
                    -6.500000,
                    -13.000000,
                    -6.500000,
                    -1.166667,
                    -2.333333,
                    4.416667,
                    6.500000,
                    16.500000,
                    15.333333,
                    -4.416667,
                    7.666667,
                ],
            ),
            (
                3,
                3,
                array(
                    [
                        [0.0, 0.0, 0.0],
                        [-0.077, 0.042, -0.067],
                        [-0.154, 0.675, 0.006],
                        [0.916, 1.177, -0.139],
                        [1.095, 0.823, -0.261],
                    ]
                ).__array__(),
                [
                    1.095000,
                    0.823000,
                    -0.261000,
                    -0.690006,
                    -0.040871,
                    -0.124105,
                    0.098690,
                    -0.004304,
                    0.146613,
                    0.024960,
                    0.044713,
                    -0.000215,
                    -0.038903,
                    0.001568,
                ],
            ),
            (
                3,
                4,
                array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 2.0, 3.0],
                        [6.0, 5.0, 4.0],
                    ]
                ).__array__(),
                [
                    6.000000,
                    5.000000,
                    4.000000,
                    -3.500000,
                    -7.000000,
                    -3.500000,
                    2.333333,
                    4.666667,
                    -0.583333,
                    3.500000,
                    0.000000,
                    2.333333,
                    0.583333,
                    1.166667,
                    1.458333,
                    2.916667,
                    -3.791667,
                    -3.208333,
                    -12.250000,
                    -9.333333,
                    -1.458333,
                    1.750000,
                    3.500000,
                    2.916667,
                    -3.791667,
                    6.708333,
                    2.625000,
                    7.291667,
                    1.750000,
                    1.750000,
                    -3.208333,
                    0.875000,
                ],
            ),
            (
                3,
                4,
                array(
                    [
                        [0.000, 0.000, 0.000],
                        [1.000, 2.000, 3.000],
                        [6.000, 5.000, 4.000],
                        [7.000, 8.000, 9.000],
                        [12.000, 11.000, 10.000],
                    ]
                ).__array__(),
                [
                    12.000000,
                    11.000000,
                    10.000000,
                    -6.500000,
                    -13.000000,
                    -6.500000,
                    -1.166667,
                    -2.333333,
                    4.416667,
                    6.500000,
                    16.500000,
                    15.333333,
                    -4.416667,
                    7.666667,
                    -0.041667,
                    -0.083333,
                    1.208333,
                    2.291667,
                    4.750000,
                    4.666667,
                    0.041667,
                    -2.250000,
                    -4.500000,
                    -11.083333,
                    -4.291667,
                    -12.291667,
                    -19.875000,
                    -22.208333,
                    -13.250000,
                    -2.250000,
                    7.791667,
                    -6.625000,
                ],
            ),
            (
                4,
                4,
                array(
                    [
                        [0.000, 0.000, 0.000, 0.000],
                        [1.000, 2.000, 3.000, 4.000],
                        [6.000, 5.000, 4.000, 3.000],
                        [7.000, 8.000, 9.000, 8.000],
                        [12.000, 11.000, 10.000, 9.000],
                    ]
                ).__array__(),
                [
                    12.000000,
                    11.000000,
                    10.000000,
                    9.000000,
                    -6.500000,
                    -13.000000,
                    -13.500000,
                    -6.500000,
                    -7.000000,
                    -0.500000,
                    -1.166667,
                    -2.333333,
                    10.500000,
                    4.416667,
                    6.500000,
                    18.583333,
                    16.500000,
                    15.333333,
                    26.666667,
                    5.166667,
                    20.833333,
                    7.750000,
                    -4.416667,
                    6.166667,
                    7.666667,
                    17.500000,
                    6.250000,
                    0.833333,
                    8.333333,
                    -6.083333,
                    -0.041667,
                    -0.083333,
                    -14.625000,
                    1.208333,
                    2.291667,
                    -19.125000,
                    4.750000,
                    4.666667,
                    -23.625000,
                    21.083333,
                    12.916667,
                    2.375000,
                    0.041667,
                    8.083333,
                    -2.250000,
                    -4.500000,
                    -18.750000,
                    -11.083333,
                    -4.291667,
                    -24.500000,
                    4.583333,
                    8.416667,
                    5.750000,
                    1.541667,
                    -12.291667,
                    -19.875000,
                    -9.958333,
                    -22.208333,
                    -13.250000,
                    -25.375000,
                    1.083333,
                    -21.458333,
                    9.125000,
                    -14.833333,
                    -16.666667,
                    -9.500000,
                    -38.250000,
                    -31.791667,
                    -19.000000,
                    -12.416667,
                    -22.458333,
                    -4.500000,
                    -2.250000,
                    -11.500000,
                    7.791667,
                    -8.166667,
                    22.666667,
                    10.166667,
                    4.250000,
                    -6.625000,
                    -15.250000,
                    -8.166667,
                    7.916667,
                    -18.458333,
                    -9.500000,
                    -14.583333,
                    -3.333333,
                    -5.125000,
                    6.708333,
                    -2.166667,
                ],
            ),
        ],
    )
    def test_build_from_path(
        self,
        num_dimensions: int,
        max_degree: int,
        path: NDArray[float32],
        expected_coefficients: list[float],
    ):
        builder = LogSignatureBuilder(max_degree, num_dimensions)
        log_signature = builder.build_from_path(path)
        coefficients = log_signature.series.coefficients
        assert len(coefficients) == len(expected_coefficients)
        for c, e_c in zip(coefficients, expected_coefficients):
            assert (c - e_c).__abs__() < 0.001, f"{c} != {e_c}"

    def test_concatenate(self):
        builder = LogSignatureBuilder(5, 5)
        log_sig_1 = builder.build_from_path(
            array([[1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0]]).__array__()
        )
        log_sig_2 = builder.build_from_path(
            array(
                [
                    [6.0, 7.0, 8.0, 9.0, 10.0],
                    [11.0, 12.0, 13.0, 14.0, 15.0],
                    [16.0, 17.0, 18.0, 19.0, 20.0],
                ]
            ).__array__()
        )
        log_sig = builder.build_from_path(
            array(
                [
                    [1.0, 2.0, 3.0, 4.0, 5.0],
                    [6.0, 7.0, 8.0, 9.0, 10.0],
                    [11.0, 12.0, 13.0, 14.0, 15.0],
                    [16.0, 17.0, 18.0, 19.0, 20.0],
                ]
            ).__array__()
        )

        concatenated_log_sig = log_sig_1.concatenate(log_sig_2)
        concat_coefficients = concatenated_log_sig.series.coefficients
        path_coefficients = log_sig.series.coefficients

        assert len(concat_coefficients) == len(path_coefficients)
        for cc, pc in zip(concat_coefficients, path_coefficients):
            assert (cc - pc).__abs__() < 0.001, f"{cc} != {pc}"

    def test_multiple_log_sigs(self):
        path = array([[1.0, 2.0], [3.0, 4.0]]).__array__()

        for degree in [1, 2, 3, 4]:
            builder = LogSignatureBuilder(degree, 2)
            _ = builder.build_from_path(path)
