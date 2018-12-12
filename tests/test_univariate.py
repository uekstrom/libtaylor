import os
import sys

import numpy as np
import pytest

sys.path.append(os.getenv('PYTAYLOR_MODULE_PATH'))

import pytaylor  # isort:skip

xs = np.linspace(-5.0, 5.0)


@pytest.mark.parametrize('x, f, Df',
                         list(zip(xs, xs + xs, np.full_like(xs, 2.0))))
def test_univariate_add(x, f, Df):
    seed = pytaylor._pytaylor.taylorD_1_1(x, 0)
    res = seed + seed
    assert res[0] == f
    assert res[1] == Df


@pytest.mark.parametrize('x, f, Df', list(zip(xs, xs * xs, 2.0 * xs)))
def test_univariate_mul(x, f, Df):
    seed = pytaylor._pytaylor.taylorD_1_1(x, 0)
    res = seed * seed
    assert res[0] == f
    assert res[1] == Df
