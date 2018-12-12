import os
import sys

import numpy as np
import pytest

sys.path.append(os.getenv('PYTAYLOR_MODULE_PATH'))

import pytaylor  # isort:skip

xs = np.linspace(-5.0, 5.0)
ys = np.linspace(-5.0, 5.0)


@pytest.mark.parametrize('x, y, f, Df_x, Df_y',
                         list(zip(xs, ys, xs + ys, np.full_like(xs, 1.0) + ys, xs + np.full_like(ys, 1.0))))
def test_bivariate_add(x, y, f, Df_x, Df_y):
    seed_x = pytaylor._pytaylor.taylorD_2_1(x, 0)
    seed_y = pytaylor._pytaylor.taylorD_2_1(y, 1)
    res = seed_x + seed_y
    assert res[0] == f
    assert res[1] == Df_x
    assert res[2] == Df_y


@pytest.mark.parametrize('x, f, Df', list(zip(xs, xs * xs, 2.0 * xs)))
def test_univariate_mul(x, f, Df):
    seed = pytaylor._pytaylor.taylorD_1_1(x, 0)
    res = seed * seed
    assert res[0] == f
    assert res[1] == Df
