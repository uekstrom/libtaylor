import os
import sys

sys.path.append(os.getenv('PYTAYLOR_MODULE_PATH'))

import pytaylor # isort:skip

def test_taylor_import():
    assert 'pytaylor' in sys.modules
