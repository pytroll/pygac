"""Test utilities."""

import numpy as np


class CalledWithArray(object):
    """Adapter for arrays in mock.assert_called_with()."""
    def __init__(self, array):
        self.array = array

    def __repr__(self):
        return repr(self.array)

    def __eq__(self, other):
        return np.all(self.array == other)
