from unittest import TestCase

import pandas as pd
import numpy as np

from tracegroomer import utils

class Test(TestCase):
    def test_compute_sums_isotopol_props(self):
        df = pd.DataFrame({"metabolite": ['a', 'a', 'b', 'b'],
              "col1": [0.3,0.7, 0.2,0.8],
              "col2": [0.7,0.3, 0.4,0.6]})
        result = utils.compute_sums_isotopol_props(df)
        self.assertEqual(result.loc['a','col1'], 1)
        self.assertEqual(result.loc['a', 'col2'], 1)
        self.assertEqual(result.loc['b', 'col1'], 1)
        self.assertEqual(result.loc['b', 'col2'], 1)