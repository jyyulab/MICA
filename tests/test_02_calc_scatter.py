#!/usr/bin/ python3

import unittest
import tempfile
import pandas as pd
from MICA.bin import calc_scatter as cs


class TestCalcScatter(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.infile = './tests/inputs/test_mi_04.h5.tmp'
        self.folder = tempfile.TemporaryDirectory()
        pd_infile = pd.HDFStore(self.infile)
        params = pd_infile['params']
        params.loc['project_name'] = self.folder.name + "/test"
        pd_infile.close()
        params.to_hdf(self.infile, 'params')

    def test_mutual_information_calculation(self):
        cs.calc_scatter(self.infile, "mi")
        pd_test = pd.HDFStore(self.folder.name + "/test_mi_04.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/calc_scatter/test_mi_04.h5")
        self.assertTrue(pd_test['mi_04'].equals(pd_ans['mi_04']))    


if __name__ == '__main__':
    unittest.main()
