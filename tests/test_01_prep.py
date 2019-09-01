#!/usr/bin/env python3

import unittest
import tempfile
import pandas as pd
from MICA.bin import prep


class TestPrep(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.infile = "./tests/inputs/PBMC_Demo_MICA_input_mini.txt"
        self.folder = tempfile.TemporaryDirectory()
        prep.prep(self.infile, self.folder.name + "/test", 250)

    def test_whole_HDF5_file(self):
        pd_test = pd.HDFStore(self.folder.name + "/test.whole.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/prep/test.whole.h5")
        self.assertTrue(pd_test['slice_0'].equals(pd_ans['slice_0']))

    def test_sliced_HDF5_file(self):
        pd_test = pd.HDFStore(self.folder.name + "/test.sliced.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/prep/test.sliced.h5")
        self.assertTrue(pd_test['slice_0'].equals(pd_ans['slice_0']) and pd_test['slice_1'].equals(pd_ans['slice_1']))


if __name__ == '__main__':
    unittest.main()
