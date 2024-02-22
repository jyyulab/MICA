#!/usr/bin/env python3

import unittest
import tempfile
import pandas as pd
import numpy.testing as np
from MICA.bin import transform


class TestTransform(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.infile = "./tests/inputs/test_dist.h5"
        self.folder = tempfile.TemporaryDirectory()
    
    def test_mi_mds_transform(self):
        transform.transform(self.infile, self.folder.name + "/test", "MDS", 20, "mi")
        pd_test = pd.HDFStore(self.folder.name + "/test_reduced.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/transform/test_reduced.h5")
        mds_test = pd_test['mds']
        mds_ans = pd_test['mds']
        for i in mds_test.index:
            arr_test = mds_test.loc[i]
            arr_ans = mds_ans.loc[i]
            np.assert_allclose(arr_test, arr_ans)


if __name__ == '__main__':
    unittest.main()

