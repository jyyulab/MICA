#!/usr/bin/env python3

import os
import unittest
import tempfile
import pandas as pd
import numpy.testing as np
from MICA.bin import merge_and_norm as merge


class TestMerge(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.infile = "./tests/inputs/merge_and_norm/test_mi_0.h5", \
                "./tests/inputs/merge_and_norm/test_mi_1.h5", \
                "./tests/inputs/merge_and_norm/test_mi_2.h5"
        self.folder = tempfile.TemporaryDirectory()
        merge.merge_and_norm(self.infile, self.folder.name + "/test", "mi")

    def test_mi_dist_metric(self):
        pd_test = pd.HDFStore(self.folder.name + "/test_dist.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/merge_and_norm/test_dist.h5")
        self.assertTrue(pd_test['mi'].equals(pd_ans['mi']))

    def test_mi_dist_metric_normalized(self):
        pd_test = pd.HDFStore(self.folder.name + "/test_dist.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/merge_and_norm/test_dist.h5")
        norm_test = pd_test['norm_mi']
        norm_ans = pd_ans['norm_mi']
        for i in norm_test.index:
            arr_test = norm_test.loc[i]
            arr_ans = norm_ans.loc[i]
            np.assert_allclose(arr_test, arr_ans)

    def test_number_of_slices_preserved(self):     
        pd_test = pd.HDFStore(self.folder.name + "/test_mi_whole.h5")
        keys = pd_test.keys()
        self.assertEqual(keys, ['/mi_0', '/mi_1', '/mi_2'])

    def test_content_of_slices_preserved(self):
        pd_test = pd.HDFStore(self.folder.name + "/test_mi_whole.h5")
        pd_ans = pd.HDFStore("./tests/answerkey/merge_and_norm/test_mi_whole.h5")
        for i in range(0, 3):
            for j in pd_test['mi_{}'.format(i)].index:
                arr_test = pd_test['mi_{}'.format(i)].loc[j]
                arr_ans = pd_ans['mi_{}'.format(i)].loc[j]
                np.assert_allclose(arr_test, arr_ans)


if __name__ == '__main__':
    unittest.main()

