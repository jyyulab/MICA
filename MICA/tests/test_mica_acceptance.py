#!/usr/bin/env python3

import unittest
import tempfile
import shlex
import subprocess
import numpy as np


class TestMICA(unittest.TestCase):
    def test_acceptance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = 'mica local -i ./tests/inputs/PBMC_Demo_MICA_input_mini.txt -p test -k 4 -o {} -sn 200'.format(tmpdir)
            exe = shlex.split(cmd)
            subprocess.check_call(exe)
            
            arr = []
            res = []
            arr_nums = []
            res_nums = []

            with open('{}/test_k4_tsne_ClusterMem.txt'.format(tmpdir), 'r') as fin:
                fin.readline()
                for line in fin:
                    line = line.split('\t')[0:3]
                    arr.append(line[0])
                    arr_nums.append(float(line[1]))
                    arr_nums.append(float(line[2]))

            with open('./tests/answerkey/cwl_local_k4_tsne_ClusterMem.txt', 'r') as ans:
                ans.readline()
                for line in ans:
                    line = line.split('\t')[0:3]
                    res.append(line[0])
                    res_nums.append(float(line[1]))
                    res_nums.append(float(line[2]))

            np_arr = np.array(arr_nums)
            np_res = np.array(res_nums)

            self.assertEqual(arr, res)
            np.testing.assert_almost_equal(np_arr, np_res, decimal=1, err_msg='Current output in MICA diverges from result by more than a decimal point')

if __name__ == '__main__':
    unittest.main()
