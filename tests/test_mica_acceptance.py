#!/usr/bin/env python3

import unittest
import tempfile
import shlex
import subprocess
# import filecmp


class TestMICA(unittest.TestCase):
    def test_acceptance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = 'mica local -i ./tests/inputs/PBMC_Demo_MICA_input_mini.txt -p test -k 4 -o {} -sn 200'.format(tmpdir)
            exe = shlex.split(cmd)
            subprocess.check_call(exe)

            arr = []
            res = []
            with open('{}/test_k4_tsne_ClusterMem.txt'.format(tmpdir), 'r') as fin:
                fin.readline()
                for line in fin:
                    line = line.split('\t')[0:3]
                    arr.append(line)

            with open('./tests/answerkey/cwl_local_k4_tsne_ClusterMem.txt', 'r') as ans:
                ans.readline()
                for line in ans:
                    line = line.split('\t')[0:3]
                    res.append(line)

            self.assertEqual(arr, res)
            # self.assertTrue(filecmp.cmp('./tests/answerkey/cwl_local_k4_tsne_ClusterMem.txt',
            # '{}/test_k4_tsne_ClusterMem.txt'.format(tmpdir)))


if __name__ == '__main__':
    unittest.main()
