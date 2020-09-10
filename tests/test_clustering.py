#!/usr/bin/env python3

import unittest
import tempfile
import cProfile
from MICA.bin import clustering as cl


class TestClustering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.infile = "./tests/inputs/test_reduced.h5"
        cls.folder = tempfile.TemporaryDirectory()
    
    def test_clustering(self):
        arr = []
        res = []
        self.maxDiff = None
        pr = cProfile.Profile()
        pr.enable()
        cl.clustering(self.infile, "mds", 5, 10, self.folder.name + "/test", "tsne", 0.1, 30, 20, 5, [20])
        pr.disable()
        pr.print_stats(sort='time')
        with open("{}/test_k5_tsne_ClusterMem.txt".format(self.folder.name), 'r') as fin:
            fin.readline()
            for line in fin:
                line = line.split('\t')[0:3]
                arr.append(line)

        with open('./tests/answerkey/clustering/test_k5_tsne_ClusterMem.txt', 'r') as ans:
            ans.readline()
            for line in ans:
                line = line.split('\t')[0:3]
                res.append(line)

        self.assertEqual(arr, res)


if __name__ == "__main__":
    unittest.main()

