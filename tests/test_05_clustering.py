#!/usr/bin python3

import unittest
import tempfile
from MICA.bin import clustering as cl

class TestClustering(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.infile = "./tests/inputs/test_reduced.h5"
        self.folder = tempfile.TemporaryDirectory()

    def test_5_clustering(self):
        arr = []
        res = []
        cl.clustering(self.infile, "mds", 5, 10, self.folder.name + "/test", "tsne", 0.1, 30, 20, 10, [20])
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
