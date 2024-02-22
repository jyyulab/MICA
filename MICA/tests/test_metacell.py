#!/usr/bin/env python3

import unittest
import pandas as pd
import networkx as nx
from MICA.lib import metacell


class TestMetaCell(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.cls_res = {'cell': [0, 1, 2, 3, 4, 5], 'label': [1, 2, 2, 2, 2, 3]}
        self.cls_res_df = pd.DataFrame(data=self.cls_res)
        self.G = nx.Graph()
        self.G.add_nodes_from([0, 1, 2, 3, 4, 5])
        self.G.add_weighted_edges_from([(0, 5, 0.4), (1, 2, 0.3), (1, 3, 0.7), (1, 4, 0.5),
                                       (2, 4, 0.9), (4, 5, 0.2)], weight='weight')
        self.out_dir = './'

    def test_pick_cls_res(self):
        metacell.pick_cells(cls_res_df=self.cls_res_df, G=self.G, method='closeness', out_dir=self.out_dir)


if __name__ == "__main__":
    unittest.main()
