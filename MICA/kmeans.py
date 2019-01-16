#!/usr/bin/env python3

import sys

import utils

utils.kmeans(
    sys.argv[1],
    sys.argv[2],
    int(sys.argv[3]),
    int(sys.argv[4]),
    sys.argv[5],
    sys.argv[6],
)
