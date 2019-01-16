#!/usr/bin/env bash

./MICA/mica.py Clust test_no_install_local ./test_data/inputs/test_local.whole.h5 \
./test_data/inputs/test_local_mi.h5 ./test_data/outputs/test_no_install_local/ test_no_install_local --k 3 4 5 6 \
--perplexity 30 --retransformation False
