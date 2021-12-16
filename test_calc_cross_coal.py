# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 15:47:07 2021

@author: Scott T. Small

"""

import msprime
print(f"msprime version is {msprime.__version__}; tested with 1.0.1")

from .calc_cross_coal import calc_cross_coalescent

def test_ccN():
    p_nodes1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    p_nodes2 = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
    p_nodes_cc = [p_nodes1, p_nodes2]
    tss = msprime.sim_ancestry(16, random_seed=234)
    ccN = calc_cross_coalescent(tss, p_nodes_cc, 10)
    assert(ccN[0] == [34, 35, 39, 41, 44, 46, 53, 53, 54, 57])
