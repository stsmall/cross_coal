# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 15:47:07 2021

@author: Scott T. Small

"""

import numpy as np

def calc_cross_coalescent(tree_seq, nodes_ls, ccN_events=10):
    """Calculate the cross coalescent of two lists of nodes.

    Parameters
    ----------
    tree_seq : Object
        object of type tskit tree seqeunce.
    nodes_ls : List
        List of node ids as integers, [[0, 1, 2],[4, 5, 6]]
    cc_N_events : int, optional
        the number of cross coalescent events to track. The default is 10.

    Returns
    -------

    ccN_rel : List
        the cross coalescent of the population from 1 ... N

    """
    ccN_ts = []
    pop_1 = nodes_ls[0]
    pop_2 = nodes_ls[1]
    iter1 = tree_seq.trees(tracked_samples=pop_1, sample_lists=True)
    iter2 = tree_seq.trees(tracked_samples=pop_2, sample_lists=True)
    p1_samples = set(pop_1)
    p2_samples = set(pop_2)
    for tree1, tree2 in zip(iter1, iter2):
        ccN_tree = []
        num_cc = 0
        used_nodes = set()
        for u in tree1.nodes(order='timeasc'):
            num_pop1 = tree1.num_tracked_samples(u)
            num_pop2 = tree2.num_tracked_samples(u)
            if num_cc < ccN_events:
                if num_pop1 > 0 and num_pop2 > 0:
                    proposed_cc1 = set(tree2.samples(u)) - p2_samples - used_nodes
                    proposed_cc2 = set(tree1.samples(u)) - p1_samples - used_nodes
                    if proposed_cc1 and proposed_cc2:
                        used_nodes |= proposed_cc1 | proposed_cc2
                        simul_cc_events = min(
                            [len(proposed_cc1), len(proposed_cc2)])
                        num_cc += simul_cc_events
                        ccN_tree.extend([u] * simul_cc_events)  # some nodes are the parent of >1 cc
        if num_cc < ccN_events:
            # padding so all trees return the same length
            ccN_tree.extend(np.repeat(np.nan, (ccN_events - num_cc)))
        ccN_ts.append(ccN_tree[:ccN_events])
    return ccN_ts
