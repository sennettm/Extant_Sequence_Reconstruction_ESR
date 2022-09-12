#!/usr/bin/env python3
#####################################################################################################
#Script Name    : print_leaf_taxa.py
#Description    : print  node labels for leaf on tree
#Author         : Michael A Sennett
#####################################################################################################

import dendropy
import copy
from dendropy.calculate import treecompare
from dendropy.calculate import treemeasure
import sys


######################################## INPUT VARIABLE #############################################
treefile = sys.argv[1] #tree that we want to know taxa of
#####################################################################################################

############################################ OUTPUT #################################################
#prints list of taxa in tree
#####################################################################################################

tree = dendropy.Tree.get(path=treefile,
                         schema="newick",
                         preserve_underscores=True)
treelen = tree.length()
leaves = tree.leaf_nodes()
leaf_num = len(leaves)
treeness = treemeasure.treeness(tree)


for node in tree.leaf_node_iter():
    brlen = node.edge_length
    print("%s" % node.taxon.label)

