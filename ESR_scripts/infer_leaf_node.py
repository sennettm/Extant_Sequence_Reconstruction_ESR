#!/usr/bin/env python3
#####################################################################################################
#Script Name    : infer_leaf_node.py
#Description    : Import a newick tree with N taxa, and export N trees where each external node is 
#                 internalized so we can infer an ancestral distribution 
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

import dendropy
import copy
from dendropy.calculate import treecompare
from dendropy.calculate import treemeasure
import sys

######################################## INPUT VARIABLE #############################################
treefile=sys.argv[1] #input tree
#####################################################################################################

############################################ OUTPUT #################################################
#${TAXA}_recon.tre which contains an internalized ${TAXA} node
#####################################################################################################


tree = dendropy.Tree.get(path=treefile,
                         schema="newick",
                         preserve_underscores=True)

tn = tree.taxon_namespace
tree.update_bipartitions(suppress_unifurcations=False)

treelen = tree.length()
leaves = tree.leaf_nodes()
leaf_num = len(leaves)
treeness = treemeasure.treeness(tree)

outfile = treefile + "_dendropy.tree"
print("%s\n" % outfile)

# remove node labels (like branch supports, which paml won't read)
for b in tree.bipartition_edge_map:
    edge = tree.bipartition_edge_map[b]
    edge.head_node.label = "";


#creates a set of modified trees, whereby each leaf is internalized one at a time
brlen = 1000.0

for node in tree.leaf_node_iter():
    ch1 = dendropy.Node(edge_length=brlen)
    ch1.taxon = tn.new_taxon("PSEUDO")
    node.add_child(ch1)
    ch2 = dendropy.Node(edge_length=brlen)
    ch2.taxon = tn.new_taxon(node.taxon.label)
    node.add_child(ch2)
    label = node.taxon.label
    node.taxon.label = ""
    outfile = label + "_recon.tre"
    print(outfile)
    print(tree.as_string("newick"))
    tree.write(path=outfile, schema="newick", unquoted_underscores=True, suppress_rooting=True)
    node.remove_child(ch1, suppress_unifurcations=False)
    node.remove_child(ch2, suppress_unifurcations=False)
    node.taxon.label = label
