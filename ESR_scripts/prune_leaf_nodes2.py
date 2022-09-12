#!/usr/bin/env python3
####################################################################################################
#Script Name    : prune_leaf_nodes2.py
#Description    : Import a newick tree and remove any extant sequences shorter than a defined length
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

import dendropy
import copy
from dendropy.calculate import treecompare
from dendropy.calculate import treemeasure
import sys

######################################## INPUT VARIABLE #############################################
treefile = sys.argv[1]      #tree to be pruned
minlen = float(sys.argv[2]) #branches to be pruned, less than or equal to

if len(sys.argv) >3:
    taxa_2_keep=sys.argv[3] #list of sequences to prevent from pruning
else:
    taxa_2_keep="NONE"
#####################################################################################################

############################################ OUTPUT #################################################
#treefile_pruned.tree, but with branches of certain length removed
#####################################################################################################

tree = dendropy.Tree.get(path=treefile,
                         schema="newick",
                         preserve_underscores=True)

treelen = tree.length()
leaves = tree.leaf_nodes()
leaf_num = len(leaves)
treeness = treemeasure.treeness(tree)

#check and see if taxa list was provided, otherwise produce empty set
if taxa_2_keep != "NONE":
    taxa_2_keep=open(taxa_2_keep,'r')
    taxa_list=[]
    for line in taxa_2_keep:
        taxa_list.append(str(line.strip("\n")))
    for item in taxa_list:
        item.replace("\'",'')
else:
    taxa_list=[]

print(taxa_list)

outfile = treefile + "_pruned.tree"

# remove node labels (like branch supports)
for b in tree.bipartition_edge_map:
    edge = tree.bipartition_edge_map[b]
    edge.head_node.label = "";

for node in tree.leaf_node_iter():
    brlen = node.edge_length
    if (brlen <= minlen):
        label=node.taxon.label
        print(label)
        if str(label) in taxa_list:
            print(label+' in taxa_2_keep')
        else:
            print( "%-16s %12.7f" % (node.taxon.label, brlen))
            parent = node.parent_node
            parent.remove_child(node)
            tree.update_bipartitions()

tree.write(path=outfile, schema="newick", unquoted_underscores=True, suppress_rooting=True)
