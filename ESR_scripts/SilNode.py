#!/usr/bin/env python3
#####################################################################################################
#Script Name    : SilNode.py
#Description    : Import a newick tree with node labels and export a newick tree without node labelsl
#Author         : Michael A Sennett
#####################################################################################################

import dendropy as dp
import sys

######################################## INPUT VARIABLE #############################################
treefile=sys.argv[1] #input tree
#####################################################################################################

############################################ OUTPUT #################################################
#treefile, but without internal node labels
#####################################################################################################

treefile=sys.argv[1]
tree=dp.Tree.get(path=treefile, schema='newick', preserve_underscores=True)
outfile=treefile
tree.write(path=outfile, schema='newick', suppress_internal_node_labels=True, unquoted_underscores=True)




