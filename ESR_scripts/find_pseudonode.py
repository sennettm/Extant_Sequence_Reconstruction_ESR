#!/usr/bin/env python3
#####################################################################################################
#Script Name    : find_pseudonode.py
#Description    : finds the previously internalized extant node
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

import dendropy
import sys

######################################## INPUT VARIABLE #############################################
treefile = sys.argv[1] #input _recon.tre
tx1 = str(sys.argv[2]) #PSEUDO or TAXA
tx2 = str(sys.argv[3]) #PSEUDO or TAXA
#####################################################################################################

############################################ OUTPUT #################################################
#Ancestral node corresponding to extant sequence we want to reconstruct
#####################################################################################################

tree = dendropy.Tree.get(path=treefile,
                         schema="newick",
                         preserve_underscores=True)

tn = tree.taxon_namespace

mrca = tree.mrca(taxon_labels=[tx1, tx2])

print(mrca.label)
