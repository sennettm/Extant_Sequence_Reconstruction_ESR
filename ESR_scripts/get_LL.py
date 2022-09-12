#!/usr/bin/env python3
#####################################################################################################
#Script Name    : get_LL.py
#Description    : calculates the log probability of some sequence for a given ancestral distribution
#Author         : Michael A Sennett
#####################################################################################################

import sys

######################################## INPUT VARIABLE #############################################
state_file=sys.argv[1]  #iqtree .state file containing ancestral distributions
gap_file=sys.argv[2]    #file so we know which sites not to infer a residue for
node=sys.argv[3]        #the node that contains distributions we need
#####################################################################################################

############################################ OUTPUT #################################################
#prints the LnP of a sequence
#####################################################################################################


#import statefile and store as a list
statefile=open(state_file, 'r')
state_lines=statefile.readlines()
statefile.close()

#import gapping file
extant_fst=open(gap_file, 'r')
seq_lines=extant_fst.readlines()
extant_fst.close()

#grab the PP dist text for node of interest
import numpy as np
node_lines=[]
for line in state_lines:
    if line.startswith("Node"+str(node)+'\t'):
        line=line.replace('\n','')
        node_lines.append(line.split('\t'))
    else:
        pass

node_vector=[]
node_vector=np.reshape(node_lines,(-1,len(node_lines[0])))
PP_vector=np.delete(node_vector, [0,1,2], 1)

#dictionary for amino acids in alphabetical order
AA_dict={"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, "H":8, "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, "T":16, "W":17, "Y":18, "V":19}

#create a single string of the fst sequence
seq=[]
for line in seq_lines:
    if line.startswith(">"):
        pass
    else:
        line=line.replace('\n','')
        seq.append(line)
seq="".join(seq)

#iterate through each element of sequence to grab probability of amino acid
ln_seq=[]
for element in range(len(seq)):
    if seq[element] == '-':
        pass
    elif seq[element] == 'X':
        ln_prob=0
        ln_seq.append(ln_prob)
    else:
        index=AA_dict.get(seq[element])
        prob=PP_vector[element][index]
        ln_prob=np.log(float(prob))
        ln_seq.append(ln_prob)

print(sum(ln_seq))
