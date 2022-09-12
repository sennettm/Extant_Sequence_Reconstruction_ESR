#!/usr/bin/env python3
#####################################################################################################
#Script Name    : extant_entropy.py
#Description    : calculate the negative of entropy (i.e. expect LnP) for an ancestral distribution
#Author         : Michael A Sennett
#####################################################################################################

import numpy as np
import sys

######################################## INPUT VARIABLE #############################################
state_file=sys.argv[1] #iqtree .state file containing ancestral distributions
node=sys.argv[2]       #node corresponding to sequence of interest
gap_file=sys.argv[3]   #reference sequence that we use to identify sites we shouldn't include
#####################################################################################################

############################################ OUTPUT #################################################
#eLnP of an ancestral distribution
#####################################################################################################

#open readable state file
tmp_state_file = open(state_file, "r")
statefile_lines = tmp_state_file.readlines()
tmp_state_file.close()

#open readable gap file
tmp_gap_file = open(gap_file, "r")
gapfile_lines = tmp_gap_file.readlines()
tmp_gap_file.close()

seq=[]
for line in gapfile_lines:
    if line.startswith(">"):
        pass
    else:
        line=line.replace('\n','')
        seq.append(line)
seq="".join(seq)

pp_vector=[]
for line in statefile_lines:
    if line.startswith('Node'+str(node)+'\t'):
        line=line.split()
        pp_vector.append(line)
    else:
        pass

#reformat data and remove first 3 columns
pp_vector=np.reshape(pp_vector,(-1,len(pp_vector[0])))
pp_vector=np.delete(pp_vector,[0,1,2],axis=1)

#calculate the entropy for each row and total the entropy of the sequence
entropy=[]
for i in range(len(pp_vector)):
    if seq[i] != '-':
        site_ent=[]
        for j in range(len(pp_vector[0])):
            if pp_vector[i][j] == '0.00000':
                site_ent.append(0)
            else:
                log=np.log(float(pp_vector[i][j]))
                prob=float(pp_vector[i][j])
                site_ent.append(prob*log)
        entropy.append(sum(site_ent))
    else:
        pass
tot_ent=sum(entropy)
print(tot_ent)




