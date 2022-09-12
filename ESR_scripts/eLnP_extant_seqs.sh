#!/bin/env bash
#####################################################################################################
#Script Name    : eLnP_extant_seqs.sh
#Description    : Find the extant distribution and calculate the expected Ln(probability) 
#                 of the sequence (eLnP), eLnp=-H   
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

#####################################################################################################

######################################## INPUT VARIABLE #############################################
BL=$1; #pruned branch length
ALGN=$2; #modified alignment
#####################################################################################################

############################################ OUTPUT #################################################
#eLnP_${BL}.txt which is a list of all eLnP for all extant sequences
#####################################################################################################

echo "Running eLnP_extant_seqs.sh";

#remove previously made file
rm eLnP_${BL}.txt

#go through and calculate the eLnP of a sequence and append to 
for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    cd ${TAXA};
    NODE=$(.././find_pseudonode.py ${ALGN}.treefile ${TAXA} PSEUDO);
    NUM=$(echo ${NODE} | sed 's/[^0-9]*//g');
    ENT=$(.././extant_entropy.py ${ALGN}.state ${NUM} ../${TAXA}.fst);
    cd ..;
    echo "${TAXA} ${ENT}" >> eLnP_${BL}.txt;
done
