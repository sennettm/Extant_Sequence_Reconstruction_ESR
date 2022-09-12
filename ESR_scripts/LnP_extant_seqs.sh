#!/bin/env bash
#####################################################################################################
#Script Name    : LnP_extant_seqs.sh
#Description    : Find the extant distribution and calculate the Ln(probability) 
#                 of the true sequence  
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

#####################################################################################################

######################################## INPUT VARIABLE #############################################
BL=$1; #pruned branch length
ALGN=$2; #modified alignment
#####################################################################################################

############################################ OUTPUT #################################################
#tru_LnP_${BL}.txt which is a list of all eLnP for all extant sequences
#####################################################################################################

echo "Running LnP_extant_seqs.sh";

rm tru_LnP_${BL}.txt;

#collate the LnP of the true sequence
for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    SEQ=${TAXA}.fst;
    cd ${TAXA};
    cp ../${TAXA}.fst .;
    NODE=$(.././find_pseudonode.py ${ALGN}.treefile ${TAXA} PSEUDO);
    NUM=$(echo ${NODE} | sed 's/[^0-9]*//g');
    LL=$(.././get_LL.py ${ALGN}.state ${SEQ} ${NUM});
    rm ${TAXA}.fst
    cd ..;
    echo "${TAXA} ${LL}" >> tru_LnP_${BL}.txt;
done
