#!/bin/env bash
#this script runs ancprobs to get the smp sequences and appends them to a recon alignment file
#####################################################################################################
#Script Name    : ms_crunch_all_leaves.sh
#Description    : Runs ancprobs to collate SMP sequences into new alignment 
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

#####################################################################################################

######################################## INPUT VARIABLES ############################################
BL="$1"; #the length of branches trimmed
ALGN=$2; #the modified alignment
#TREE=$3; #the modified tree
#####################################################################################################

########################################### OUTPUT ##################################################
#residue stats as .out file
#recon alignment as .2m file
#####################################################################################################

echo "Running ms_crunch_all_leaves.sh";

#seqcon -E ${ALGN}; #makes algnment into individ fasta files
./ms_seqcon.py -E --MSA ${ALGN} --form1 fasta;

rm -f ${ALGN}_recon_${BL}.a2m; #rm previously made recon algnmt

#iterate through all trees, grab taxa & mrca node num
#pull gapped anc with ancprobs as out file and append previously created smp fasta to .a2m file 
for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    NODE=$(./find_pseudonode.py ${TAXA}/${ALGN}.treefile ${TAXA} PSEUDO);
    NUM=$(echo ${NODE} | sed 's/[^0-9]*//g');
    ./ms_ancprobs.py -N ${NUM} -F ${TAXA}/${ALGN}.state -G ${TAXA}.fst > ${TAXA}_${NUM}_ancprobs.out;
    #ancprobs -p iqtree -n ${NUM} -f ${TAXA}/${ALGN}.state -g ${TAXA}.fst > ${TAXA}_${NUM}_ancprobs.out;
    egrep '^SITE:' ${TAXA}_${NUM}_ancprobs.out | awk '{print $2 "\t" $3 "\t" $6}' > ${TAXA}_${NUM}_probs.txt;
    FILIN=$(cat ${TAXA}.fst| wc -l);
    LINE=$(expr ${FILIN} - 1);
    echo ">${TAXA}_recon" >> ${ALGN}_recon_${BL}.a2m;
    egrep -A ${LINE} "smp_gaps" ${TAXA}_${NUM}_ancprobs.out | tail -n${LINE} >> ${ALGN}_recon_${BL}.a2m;
done

