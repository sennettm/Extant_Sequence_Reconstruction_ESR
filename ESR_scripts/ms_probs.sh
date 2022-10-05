#!/usr/bin/env bash
#####################################################################################################
#Script Name    : ms_probs.sh
#Description    : Generates the gapped alignment of SMP reconstructions and calculates the actual fraction
#                 correct of the SMP sequence
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

######################################## INPUT VARIABLE #############################################
BL=$1; #branch length pruning
ALGN=$2; #modified alignment
#####################################################################################################

############################################ OUTPUT #################################################
#realprobs.txt file that has list of smp fraction correct
#aveprobs.txt file that has list of smp avg post prob
#recon and true alignments
#####################################################################################################

echo "Running ms_probs.sh";

#in case you run this more than once
rm realprobs_${BL}.txt;
rm aveprobs_${BL}.txt;
rm ${ALGN}_extant_${BL}.a2m;

#grabbing all the avg post prob from the ancporbs outfile
for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    cat ${TAXA}.fst >> ${ALGN}_extant_${BL}.a2m;
    PROB=$(ls *probs.out | egrep "\b${TAXA}_");
    egrep -m 1 'ave prob' ${PROB} | awk '{print $3}' >> aveprobs_${BL}.txt;
done

#make sure the alingments are in fasta format
#seqcon -f fasta ${ALGN}_extant_${BL}.a2m
#seqcon -f fasta ${ALGN}_recon_${BL}.a2m

EALGN=${ALGN}_extant_${BL}.a2m;
RALGN=${ALGN}_recon_${BL}.a2m;

#applies gaps from the extant alignment to the recon alignment
RALGNG=$(./ms_seqcon.py -G --MSA ${EALGN} ${RALGN} --form1 fasta);

#gets fraction correct between extant and recon alignment
./ms_seqcon.py -C --MSA ${EALGN} ${RALGNG} --form1 fasta;
mv fraction_correct.txt realprobs_${BL}.txt;
