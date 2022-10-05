#!/usr/bin/env bash
#####################################################################################################
#Script Name    : LnP_SMP_seqs.sh
#Description    : Pull the Ln(probability) of the SMP sequence  
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

#####################################################################################################

######################################## INPUT VARIABLE #############################################
BL=$1; #pruned branch lengths
#####################################################################################################

############################################ OUTPUT #################################################
#SMP_LnP_${BL}.txt which is a list of all SMP sequence LnPs for all extant sequences
#####################################################################################################

echo "Running LnP_SMP_seqs.sh";

rm SMP_LnP_${BL}.txt;

for TAXA in *_recon.tre; do
	SEQ=${TAXA%_recon.tre};
	LL=$(egrep -m 1 "log prob" ${SEQ}_*.out | awk '{print $3}')
	echo "${SEQ}" "${LL}" >> SMP_LnP_${BL}.txt;
done
