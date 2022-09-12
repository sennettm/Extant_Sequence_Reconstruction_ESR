#!/bin/env bash
#####################################################################################################
#Script Name    : ms_infer_leaves.sh
#Description    : Creates a directory for each extant sequence and in that directory infers the 
#                 conditional distribution for that extant sequence 
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

####################################### INPUT VARIABLE ##############################################
ALGN=$1; #the modified allignment w/ pseudo seq
TREE=$2; #the modified tree
PARAM=$3; #the .iqtree file with model parameters
NDS=$4; #the number of nodes available to run reconstructions
#####################################################################################################

########################################### OUTPUT ##################################################
#extant folder to accomodate all the iqtree files
#.state file that has the extant recon post prob dist
#smp seq .txt file (not used in final analysis)
#####################################################################################################

echo "Running ms_infer_all_leaves.sh";

#function that runs script to reconstruct ancestors for a given alignment, tree, and set of model parameters

N=$(($NDS/4)); #the number of parallel IQ-Tree recons at once

echo "Number of runs in parallel: " "$N";

function do_asr()
{
    .././iqtree_recon_asr.sh ${ALGN} ${TAXA}_recon.tre ../${PARAM};
}

#make the extant taxa folder to perform reconstructions

for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre}; #getting taxa name
    echo ${TAXA};
    rm -rf ${TAXA};
    mkdir ${TAXA};
    cp ${TAXA}_recon.tre ${TAXA}/.; #cp tree to new direct 
    cp ${ALGN} ${TAXA}/.; #cp algnment to new direct
done

#perform ESR for each extant sequence in parallel

for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    cd ${TAXA}; 
    (  
        do_asr
        sleep 3;
    ) &
    if (( $(jobs -r -p | wc -l) > $N));
    then
        wait -n;
    fi
    cd ..;    
done
wait

#this pulls the smp recon from the iqtree state file using ancprobs
#note this sometimes this doesn't agree with the iqtree smp because its gaps are not quite right

for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    cd ${TAXA};
    NODE=$(.././find_pseudonode.py ${ALGN}.treefile ${TAXA} PSEUDO);
    NUM=$(echo ${NODE} | sed 's/[^0-9]*//g');
    echo '>'${TAXA}_${NUM} >> ${TAXA}_smp_recon.fasta;
    egrep -w "${NODE}" ${ALGN}.state | awk '{print $3}' >tmp; paste -s -d '\0' tmp >> ${TAXA}_smp_recon.fasta
    rm tmp;
    cd ..;
done
