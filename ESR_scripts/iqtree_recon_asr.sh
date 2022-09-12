#!/bin/env bash
#####################################################################################################
#Script Name    : iqtree_recon_asr.sh
#Description    : performs ancestral sequence reconstruction using iqtree
#Author         : Michael A Sennett
#####################################################################################################

###################################### INPUT VARIABLES ##############################################
ALN=$1; #the alignment that has the pseudo sequence appended
TRE=$2; #the _recon.tre file that has the additional daughter nodes internalizing an extant seq
PRM=$3; #the .iqtree file with the equil freq and alpha 
#####################################################################################################

#########################################  OUTPUT  ##################################################
#ALN.state, .iqtree, .tree, .log files from IQ-Tree in ${TAXA} subfolder
#####################################################################################################

MODEL=$(egrep "Model of substitution:" ${PRM} | awk '{print $4}'); #model

IFS='+' read -a ARR <<< "$MODEL"; #creates an array containing  each model param

SM=${ARR[0]};

REC_MDL=$(echo $SM);

for p in ${ARR[@]}; do
	if [ $p == 'FO' ] || [ $p == 'F' ]; then
		FREQ=$(egrep -A 21 "State frequencies:" ${PRM} | tail -20 | awk '{print $3}'| paste -s -d ',')
		FR=$(echo "F{"$FREQ"}")
		REC_MDL="$REC_MDL+$FR";
		break
	elif [ $p == 'FQ' ]; then
		REC_MDL="$REC_MDL+FQ";
		break
	fi
done

for p in ${ARR[@]}; do
	if [ $p == 'I' ]; then
		INV=$p;
		PROP=$(egrep "Proportion of invariable sites:" ${PRM} | sed 's/[^0-9.]*//g');
		echo "$PROP";
		REC_MDL="$REC_MDL+$INV";
		break
	else
		PROP=-1
	fi
done

for p in ${ARR[@]}; do
	if [[ $p == G* ]]; then
		GAM=$p;
		ALPHA=$(egrep "alpha" ${PRM} | sed 's/[^0-9.]*//g');
		REC_MDL="$REC_MDL+$GAM";
		break
	else
		ALPHA=-1
	fi
done
		
echo "$REC_MDL";
echo "$PROP";
echo "$ALPHA";

if [ $PROP != -1 ] && [ $ALPHA != -1 ]; then
	nice iqtree -s ${ALN} -te ${TRE} -st AA -nt 4 -m $REC_MDL -a $ALPHA -i $PROP -blfix -asr;
elif [ $PROP != -1 ]; then
	nice iqtree -s ${ALN} -te ${TRE} -st AA -nt 4 -m $REC_MDL -i $PROP -blfix -asr;
elif [ $ALPHA != -1 ]; then
	nice iqtree -s ${ALN} -te ${TRE} -st AA -nt 4 -m $REC_MDL -a $ALPHA -blfix -asr;
else
	nice iqtree -s ${ALN} -te ${TRE} -st AA -nt 4 -m $REC_MDL -blfix -asr;
fi


