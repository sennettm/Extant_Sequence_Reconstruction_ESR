#!/bin/env bash
#####################################################################################################
#Script Name    : ESR_v0.1.sh
#Description    : Given an aligned set of sequences, a phylogenetic tree, and an evolutionary model
#                 a conditional probability distribution for each extant sequence will be calculated
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

#####################################################################################################

usage() {
	echo "Usage: $0 [ -b branch length ] [ -k txt file of taxa to keep ] \
	< -a alignment .a2m format > < -t tree > < -i .iqtree file > [-n number of cores available]" \
		1>&2
}

exit_abnormal() {
	usage
	exit 1
}

#####################################################################################################
#Input Variables:  b/BL    is the branch length that will be pruned, 
#                          if 0 then no branches are trimmed
#                          default is 0
#                  k/TXKP  are the taxa that will be kept if the tree is pruned
#                          default is NONE
#                  a/ALGN  is the alignment that will be reconstructed,
#                          alignment must be in .a2m or .fasta format
#                  t/TREE  is the phylogenetic tree associated with the alignment,
#                          treefile must be in newick format
#                  i/PARAM is the .iqtree file that contains the model estimates
#                  n/NDS   the total number of threads available to run analysis,
#                          so the process can be  parallelized
#####################################################################################################

BL=0;

while getopts ":b:k:a:t:h:i:n:" options; do
	case "${options}" in
		b)
			BL=${OPTARG}
			;;
		k)
			TXKP=${OPTARG:-NONE}
			;;
		a)
			ALGN=${OPTARG}
			;;
		t)
			TREE=${OPTARG}
			;;
		i)
			PARAM=${OPTARG}
			;;
		n)
			NDS=${OPTARG}
			;;
		h)
			usage
			;;
		*)
			exit_abnormal
			;;
	esac
done


echo "$BL" "$TXKP" "$ALGN" "$TREE" "$PARAM" "$NDS";
sleep 2;

#####################################################################################################
#
#remove any previously made folders/files 
#
#####################################################################################################
FLD=$(egrep ">" $ALGN | sed 's/>//');

rm -rf ${FLD};

rm -rf *.fst *.tre *.out *.phy *_${BL}.a2m;

rm -rf *_node_*.txt;
#####################################################################################################
#
#modify tree and alignment for processing
#
#####################################################################################################
#prune the tree of interest by branch length = BL
./SilNode.py ${TREE};
./prune_leaf_nodes2.py ${TREE} ${BL} ${TXKP};

#rename the pruned tree to incorporate the BL trim
cp ${TREE}_pruned.tree ${TREE}_pruned_${BL}.tree;
PTREE=${TREE}_pruned_${BL}.tree;

#create all the trees required for extant reconstruction
./infer_leaf_node.py ${PTREE};

#create a list of all the taxa in the pruned tree
./print_leaf_taxa.py ${PTREE} | tee taxa_${BL}.txt

#make sure the names are all the same
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1' taxa_${BL}.txt

#rename the alignment to correspond to the BL
#seqcon -E ${ALGN};
./ms_seqcon.py --MSA ${ALGN} --form1 fasta -E;
touch ${ALGN%.a2m}_${BL}.a2m;
PALGN=${ALGN%.a2m}_${BL}.a2m;
for TAXA in *_recon.tre; do
	cat ${TAXA%_recon.tre}.fst >> ${PALGN}
done

#remove any gap only columns and format as fasta
#seqcon -d -f fasta ${PALGN};
PALGN=$(./ms_seqcon.py --MSA ${PALGN} --form1 fasta -D);
#PALGN=${PALGN%.a2m}_del.fasta;
#./ms_seqcon.py --MSA ${PALGN} --form1 fasta -R;

#rename the algnment to rm _sc from name
#cp ${PALGN%.a2m}_sc.a2m ${PALGN};

#reformat alignment to phyml format
#seqcon -f phyml ${PALGN};
./ms_seqcon.py --MSA ${PALGN} --form1 fasta --form2 phylip -R;

#rename the algnment to remove _sc
#cp ${PALGN%.a2m}_sc.phy ${PALGN%.a2m}.phy;
#PHYALGN=${PALGN%.a2m}.phy
PHYALGN=${PALGN%.fasta}.phylip

#determine the number of taxa and the sequence length
head -1 ${PHYALGN}

#creates and appends "PSUEDO" to .fasta file that will be needed for iqtree
echo ">PSEUDO" >> ${PALGN}


#append a L to algnmt file so we have a PSEUDO followed by poly leucine sequence
#each tree from infer_leaf_node.py script has a "PSEUDO" daughter node
SEQLEN=$(head -1 ${PHYALGN} | awk '{print $2}')
for (( i=0; i < $SEQLEN; ++i )); do
    printf "%c" "L" >> ${PALGN}
done

echo >> ${PALGN}
#####################################################################################################
#
#Perform ESR on modified alignments and trees, then pull some useful information
#
#####################################################################################################
#./infer_leaves.sh baliphyalignmentwithancs.fasta
./ms_infer_leaves.sh ${PALGN} ${PTREE} ${PARAM} ${NDS};

#make fasta algnt w/ gapped recons
./ms_crunch_all_leaves.sh ${BL} ${PALGN};

#calc the seq ID and avg pp for each sequence
./ms_probs.sh ${BL} ${PALGN} | tee probs_${BL}.txt;

#calculate the expected probability for each reconstructed sequence
./eLnP_extant_seqs.sh ${BL} ${PALGN};

#calculate the conditional probability for each true sequence
./LnP_extant_seqs.sh ${BL} ${PALGN};

#calculate the conditional probability for each SMP sequence
./LnP_SMP_seqs.sh ${BL};
#####################################################################################################
#
# clean up
#
#####################################################################################################
rm matlab_* twist_* node_*;

for TAXFILE in *_recon.tre; do
    TAXA=${TAXFILE%_recon.tre};
    mv ${TAXA}_recon.tre ${TAXA}/.;
    mv ${TAXA}_*_ancprobs.out ${TAXA}/.;
    mv ${TAXA}_*_probs.txt ${TAXA}/.;
    mv ${TAXA}.fst ${TAXA}/.;
done
#####################################################################################################

echo "Done with ESR";
