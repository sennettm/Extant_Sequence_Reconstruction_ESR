# Extant Sequence Reconstruction (ESR)
**still in beta** 
## Overview

**What is ESR?**

Extant sequence reconstruction (ESR) is a workflow to calculate the conditional probability distribution of an extant sequence given an alignment, a phylogenetic tree, and an evolutionary model.

**What is the motivation behind ESR?**

This is useful for researchers who use ancestral sequence reconstruction (ASR) and have ever wondered whether a model parameter is actually improving your ancestral sequence reconstructions. ESR provides a mechanism for researchers to evaluate sequence reconstructions, since you can compare a sequence reconstruction to the true sequence. ESR is tool to provide a gold-standard to benchmark sequence reconstructions.

**How does ESR work?**

ESR works by applying a leave-one-out cross-validation to each sequence in a phylogenetic tree and alignment.  

<img width="700" alt="image" src="https://user-images.githubusercontent.com/111892527/186263175-50b87311-8f82-41c4-97ca-de61cababddd.png">

In short, extant nodes are treated like unknown ancestral nodes. A conditional probability distribution is calculated for an extant node, from which a sequence can be generated. This reconstructed sequence and its properties can be compared to the true sequence and its properties. 

**What can we do with ESR?**

The side-chains from most probable sequence reconstructions using the LG substitution matrix are more chemically similar to the true sequence than the most probable reconstructions from the Poisson, according to the Grantham distance.

<img width="468" alt="image" src="https://user-images.githubusercontent.com/111892527/187996511-9b6ad4fe-a755-4671-8fe3-4768de32a47d.png">

Blue are Grantham distances measured per SMP mistake from reconstructions generated with the Poisson substitution matrix, whereas orange are from the LG
substitution matrix. A lower Grantham distance means the mistake is more chemically similar to the true residue.


In addition, we have used ESR to visualize the relative uncertainty of a probability distribution to demonstrate where information is lacking in a phylogenetic tree. 

<img width="468" alt="image" src="https://user-images.githubusercontent.com/111892527/186267196-de75a0f4-2dc9-4665-8c44-554634edffc0.png">

The tree is colored based on relative eLL for each node. Red indicates there is relatively more uncertainty and blue indicates relatively less uncertainty.  

## Implementing

**Requirements for ESR**

IQ-Tree: http://www.iqtree.org/#download
```
sudo apt-get install iqtree
```
***Please note ESR is currently only compatible with IQ-Tree files.***

Dendropy: https://dendropy.org/downloading.html
```
python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git
```

NumPy: https://numpy.org/install/
```
pip install numpy
```

Biopython: https://biopython.org/wiki/Download
```
pip install biopython
```

Download the ESR_script folder and Test_Data folder.

## Output

realprobs_#.txt contains the actual fraction correct  
aveprobs_#.txt contains the average probability of each SMP sequence  
tru_LnP_#.txt contains the true sequence ln-probability 
SMP_LnP_#.txt contains the SMP sequence ln-probability  
eLnP_#.txt contains the expected ln-probability of the sequence distribution  
A folder containing the iqtree output files, including the .state file which contains the sequence distribution  
something_recon_0.a2m & something_extant_0.a2m, are MSA files of the SMP and true sequences respectively  

## Running

Copy the scripts and one of the test data sets into some test folder. 

Usage
```
./ESR_v0.1.sh -h 
```

Example: Perform ESR on an Apicomplexan L/MDH dataset on a computer with with 8 threads
```
./ESR_v0.1.sh -a Apico2020_seqs.fasta_mafft -t Apico2020_seqs.fasta_mafft.treefile -i Apico2020_seqs.fasta_mafft.iqtree -n 8
```

## Common Errors

Sequence names cannot contain ".1", as is common in an accession number. You will have to replace all "." with "_".

If sequences contain very short branch lengths, then you may encounter the following:
```
.././get_LL.py:69: RuntimeWarning: divide by zero encountered in log
  ln_prob=np.log(float(prob))
  ```
You will have to remove the offending sequence.
