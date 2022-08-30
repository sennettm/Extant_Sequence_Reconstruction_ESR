# Extant Sequence Reconstruction (ESR)

## Overview

**What is ESR?**

Extant sequence reconstruction (ESR) is a workflow to calculate the conditional probability distribution of an extant sequence given an alignment, a phylogenetic tree, and an evolutionary model.

**What is the motivation behind ESR?**

This is useful for researchers who use ancestral sequence reconstruction (ASR) and have ever wondered whether a model parameter is actually improving your ancestral sequence reconstructions. ESR provides a mechanism for researchers to evaluate sequence reconstructions, since you can compare a sequence reconstruction to the true sequence. ESR is tool to provide a gold-standard to benchmark sequence reconstructions.

**How does ESR work?**

ESR works by applying a leave-one-out cross-validation to each sequence in a phylogenetic tree and alignment.  

<img width="468" alt="image" src="https://user-images.githubusercontent.com/111892527/186263175-50b87311-8f82-41c4-97ca-de61cababddd.png">

In short, extant nodes are treated like unknown ancestral nodes. A conditional probability distribution is calculated for an extant node, from which a sequence can be generated. This reconstructed sequence and its properties can be compared to the true sequence and its properties. 

**What can we do with ESR?**

We have determined that using the LG substitution matrix instead of the Poisson substitution matrix may slightly increase the number of incorrect residues in the most probable sequence reconstructions, but the most probable sequence reconstructions from LG are still more chemically similar to the true sequence than the most probable reconstructions from the Poisson.

In addition, we have used ESR to visualize the relative uncertainty of a probability distribution to demonstrate where information is lacking in a phylogenetic tree. 

<img width="468" alt="image" src="https://user-images.githubusercontent.com/111892527/186267196-de75a0f4-2dc9-4665-8c44-554634edffc0.png">

To improve the likelihood of this tree we should include more closely related sequences to those colored in red, since it will help reduce the uncertainty of sequences in those clades.

## Implementing

**Requirements for ESR**

IQ-Tree: http://www.iqtree.org/#download
```
sudo apt-get install iqtree
```

Dendropy: https://dendropy.org/downloading.html
```
python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git
```

NumPy: https://numpy.org/install/
```
pip install numpy
```

seqcon: 

ancprobs:

## Running

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
