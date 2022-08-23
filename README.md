# Extant Sequence Reconstruction (ESR)
Calculate conditional probability distributions of extant sequences given a phylogenetic tree, sequence alignment, and evolutionary model.
Extant Sequence Reconstruction
Extant sequence reconstruction (ESR) is a workflow to calculate the conditional probability distribution of an extant sequence given an alignment, a phylogenetic tree, and an evolutionary model.

ESR works by applying a leave-one-out cross-validation to each sequence in a phylogenetic tree. 

 

This is useful for researchers who use ancestral sequence reconstruction (ASR). ESR provides a mechanism for researchers to evaluate sequence reconstructions, since you can compare a sequence reconstruction to the true sequence. For example, we have demonstrated that sequence reconstructions from a model using the LG substitution matrix results in incorrect residues that are more chemically similar to the true residue than if the model uses a Poisson substitution matrix.

The goal of ESR is tool to provide a gold-standard to benchmark sequence reconstructions.
![image](https://user-images.githubusercontent.com/111892527/186262853-7f59a0fa-8e8d-41c7-9cd0-79534db39f2c.png)
