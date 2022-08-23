# Extant Sequence Reconstruction (ESR)
Extant sequence reconstruction (ESR) is a workflow to calculate the conditional probability distribution of an extant sequence given an alignment, a phylogenetic tree, and an evolutionary model.

ESR works by applying a leave-one-out cross-validation to each sequence in a phylogenetic tree.  

<img width="468" alt="image" src="https://user-images.githubusercontent.com/111892527/186263175-50b87311-8f82-41c4-97ca-de61cababddd.png">

In short, extant nodes are treated like ancestral nodes. A conditional probability distribution is calculated for an extant node and a sequence can be generated. This reconstructed sequence and its properties can be compared to the true sequence and its properties.

This is useful for researchers who use ancestral sequence reconstruction (ASR) and have ever wondered whether a model parameter is actually improving your ancestral sequence reconstructions. ESR provides a mechanism for researchers to evaluate sequence reconstructions, since you can compare a sequence reconstruction to the true sequence. For example, we have demonstrated that sequence reconstructions from a model using the LG substitution matrix results in incorrect residues that are more chemically similar to the true residue than if the model uses a Poisson substitution matrix.

The goal of ESR is tool to provide a gold-standard to benchmark sequence reconstructions.

