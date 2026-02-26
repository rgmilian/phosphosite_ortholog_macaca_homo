# phosphosite_ortholog_macaca_homo
Identification of Human Orthologous Phosphorylation Sites

To identify human orthologous phosphorylation sites corresponding to Macaca mulatta phosphoproteomics data, we used the following bioinformatic pipeline integrating ortholog mapping and sequence alignment. 
For each macaque protein, we queried the UniProt REST API to retrieve gene name annotations and identify the corresponding human ortholog (UniProt Consortium, 2025). Human orthologs were identified by searching for reviewed Swiss-Prot entries (organism ID: 9606) matching the macaque gene name.

Full-length protein sequences for both macaque and human orthologs were retrieved from UniProt in FASTA format. To map phosphorylation site positions between species, we performed global pairwise sequence alignment using the Biopython Bio.Align.PairwiseAligner module (Cock et al., 2009). The peptide sequences containing phosphorylated residues were first localized within the full-length macaque protein sequence. The corresponding amino acid position in the macaque protein was then mapped to the human ortholog through the global alignment, accounting for gaps and sequence variations. The aligned human position was verified by confirming the presence of the same phosphorylatable amino acid (serine, threonine, or tyrosine) at the mapped location.

This approach enabled systematic identification of conserved phosphorylation sites across species while accounting for insertions, deletions, and sequence divergence between macaque and human proteins. All bioinformatic analyses were performed in Python using pandas (v3.13.5), Biopython (v1.86), and the requests library for API queries.

References

The UniProt Consortium , UniProt: the Universal Protein Knowledgebase in 2025, Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D609–D617, https://doi.org/10.1093/nar/gkae1010

Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. Epub 2009 Mar 20. PMID: 19304878; PMCID: PMC2682512.
