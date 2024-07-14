This workflow quickly identifies MAGs representing species that have no representative genome in GTDB, referred to as novel species-level MAGs (nMAGs). It adds these MAGs to the reference database before estimating community composition for of set of metagenomes. Of note, this worflow does not provide taxonomic annotations of nMAGs, which can be achieved independently through GTDB-Tk.

Looking for the [Poster presented at ISMB24](figures/ronj2303_ISMB2024.pdf) ? Special thanks to Sarah Ishak for the immense sample collection & extraction effort. We have [a preprint](https://www.biorxiv.org/content/10.1101/2023.04.06.535926v2) out that uses the methods presented here, check it out. 

1. Quality assessment with `CheckM 2`
2. Full pairwise comparisons with GTDB reference database with `skANI`
3. nMAGs identification: Quality score (QS) > 0.50 without ANI >= 95%
4. Better MAGs identification: MAGs with ANI >= 95% are kept if best reference genome has lower QS
5. nMAGs added to reference database and better MAGs substituted for reference genome
6. Community composition estimation of metagenomic samples with `Sourmash gather`
7. Containment and alpha-diversity comparisons between default vs. improved reference database `ggplot2`
