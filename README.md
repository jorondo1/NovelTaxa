This workflow quickly identifies MAGs representing species that have no representative genome in GTDB, referred to as novel species-level MAGs (nMAGs). It adds these MAGs to the reference database before estimating community composition for of set of metagenomes. Of note, this worflow does not provide annotations of nMAGs, which can be achieved independently through GTDB-Tk.

1. Quality assessment with `CheckM 2`
2. Full pairwise comparisons with GTDB reference database with `skANI`
3. nMAGs identification : Quality score > 0.50 with no ANI >= 95%
4. Community composition estimation of metagenomic samples with `Sourmash gather`
5. Containment and alpha-diversity comparisons before / after adding nMAGs with `ggplot2`

Example usage :

```
export MAG_DIR=$ILAFORES/analysis/projet_PROVID19/Saliva/MAG_analysis/drep_genomes/dereplicated_genomes
export SAM_LIST=$ILAFORES/analysis/projet_PROVID19/Saliva/preproc/preprocessed_reads.sample.tsv
bash scripts/Novelty_minimal_pipeline.sh -o Provid19 -g r214 -m $MAG_DIR -s $SAM_LIST
```
