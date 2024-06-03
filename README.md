export MAG_DIR=$ILAFORES/analysis/projet_PROVID19/Saliva/MAG_analysis/drep_genomes/dereplicated_genomes
export SAMPLE_DIR=$ILAFORES/analysis/projet_PROVID19/Saliva/preproc
bash scripts/Novelty_minimal_pipeline.sh -o Provid19 -g r214 -m $MAG_DIR -s $SAMPLE_DIR
