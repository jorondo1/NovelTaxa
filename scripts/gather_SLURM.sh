#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/NovelSpecies/Moss_test
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/NovelSpecies/Moss_test/sourmash/logs/sourmash-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=15G
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -A def-ilafores
#SBATCH -J sourmash

newgrp def-ilafores
export ILAFORES=${ANCHOR}/${ILAFORES}
export OUTDIR=${ANCHOR}/${OUTDIR}
export SM_DB=${ANCHOR}/${DB}/sourmash_db/gtdb-rs214-reps.k31.zip
export MAGs_IDX=${ANCHOR}/${MAGs_IDX}
export SAMPLE_DIR=${ANCHOR}/${SAMPLE_DIR}
export sourmash="singularity exec --writable-tmpfs -e -B ${ANCHOR}/home:${ANCHOR}/home ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"
echo "$sourmash"
echo "loading env"
module load StdEnv/2020 apptainer/1.1.5

# copy container
# cp ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.7.0.sif $tmp/

# Define sample
export SAMPLE_NUM=$(cat ${SAMPLE_DIR}/clean_samples.tsv | awk "NR==$SLURM_ARRAY_TASK_ID")
export SAMPLE=$(echo -e "$SAMPLE_NUM" | cut -f1)
export SIG="${OUTDIR}/sourmash/sketches/${SAMPLE}.sig"
export FQ_DIR=${SAMPLE_DIR}/${SAMPLE}

echo "Executing pipeline on sample ${SAMPLE} !"

if [[ ! -f $SIG ]]; then
	echo "Sketch metagenomes"
$sourmash sketch dna -p k=31,scaled=1000,abund --merge ${SAMPLE} -o $SIG \
	${FQ_DIR}/${SAMPLE}_paired_1.fastq ${FQ_DIR}/${SAMPLE}_paired_2.fastq \
	${FQ_DIR}/${SAMPLE}_unmatched_1.fastq ${FQ_DIR}/${SAMPLE}_unmatched_2.fastq
else
	echo "Metagenome sketches found. Skipping..."
fi

if [[ ! -f ${OUTDIR}/sourmash/${SAMPLE}_${$GTDB_V}_gather.csv ]]; then
	echo "Gather against the gtdb index"
	$sourmash gather $SIG ${SM_DB} -o ${OUTDIR}/sourmash/${SAMPLE}_{$GTDB_V}_gather.csv
else
	echo "Gather output found. Skipping..."
fi

if [[ ! -f ${OUTDIR}/sourmash/${SAMPLE}_custom_gather.csv ]]; then
	echo "Gather again but add the novel MAGs"
	$sourmash gather $SIG ${MAGs_IDX} ${SM_DB} -o ${OUTDIR}/sourmash/${SAMPLE}_custom_gather.csv
else
	echo "Gather output found. Skipping"
fi

echo "Done !"