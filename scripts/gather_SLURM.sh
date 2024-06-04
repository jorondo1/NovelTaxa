#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/NovelTaxa/
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/NovelTaxa/logs/sourmash-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=15G
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -A def-ilafores
#SBATCH -J sourmash

newgrp def-ilafores
export ILAFORES=${ANCHOR}/${ILAFORES}
export OUTDIR=${ANCHOR}/${OUTDIR}/sourmash
export SM_DB=${ANCHOR}/${DB}/sourmash_db/gtdb-rs214-reps.k31.zip #atm hardcoded because not produced locally
export MAGs_IDX=${ANCHOR}/${MAGs_IDX}
export SM_SK=${ANCHOR}/${SM_SK}
export sourmash="singularity exec --writable-tmpfs -e -B ${ANCHOR}/home:${ANCHOR}/home ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"
echo "$sourmash"
echo "loading env"
module load StdEnv/2020 apptainer/1.1.5

# copy container
# cp ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.7.0.sif $tmp/

export SAM_NUM=$(cat ${ANCHOR}/${SAM_LIST} | awk "NR==$SLURM_ARRAY_TASK_ID")
export SAMPLE=$(echo -e "$SAM_NUM" | cut -f1)
export FQ_P1=$(echo -e "$SAM_NUM" | cut -f2)
export FQ_P2=$(echo -e "$SAM_NUM" | cut -f3)
export FQ_U1=$(echo -e "$SAM_NUM" | cut -f4)
export FQ_U2=$(echo -e "$SAM_NUM" | cut -f5)
export SIG="${SM_SK}/${SAMPLE}.sig"

echo "Executing pipeline on sample ${SAMPLE} !"

if [[ ! -f $SIG ]]; then
	echo "Sketch metagenomes"
$sourmash sketch dna -p k=31,scaled=1000,abund --merge ${SAMPLE} -o $SIG $FQ_P1 $FQ_P2 $FQ_U1 $FQ_U2
else
	echo "Metagenome sketches found. Skipping..."
fi
echo ${OUTDIR}/sourmash/${SAMPLE}_${GTDB_V}_gather.csv

if [[ ! -f ${OUTDIR}/${SAMPLE}_${GTDB_V}_gather.csv ]]; then
	echo "Gather against the gtdb index"
	$sourmash gather $SIG ${SM_DB} -o ${OUTDIR}/${SAMPLE}_${GTDB_V}_gather.csv
else
	echo "Gather output found. Skipping..."
fi

if [[ ! -f ${OUTDIR}/${SAMPLE}_custom_gather.csv ]]; then
	echo "Gather again but add the novel MAGs"
	$sourmash gather $SIG ${MAGs_IDX} ${SM_DB} -o ${OUTDIR}/${SAMPLE}_custom_gather.csv
else
	echo "Gather output found. Skipping"
fi

echo "Done !"