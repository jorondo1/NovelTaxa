export PARENT=${ILAFORES}/analysis/NovelSpecies
export ILL_PIPELINES=${ILAFORES}/programs/ILL_pipelines
export MAG_DIR=$ILAFORES/analysis/boreal_moss/MAG_analysis
export DB=${ILAFORES}/ref_dbs
cd $PARENT

# Create a directory with your dereplicated MAGs and their checkm output
mkdir -p MAGs out/GTDBTk
cp -r $MAG_DIR/drep_genomes/dereplicated_genomes/* MAGs/
cp -r $MAG_DIR/drep_genomes/checkm_final_set.tsv MAGs/
find MAGs -type f -name '*.fa' | sed 's#.*/##' > MAGs/List.txt

### QUALITY SCORE ASSIGNMENT
cat MAGs/checkm_final_set.tsv | cut -f1,12,13 | \
	awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $2 - 5*$3}' | awk '$4 >= 50.00' > Quality_MAGs.txt

#############################
### GTDB TAXONOMIC ASSIGNMENT
#############################
export GTDB=$ILAFORES/ref_dbs/GTDB

module load apptainer
singularity exec --writable-tmpfs --env GTDBTK_DATA_PATH=${GTDB}/release220/ \
-B /home:/home -e ${ILAFORES}/programs/ILL_pipelines/containers/gtdbtk.2.4.0.sif \
gtdbtk classify_wf --cpus 32 --genome_dir ${PARENT}/MAGs \
--out_dir ${PARENT}/out/GTDBTk --mash_db ${GTDB}/release220/ --extension fa

# Make a light version of https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $6, $7}' bac120_metadata_r220.tsv > bac120_metadata_r220_short.tsv 

###############################
### AVERAGE NUCLEOTIDE IDENTITY 
###############################

# Create a sketch of the GTDB genome reps
export skani=/home/def-ilafores/programs/skani/skani
if [[ ! -d "$GTDB"/gtdb_skani_database_ani ]]; then
	cd $GTDB
	find gtdb_genomes_reps_r220/ -name '*.fna.gz' > gtdb_file_names.txt
	"$skani" sketch -l gtdb_file_names.txt -o gtdb_skani_database_ani -t 72
	cd $PARENT
fi 

# Compute ANI
"$skani" search -d $GTDB/gtdb_skani_database_ani -o out/ANI_results.txt -t 72 --ql MAGs/List.txt 

cat out/ANI_results.txt

#################
### MASH DISTANCE USING ENSEMBLBACTERIA57
#################
module load StdEnv/2020 gcc/9.3.0 mash/2.3

if [[ ! -f $ILAFORES/EnsemblBacteria57/EnsemblBacteria57.msh ]]; then
	cd $ILAFORES/EnsemblBacteria57
	find -type f -name '*.fa.gz' | sed 's#.*/##' > genomes_list.txt
	mash sketch -p 48 -l genomes_list.txt -o EnsemblBacteria57
	cd $PARENT
fi

# Evaluate distances between our bins and all genomes from the list 
cd MAGs
mash dist -p 48 -d 0.05 $ILAFORES/EnsemblBacteria57/EnsemblBacteria57.msh -l List.txt > ../out/mash_dist.tab

# parse Mash output to identify new species 
python3 $ILAFORES/programs/parse_mash.py \
-m mash_dist_result.tab -o Mash_processed \
	-f $(echo "$DREP/dereplicated_genomes" | tr " " "\n") --discard-repeat-strains

###########
### CHECKM2
###########
mkdir -p checkm2
singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} \
${ILL_PIPELINES}/containers/checkm2.1.0.2.sif \
checkm2 predict --threads 48 \
--database_path ${DB}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
--input ${MAG_DIR} --extension .fa --output-directory "$PARENT"/out/checkm2

tail -n +2 out/checkm2/quality_report.tsv | cut -f1,2,3 | \
	awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $2 - 5*$3}' | \
	awk '$2 >= 50.00 && $4 >= 50.00 && $3 < 5.00' > Quality_MAGs_checkm2.txt



SOURMASH="/home/def-ilafores/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"
