
assembly_stats() {
MAG_DIR="$PARENT_DIR/MAG_analysis"
DREP="$MAG_DIR/drep_genomes"
find $DREP/dereplicated_genomes/ -type f -name '*.fa' | \
xargs -n1 basename | \
grep -wf - $DREP/data/checkM/checkM_outdir/results.tsv > $DREP/checkm_final_set.tsv
                
# check average completeness & contamination
MEAN_COMP=$(cut -f12 $DREP/checkm_final_set.tsv | awk '{x+=$0}END{print x/NR}')
MEAN_CONT=$(cut -f13 $DREP/checkm_final_set.tsv | awk '{x+=$0}END{print x/NR}')

# sd completeness and contamination
SD_COMP=$(cut -f12 $DREP/checkm_final_set.tsv | awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}')
SD_CONT=$(cut -f13 $DREP/checkm_final_set.tsv | awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}')

echo "mean $MEAN_COMP sd $SD_COMP completion"
echo "mean $MEAN_CONT sd $SD_CONT contamination"
echo "$(tail -n +2 $DREP/checkm_final_set.tsv | wc | awk '{print $1}') total bins"
find $DREP/dereplicated_genomes -type f -name "*.fa" > $DREP/bins_list.txt

}

# Concatenate samples 
cat_samples() {
    local metadata_file="$1"

    # Extract unique combinations of COM and MS from metadata
    cut -f2,3 "$metadata_file" | sort -u > unique_combinations.txt

    # Loop through the unique combinations and concatenate the files
    while read -r COM MS; do
        # Create an array to store the filenames matching the combination
        file_array_1=()
        file_array_2=()
        # Loop through SAM values to find files matching the combination
        while read -r SAM; do
            file_array_1+=( "preproc/$SAM/${SAM}_paired_1.fastq" )
            file_array_2+=( "preproc/$SAM/${SAM}_paired_2.fastq" )
        done < <(grep -P "$COM\t$MS" "$metadata_file" | cut -f1)

        # Concatenate the matching files to create the output file
        cat "${file_array_1[@]}" > cat_reads/"${COM}_${MS}_1.fastq"
        cat "${file_array_2[@]}" > cat_reads/"${COM}_${MS}_2.fastq"
    done < unique_combinations.txt
}

novel_genomes() {
MAG_DIR="$PARENT_DIR/MAG_analysis"
DREP="$MAG_DIR/drep_genomes"

export GTDB_TSV=$MAG_DIR/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv

cat $MAG_DIR/drep_genomes/checkm_final_set.tsv | cut -f1,12,13 | \
        awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $2 - 5*$3}' | awk '$4 > 50.00' > MAGs_QS.txt

# Evaluate distances between our bins and all genomes from the list 
mash dist -p 48 -d 0.05 "mash_sketch.msh" -l $DREP/bins_list.txt > mash_dist_result.tab

# parse Mash output to identify new species 
python3 $ILAFORES/programs/parse_mash.py \
-m mash_dist_result.tab -o Mash_processed \
        -f $(echo "$DREP/dereplicated_genomes" | tr " " "\n") --discard-repeat-strains

# Filter-out non-novel species (mash distance), keep only unidentified GTDB species, 
# then compute quality score COMPLETENESS - 5 * CONTAMINATION 
find Mash_processed/New_species/ -type f -exec basename {} \; | sed 's/\.fa//' | \
        grep -wf - $GTDB_TSV | grep -w 's__' - | \
        cut -f1 | awk '{print $0 ".fa"}' | grep -wf - $MAG_DIR/drep_genomes/checkm_final_set.tsv | cut -f1,12,13 | \
        awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $2 - 5*$3}'> novel_quality_scores.txt

# Only keep >50 QS and create a list of bins and corresponding taxonomy
:> novel_MAGs.txtq
cat novel_quality_scores.txt | awk '$4 > 50.00' | \
        cut -f1 | sed 's/\.fa//' | grep -wf - $GTDB_TSV | \
        cut -f1,2 | sed 's/\t/,,/g' | sed 's/;/,/g' >> novel_MAGs.txt
}

# GTDB taxonomy has commas in it, this messes up the column recognition we need in eval_cont(). Fixing it:
fix_gtdb() {
export gather=$PWD/"${1}"
for file in $gather/*.csv; do
        awk -F'"' -v OFS='"' '{ for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done
}

eval_cont() {
	
if [[ -z ${1} ]]; then
	echo 'Missing positional argument!'
	echo '1: directory containing gather results'
	exit 1
fi

# columns of interest

local SM_DIR=$(realpath "${1}")
q_index=$(awk -v RS=',' '/query_name/{print NR; exit}' ${SM_DIR}/*_gather.csv)
c_index=$(awk -v RS=',' '/f_unique_weighted/{print NR; exit}' ${SM_DIR}/*_gather.csv)

# Compile containment by run
gather_files=($(find ${SM_DIR} -maxdepth 1 -type f -name '*_gather.csv'))
:> ${1}/cntm_sum.txt
for file in "${gather_files[@]}"; do
        RUN_CNTM=$(cut -d',' -f $c_index $file | tail -n +2 | awk '{sum+=$1;} END{print sum;}')
        RUN_ID=$(cut -d',' -f $q_index $file | head | tail -n 1)
        echo $RUN_ID
        DB_ID=$(echo $file | sed "s|.*${RUN_ID}_||" | sed 's/_gather\.csv//')
        echo -e "${RUN_ID}\t${DB_ID}\t${RUN_CNTM}" >> ${1}/cntm_sum.txt
done
sort -o ${1}/cntm_sum.txt -k1,1 -k2,2 ${1}/cntm_sum.txt  

# List dbs 
dbs=($(awk '{if (!seen[$2]++) print $2}' ${1}/cntm_sum.txt))

# Compile overall containment by db type
echo -e "ref_db\tcntm_avg\tcntm_sd"
for i in "${dbs[@]}"; do
        CNTM_AVG=$(grep -w "${i}" "${1}/cntm_sum.txt" | cut -f3 - | awk '{x+=$0}END{print x/NR}') # compute average run containment
        CNTM_SD=$(grep -w "${i}" "${1}/cntm_sum.txt" | cut -f3 - | awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}') # and standard deviation
        echo -e "$i\t${CNTM_AVG}\t${CNTM_SD}"
done
}
 
 
