#!/bin/bash

# https://meme-suite.org/meme/doc/ame.html

ALL_PROMOTERS="mm10_promoter_500bp.fa"
MODULE_CSV="HM_ModGene.csv"
EVAL_THRESHOLD="1"

mkdir -p mod_fasta
mkdir -p AME_out
# Module ID
modules=$(awk -F',' 'NR>1 {print $1}' $MODULE_CSV | sort | uniq)  # skip header!

for mod in $modules; do
    echo "Processing Module: $mod"
    awk -F',' -v m="$mod" '$1 == m {print $2}' $MODULE_CSV > mod_fasta/Module${mod}.genes
    seqkit grep -i --id-regexp "^(.+?)::" -f mod_fasta/Module${mod}.genes $ALL_PROMOTERS -o mod_fasta/Module${mod}.fa
done

for mod in $modules; do  
    if [[ "$mod" == "0" ]]; then
        continue
    fi
    # Do not conduct motif analysis if there are too few gene
    count=$(grep -c ">" mod_fasta/Module${mod}.fa)
    if [ "$count" -lt 50 ]; then
        echo "Skipping $mod: Too few sequences ($count)."
        continue
    fi
    echo "Submitting Module: $mod"
        # Do not conduct motif analysis if there are too few genes
    sbatch AME.slurm "$mod" "mod_fasta/Module${mod}.fa" "$EVAL_THRESHOLD" "AME_out"
done

echo "All modules processed and jobs submitted."