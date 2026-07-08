#!/bin/bash
awk '$3=="gene"' /export/genomic/mm10/gencode.vM23.annotation.gtf | \
awk 'BEGIN{OFS="\t"}
{
gene_name=""; match($0, /gene_name "([^"]+)"/, a); gene_name=a[1];
if($7=="+"){print $1,$4-1,$4,gene_name,$6,$7}
else if($7=="-"){print $1,$5-1,$5,gene_name,$6,$7}}' > mm10_TSS.bed