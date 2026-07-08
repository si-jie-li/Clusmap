#!/bin/bash
samtools faidx mm10.fa
bp=$1
bedtools slop -i mm10_TSS.bed -g mm10.fa.fai -b $bp> mm10_promoter_${bp}bp.bed