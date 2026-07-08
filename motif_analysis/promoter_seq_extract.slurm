#!/bin/bash
prefix=$1
bedtools getfasta -fi mm10.fa -bed ${prefix}.bed -s -name > ${prefix}.fa
