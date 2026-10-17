#!/bin/sh
##script to run vcftools to estimate the allele frequencies in South Asian dataset

file=$1 #vcf file
ind=$2 #list of individuals 

vcftools --vcf ${file} --freq --keep ${ind} --out ${file}_out #from the outfile grep the snp of interest to get allele frequency
