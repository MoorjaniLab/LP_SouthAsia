#!/bin/sh

for chr in {1..22}
do
        # Decode genetic map (available at https://www.decode.com/addendum ) converted in PLINK format: 
        map_plinkformat="/decode_recombinationmap/chr${chr}_hg38_genetic_map.txt" 

        ##Need phased vcf
        vcf="/PhasedVcf/chr${chr}.myvcf.gz"

        #hap-ibd can be download here : https://github.com/browning-lab/hap-ibd      
        java -jar hap-ibd.jar  gt=${vcf} map=${map_plinkformat} out=hapIBD_out_chr${chr}   min-seed=0.5  max-gap=1000 min-extend=0.5 min-output=0.5 min-markers=100 min-mac=2 nthreads=1
done
