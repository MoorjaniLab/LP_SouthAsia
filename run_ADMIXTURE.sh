#!/bin/sh
##Script to run ADMIXTURE 
#softwarde can be download here : https://dalexander.github.io/admixture/download.html

file=$1 #in eigenstrat format

for K in {2..8} #run ADMIXTURE for different K
do
  admixture_linux-1.3.0/admixture --cv=10 ${file}".bed" ${K} 
done
