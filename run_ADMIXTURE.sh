#!/bin/sh
##Script to run ADMIXTURE 
#softwarde can be download here : https://dalexander.github.io/admixture/download.html

file=$1 #in eigenstrat format, input the file when posting this run script

for K in {2..8} #run ADMIXTURE for K+2 through 8 with 10-fold cross validation
do
  admixture_linux-1.3.0/admixture --cv=10 ${file}.bed ${K} 
done
