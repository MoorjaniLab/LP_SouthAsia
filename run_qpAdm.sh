#!/bin/sh

##Script to create a parfile and run qpAdm 
##Software part of the ADMIXTOOLS package, which can be download here: https://github.com/DReichLab/AdmixTools

ind=$1 #what individual, or population, on which you would like to run this analysis
file=$2 #eigenstrat file - more about the format : https://reich.hms.harvard.edu/software/InputFileFormats


  par="par_${ind}" #defining the parfile name
	list="list_${ind}.txt" #list of the individuals to test and the source populations - in this case, we use Indus Periphery West (IPWest), Central Steppe MLBA, and Onge (as SAHG proxy)
  echo -e "${ind}\nIPWest\nCentral_Steppe_MLBA\nOnge" > ${list}
	
	echo "DIR: /your_directoyr" >${par}
  echo "S1: ${file}" >>${par}
  echo "indivname: S1.ind" >>${par}
	echo "snpname: DIR/S1.snp" >>${par}
  echo "genotypename: DIR/S1.geno" >>${par}
  echo "popright: right_pop_list.txt" >>${par} #list of outgroup populations, or "right populations"
	echo "popleft: ${list}" >>${par}
	echo "hashcheck: NO" >>${par}
	echo "details: YES" >>${par}	
	echo "allsnps: NO" >>${par} 
	echo "inbreed: NO" >>${par}	
	
  AdmixTools/bin/qpAdm -p ${par}
  echo "done" ${par}
