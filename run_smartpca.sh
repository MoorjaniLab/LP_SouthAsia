#!/bin/sh

##Script to create a parfile and run SMARTPCA, part of the EIGENSOFT package
##package and software can be download here: https://github.com/DReichLab/EIG

inputfile=$1 #in eigenstrat formatl (more about the format : https://reich.hms.harvard.edu/software/InputFileFormats)
pop_to_consider=$2 #list of populations on which you would like to base your PCA
pop=${pop_to_consider}

	# create parfile
echo "DIR: ${DIR} " > par_${pop} #define your directory where your files are
echo "SSS: ${inputfile}   " >>par_${pop} #base of the eigenstrat file names
echo "S1:  ${pop} " >> par_${pop}
echo "genotypename: DIR/SSS.geno" >>par_${pop} 
echo "snpname: DIR/SSS.snp" >>par_${pop}
echo "indivname: DIR/renamed.ind" >>par_${pop}
echo "evecoutname: S1.evec" >>par_${pop}
echo "evaloutname: S1.eval" >>par_${pop}
echo "nsnpldregress: 2       " >>par_${pop}
echo "snpweightoutname: S1_SNPWeightout.txt " >>par_${pop}
echo "numoutliter: 0" >>par_${pop}
echo "numoutleigs: 10" >>par_${pop}
echo "numoutevec: 10" >>par_${pop}
echo "pubmean: YES" >>par_${pop}
echo "phylipname: S1.phyl" >>par_${pop}
echo "fastmode: NO" >>par_${pop}
echo "fstonly: NO" >>par_${pop}
echo "poplistname: ${pop_to_consider}.txt" >> par_${pop}
echo "tracefile: ${pop}.trace" >>par_${pop}
# Run smartpca
smartpca -p par_${pop} #include your path to smartpca

echo "done"
