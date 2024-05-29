#!/bin/bash

## PoolSNP

## Author: Martin Kapun
## test if Shell is indeed BASH

###############################################
######### READ COMMANDLINE ARGUMENTS ##########
###############################################

for i in "$@"
do
  case $i in
    mpileup=*)
      mpileup="${i#*=}"
      ;;
    reference=*)
      ref="${i#*=}"
      ;;
    min-cov=*)
      mic="${i#*=}"
      ;;
    max-cov=*)
      mac="${i#*=}"
      ;;
    min-count=*)
      mico="${i#*=}"
      ;;
    min-freq=*)
      mif="${i#*=}"
      ;;
    miss-frac=*)
      mis="${i#*=}"
      ;;
    base-quality=*)
      bq="${i#*=}"
      ;;
    names=*)
      names="${i#*=}"
      ;;
    jobs=*)
      jobs="${i#*=}"
      ;;
    output=*)
      out="${i#*=}"
      ;;
    badsites=*)
      BS="${i#*=}"
      ;;
    allsites=*)
      AS="${i#*=}"
      ;;
    *)
      # unknown option
      ;;
  esac
done

###############################################
######## TEST IF GNU PARALLEL INSTALLED #######
###############################################

echo "Test if GNU parallel installed"
echo "*********************"
echo ""

command -v parallel >/dev/null 2>&1 || { echo >&2 "PoolSNP requires GNU parallel but it's not installed. Please check here: https://www.gnu.org/software/parallel/ for more details. Aborting."; exit 1; }

echo "done"
echo ""
###############################################
######## TEST IF ALL PARAMETERS ARE SET #######
###############################################

mpileup=${mpileup}
out=${out}
ref=${ref}
names=${names}

help='''
********************************
************ HELP **************
********************************


PoolSNP v. 1.05 - 13/11/2017

A typcial command line looks like this:

sh ~/PoolSNP.sh \
  mpileup=~/input.mpileup \       ## The input mpileup
output=~/output \               ## The output prefix
reference=~/reference.fa \      ## The reference FASTA file
names=sample1,sample2,sample3 \ ## A comma separated list of samples names according to the order in the mpileup file
min-cov=10 \                    ## sample-wise minimum coverage
max-cov=0.95 \                    ## Either the maximum coverage percentile to be computed or an input file
min-count=20 \                  ## minimum alternative allele count across all populations pooled
min-freq=0.001 \                ## minimum alternative allele frequency across all populations pooled
miss-frac=0.1 \                 ## maximum allowed fraction of samples not fullfilling all parameters
base-quality 15 \               ## minimum base quality for every nucleotide
jobs=10                         ## number of parallel jobs/cores used for the SNP calling

Please see below, which parameter is missing:
*******************************
'''

if [ -z "$mpileup" ]; then echo "$help\nmpileup=~/input.mpileup is missing: The full path to the input mpileup needs to be specified"; exit 1 ; fi
if [ -z "$out" ]; then echo "$help\noutput=~/output is missing: The full path to the output file prefix needs to be defined. The script will automatically append file extensions"; exit 2 ; fi
if [ -z "$names" ]; then echo "$help\nComma separated list of names correpsonding to the samples in the mpileup file"; exit 3; fi
if [ -z "$ref" ]; then echo "$help\nThe full path to the reference FASTA genome is missing"; exit 3; fi
if [ -z "$mac" ]; then echo "$help\nEither provide an upper quantile cutoff (e.g. 0.95 to exclude all reads lager than the 95% coverage percentile) or input file"; exit 4; fi
if [ -z "$mic" ]; then mic=10; fi
if [ -z "$mis" ]; then mis=0.1; fi
if [ -z "$mif" ]; then mif=0.01; fi
if [ -z "$mico" ]; then mico=20; fi
if [ -z "$bq" ]; then bq=15; fi
if [ -z "$jobs" ]; then jobs=1; fi
if [ -z "$BS" ]; then BS=1; fi
if [ -z "$AS" ]; then AS=0; fi

## change to home directory of scripts
BASEDIR=$(dirname $0)
cd $BASEDIR

## replace commas with tabs of name list
new=$(echo $names | tr ',' '\t')

## create temporary directory
mkdir -p $out/temp

###############################################
############### CREATE VCF HEADER #############
###############################################

## create VCF header
echo -e """##fileformat=VCFv4.2
##fileDate=$(date +%d'/'%m'/'%y)
##Source=PoolSnp-1.05
##Parameters=<ID=MinCov,Number=$mic,Type=Integer,Description=\"Minimum coverage per sample\">
##Parameters=<ID=MaxCov,Number=$mac,Type=Integer,Description=\"Maximum chromosome- and sample-specific maximum coverage; Either a precomuted file or the maximum percentile cutoff, eg. 0.95 to consider only reads within the 95% coverage percentile\">
##Parameters=<ID=MinCount,Number=$mico,Type=Integer,Description=\"Minimum alternative allele count across all samples pooled\">
##Parameters=<ID=MinFreq,Number=$mif,Type=Float,Description=\"Minimum alternative allele frequency across all samples pooled\">
##Parameters=<ID=MaximumMissingFraction,Number=$mis,Type=Float,Description=\"Maximum fraction of samples allowed that are not fullfilling all parameters\">
##Parameters=<ID=BaseQual,Number=$bq,Type=Integer,Description=\"Minimum PHRED scaled base quality\">
##Reference=$ref
##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >=$bq\">
##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference Counts\">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Counts\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=FREQ,Number=1,Type=Float,Description=\"Variant allele frequency\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$new""" > $out/temp/header.txt

###############################################
#### calculate max coverage when necessary ####
###############################################

# test if $mac is not a file else continue
if ! [ -f $mac ]

then

# if $mac is a decimal, calculate the coverage threshold for each chromsomal arm and sample smaller than the cutoff percentile, for example, if you provide 0.95, the maximum coverage of the 95% coverage percentile will be calculated. This is REALLY computationally demanding. Thus, to parallelize to the extend possible, I split the mpileup by contigs and execute these in parallel.

mkdir -p $out/temp/cov/mpileups

echo "Splitting Mpileup by contigs"
echo "*********************"
echo ""


# test if mpileup is gzipped and split the mpileup by contig name. AWK is the fastest, I tried a lot!!
if [[ $mpileup =~ \.gz$ ]]

then

gunzip -c $mpileup | awk -v out="$out/temp/cov/mpileups/" -v file="" -v i="" '{if (i!=$1) {close(file); close(command); if (i!="") {print i" done"}} else {file=out$1".mpileup.gz" ; command= "gzip > "file  ; print | command; }; i=$1 }'

else

awk -v out="$out/temp/cov/mpileups/" -v file="" -v i="" '{if (i!=$1) {close(file); close(command); if (i!="") {print i" done"}} else {file=out$1".mpileup.gz" ; command= "gzip > "file  ; print | command; }; i=$1 }' $mpileup

fi

mkdir -p $out/temp/cov/cutoffs

## read FASTA contig names and run Python script "max-cov.py" to calculate coverage thresholds

## Note that ALL contigs with special characters (|/\,) will be ignored

# test if Reference FASTA is gzipped

if [[ $ref =~ \.gz$ ]]

then

gunzip -c $ref | grep '^>' | awk '$1!~/\|/ && $1!~/\// && $1!~/\\/ && $1!~/,/' | awk -v inp="$out" -v cut="$mac" ' {print "python3 scripts/max-cov.py --mpileup "inp"/temp/cov/mpileups/"substr($1,2)".mpileup.gz --cutoff "cut" --contig "substr($1,2)" --out "inp"/temp/cov/cutoffs/"substr($1,2)".txt" }' > $out/temp/contignames.txt

else

grep '^>' $ref | awk '$1!~/\|/ && $1!~/\// && $1!~/\\/ && $1!~/,/' | awk -v inp="$out" -v cut="$mac" '{print "python3 scripts/max-cov.py --mpileup "inp"/temp/cov/mpileups/"substr($1,2)".mpileup.gz --cutoff "cut" --contig "substr($1,2)" --out "inp"/temp/cov/cutoffs/"substr($1,2)".txt" }' > $out/temp/contignames.txt

fi


echo "done"
echo ""
echo "Max coverage percentile $mac calculation started"
echo "*********************"
echo ""

## start parallel calulation of maximum coverage

cat $out/temp/contignames.txt | parallel --no-notice --jobs $jobs

# merge files into single coverage file
for f in $out/temp/cov/cutoffs/*
do
cat "$f" >> $out-cov-$mac.txt
done

echo "done"
echo ""

## once the calulations are done, point to the new coverage cutoff file
mac=$out-cov-$mac.txt

fi

###############################################
################## RUN POOLSNP ################
###############################################

echo "PoolSNP started"
echo "*********************"
echo ""

# test if Mpileup is comressed

if [[ $mpileup =~ \.gz$ ]]

then

## run PoolSNP with GNU parallel
gunzip -c $mpileup | parallel \
-k \
--pipe \
-j $jobs \
--no-notice \
--cat python3 scripts/PoolSnp.py \
--mpileup  {} \
--min-cov $mic \
--max-cov  $mac \
--min-freq  $mif \
--miss-frac  $mis \
--min-count  $mico \
--base-quality  $bq \
--allsites $AS \
>  $out/temp/SNPs.txt

else

parallel \
-k \
-j  $jobs \
--pipepart \
--no-notice \
-a  $mpileup \
--cat python3 scripts/PoolSnp.py \
--mpileup  {} \
--min-cov $mic \
--max-cov  $mac \
--min-freq  $mif \
--miss-frac  $mis \
--min-count  $mico \
--base-quality  $bq \
--allsites $AS \
>  $out/temp/SNPs.txt
fi

echo "done"
echo ""

###############################################
##### concatenate and remove temp files #######
###############################################

cat $out/temp/header.txt $out/temp/SNPs.txt | gzip > $out.vcf.gz

rm -r $out/temp

###############################################
## export file containing sites with bad cov ##
###############################################


if [[ $BS = 1 ]]

echo "BadSites caluclation started"
echo "*********************"
echo ""


then

if [[ $mpileup =~ \.gz$ ]]

then

gunzip -c $mpileup | parallel \
-k \
--pipe \
-j $jobs \
--no-notice \
--cat python3 scripts/bad-sites.py \
--mpileup  {} \
--min-cov $mic \
--max-cov  $mac \
--base-quality  $bq \
| gzip >  $out\_BS.txt.gz

else

parallel \
-k \
-j  $jobs \
--pipepart \
--no-notice \
-a  $mpileup \
--cat python3 scripts/bad-sites.py \
--mpileup  {} \
--min-cov $mic \
--max-cov  $mac \
--base-quality  $bq \
| gzip >  $out\_BS.txt.gz

fi

echo "done"
echo ""

fi
