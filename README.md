# PoolSNP

## Heuristic SNP caller for pooled sequencing data

### Requirements

[GNU parallel](https://www.gnu.org/software/parallel/) 
Python 2.7x

### Description

PoolSNP is a heuristic SNP caller, which uses an MPILEUP file and a reference genome in FASTA format as inputs. NOTE, that the FASTA headers may NOT contain any special characters, such as "/\|,:", or else they will be ignored. Heuristic parameters to be passed to the shellscript are: 

* **Minimum coverage** across all libraries; e.g. *min-cov=10* will only consider positions with a minimum coverage >10 

* **Maximum coverage** is calculated for every library and chromosomal arm as the percentile of a coverage distribution; e.g. *max-cov=0.98* will only consider positions within the 98% coverage percentile for a given sample and chromosomal arm; This cutoff can either be newly calculated by providing the threshold percentile or by providing the full path to a max-coverage file. This needs to be tab-delimited with two columns, where the first contains the chromosome name and the second a comma-separated list of coverage-thresholds corresponding in length to the number of samples in the mpileup file. For example, if a max-cov threshold of 100 and 200 should be set for the first three and the last two samples in the mpileup file for chromosome 2L, a typical line in the max-cov file should look like this: *2L 100,100,100,200,200*

* **Minimum allele count** of a minor allele across all libraries combined; e.g. *min-count=10* will only consider a position as polymorphic if the cummulative count of a minor allele across all samples in the input mpileup is equal or larger than the threshold.

* **Minimum allele frequency** of a minor allele across all libraries combined; e.g. *min-freq=10* will only consider a position as polymorphic if the relative frequency of a minor allele calculated from cummulative counts across all samples in the input mpileup is equal or larger than the threshold. 

* **Missing fraction**, is the maximum percentage of libraries that are allowed to NOT full-fill all above criteria; e.g. if *miss-frac=0.2* is set, the script will report sites even if 20% of all samples (e.g. 2 out of 10) do not fullfill the coverage criteria.

PoolSNP creates multiple output files:
* A gzipped **VCF file** (v.4.2) containing  allele counts and frequencies for every position and library
* A **Max-coverage file** containing the maximum coverage thressholds for all chromosomal arms and libraries in the mpileup file (separarted by a column)
* Optionally a **Bad-sites file** (by setting the parameter BS=1), which contains a list of (variable and invariable) sites that did not pass the SNP calling criteria. This file can be used to weight windows for the calulation of population genetic estimators with [PoolGEN](https://github.com/capoony/repo/PoolGen/readme.md).

PoolSNP is a shell script using GNU parallel to utilize multiple threads. Parameters need to be passed to the shell script and all necessary steps will be processed serially (or in parallel, whenever possible). The three python scripts being part of the PoolSNP pipeline require Python2.7 and can also be used as standalone scripts. Please see the documentation within each of these scripts by typing: 

```bash
python2.7 script.py -h
```

PoolSNP has been tested on Mac OSX (10.11) and Linux Ubuntu (16.10). The shellscript only works with a BASH shell and requires Python 2.7 and GNU parallel to be installed in Path

to get more help on the different parameters, execute the shellscript without parameters

### A typical command-line for an mpileup file with four samples looks like this:

```bash
bash PoolSNP.sh   \
mpileup=input.mpileup.gz \
reference=reference.fasta.gz \
names=Sample1,Sample2,Sample3,Sample4 \
max-cov=0.98 \
min-cov=10 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=22 \
BS=1 \
output=output-file
```
