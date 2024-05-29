import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup
import gzip

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = """python3 %prog \
      --mpileup data.mpileup \
      --min-cov 10 \
      --max-cov data.cov \
      --min-count 10 \
      --min-freq 0.01 \
      --base-quality-threshold 15 \
      --names name1,name2 \
      --coding 1.8 \
      > output.vcf"""
parser = OptionParser(usage=usage)
helptext = """

H E L P :
_________

Print a string of 0 (not passing all critera) and 1 (passing all critera), where the position refers to the position in the MPILUP file for each SNP separately
"""
group = OptionGroup(parser, helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--mpileup", dest="m", help="A mpileup file")
parser.add_option("--min-cov", dest="minc",
                  help="The minimum coverage threshold: e.g. 10", default=10)
parser.add_option("--max-cov", dest="max",
                  help="An input file with precomputed coverage thresholds")
parser.add_option("--base-quality-threshold", dest="b",
                  help="The Base-quality threshold for Qualities encoded in Sanger format (Illumina 1.8 format)", default=15)
parser.add_option("--coding", dest="c",
                  help="the Illumina FASTQ quality coding", default=1.8)

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################

def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def keywithmaxvalue(D):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash = d(list)
    for k, v in D.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]


def splitter(l, n):
    ''' This generator function returns equally sized cunks of an list'''
    # credit: Meric Lieberman, 2012
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i + n]


def extract_indel(l, sign):
    ''' This function returns an Indel from a sequence string in a pileup'''
    position = l.index(sign)
    numb = ""
    i = 0
    while True:
        if l[position + 1 + i].isdigit():
            numb += l[position + 1 + i]
            i += 1
        else:
            break

    seqlength = int(numb)
    sequence = l[position:position + i + 1 + seqlength]
    indel = sequence.replace(numb, "")

    return sequence, indel

################################## parameters ########################################


data = options.m
minimumcov = int(options.minc)
baseqthreshold = int(options.b)
phred = float(options.c)


############################ calculate PHRED cutoff  #############################

# calculate correct PHRED score cutoff: ASCII-pc

if phred >= 1.0 and phred < 1.8:
    pc = 64
else:
    pc = 33

############################ get MAX coverage threshold  #############################
maximumcov = d(list)
for l in open(options.max, "r"):
    if l.startswith("#") or l.startswith("calculating"):
        continue
    k, v = l.split("\t")
    maximumcov[k] = [int(x) for x in v.split(",")]
# print maximumcov
############################ parse MPILEUP ###########################################

# parse mpileup and store alternative alleles:

for line in load_data(data):
    if len(line.split("\t")) < 2:
        continue

    k = line[:-1].split('\t')
    CHR, POS, REF = k[:3]

    # only keep chromosomal arms with maximum coverage threshold
    if CHR not in maximumcov:
        # print CHR
        continue

    div = list(splitter(k, 3))
    libraries = div[1:]
    # loop through libraries
    totalalleles = d(int)
    alleles = d(lambda: d(int))
    covtest = []

    for j in range(len(libraries)):
        alleles[j]
        nuc = libraries[j][1]
        qualities = libraries[j][2]

        # test if seq-string is empty
        if nuc == "*":
            covtest.append(1)
            continue

        # find and remove read indices and mapping quality string
        nuc = re.sub(r'\^.', r'', nuc)
        nuc = nuc.replace('$', '')
        cov = 0

        # find and remove InDels
        while "+" in nuc or "-" in nuc:
            if "+" in nuc:
                insertion, ins = extract_indel(nuc, "+")
                nuc = nuc.replace(insertion, "")
            else:
                deletion, dele = extract_indel(nuc, "-")
                nuc = nuc.replace(deletion, "")

        # test for base quality threshold (if below: ignore nucleotide)
        # print len(nuc),len(qualities)
        nuc = "".join([nuc[x] for x in range(len(nuc)) if ord(
            qualities[x]) - pc >= baseqthreshold])
        nuc = "".join([nuc[x] for x in range(len(nuc)) if nuc[x] != "*"])

        # ignore if coverage is below or above thresholds after filtering for 1) InDels and 2) base-quality
        if len(nuc) < minimumcov or len(nuc) > maximumcov[CHR][j]:
            covtest.append(1)
        else:
            covtest.append(0)

    if len(set(covtest)) != 1:
        print(CHR + "\t" + POS + "\t" + "".join([str(x) for x in covtest]))
        continue
    if len(set(covtest)) == 1 and covtest[0] == 1:
        print(CHR + "\t" + POS + "\t" + "".join([str(x) for x in covtest]))
        continue
