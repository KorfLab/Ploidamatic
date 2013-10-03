#!/usr/bin/env python3
import argparse, os, re

parser = argparse.ArgumentParser(description='genes_within_cnvs.py takes in a list of CNV loci and a GFF file and outputs \
                                 all non-transposable element genes located within those loci')
parser.add_argument('-d', default="SuecicaProject/Sue1Set5/", type=str, help='Directory where duplicated_loci.txt is located \
                    and where CNV genes will be output', metavar='Directory')
parser.add_argument('-g', default="../Genomes/Thalyrata.gff", type=str, help='Path to GFF file', metavar='GFF')
args = parser.parse_args()
indir = args.d
gff_file = args.g

i_file = os.path.join(indir, "cnv_loci.txt")
cnv_loci = {}
with open(i_file) as infile:
    for line in infile:
        line = line.split()
        chrom, start, end, state = line
        start = int(start)
        end = int(end)
        if chrom not in cnv_loci:
            cnv_loci[chrom] = []
        cnv_loci[chrom].append((start, end, state))

pattern = r'Name=(AT\d{1}G\d{5})'
recomp = re.compile(pattern)
chrom_accept = ['AtChr1','AtChr2','AtChr3','AtChr4','AtChr5']
results = {}
with open(gff_file) as infile:
    for line in infile:
        line = line.split()
        chrom = line[0]
        feature = line[2]
        if feature == "gene" and chrom in chrom_accept:
            start = int(line[3])
            end = int(line[4])
            desc = line[8]
            match = recomp.search(desc)
            gene = match.group(1)
            for (ref_start, ref_end, state) in cnv_loci[chrom]:
                if (start > ref_start and start < ref_end) or (end > ref_start and end < ref_end):
                    match = recomp.search(desc)
                    gene = match.group(1)
                    if state not in results:
                        results[state] = []
                    results[state].append(gene)
                    break

for state in results:
    o_file = os.path.join(indir, state + "_genes.txt")
    with open(o_file, 'w') as outfile:
        for gene in results[state]:
            outline = gene + "\n"
            outfile.write(outline)