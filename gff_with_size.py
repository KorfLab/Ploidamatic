#!/usr/bin/env python3
import argparse, os

parser = argparse.ArgumentParser(description='gff_with_size.py will add a column containing the size of a CNV into a GFF file')
parser.add_argument('-d', default='ploidamatic_scratch', type=str, help='', metavar='Directory')
args = parser.parse_args()
indir = args.d
i_path = os.path.join(indir, "seq.gff")
o_path = os.path.join(indir, "seq_size.gff")

with open(i_path) as infile:
    with open(o_path, 'w') as outfile:
        for line in infile:
            line = line.strip()
            line = line.split()
            size = int(line[4]) - int(line[3])
            line = "\t".join(line) + "\t" + str(size) + "\n"
            outfile.write(line)