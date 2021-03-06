#!/usr/bin/env python3
import argparse, os

parser = argparse.ArgumentParser(description='gff_to_bed.py converts a gff file containing CNV information into a BED file that can \
											  be viewed and interpreted correctly as CNV information within IGV Viewer. This script \
											  also converts -Z and -H states into the parent state')
parser.add_argument('-d', default='SuecicaProject/Sue1Set5/', type=str, help='Directory containing seq.gff', metavar='Directory')
args = parser.parse_args()
indir = args.d
i_path = os.path.join(indir, "seq.gff")
o_path = os.path.join(indir, "seq.bed")

prev_chrom = ""
prev_state = ""
real_start = 0
real_end = 0
with open(o_path, 'w') as outfile:
	with open(i_path) as infile:
		for line in infile:
			line = line.split("\t", 5)
			chrom, stochhmm, state, start, end, extra = line
			if prev_chrom == "":
				prev_chrom = chrom
				prev_state = state[:1]
				real_start = start
			if chrom != prev_chrom:
				if prev_state != "2":
					outline = "\t".join([prev_chrom, str(int(real_start) - 1), str(int(real_end) - 1), prev_state]) + "\n" # BED uses 0-based position
					outfile.write(outline)
				prev_chrom = chrom
				prev_state = state[:1]
				real_start = start
			if state[:1] != prev_state:
				if prev_state != "2":
					outline = "\t".join([prev_chrom, str(int(real_start) - 1), str(int(real_end) - 1), prev_state]) + "\n" # BED uses 0-based position
					outfile.write(outline)
				prev_state = state[:1]
				real_start = start
			real_end = end
	if prev_state != "2":
		outline = "\t".join([prev_chrom, str(int(real_start) - 1), str(int(real_end) - 1), prev_state]) + "\n" # BED uses 0-based position
		outfile.write(outline)