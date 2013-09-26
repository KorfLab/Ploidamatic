#!/usr/bin/env python3
import argparse, os

parser = argparse.ArgumentParser(description='gff_to_seq.py converts a gff file containing CNV information into a SEG file that can \
											  be viewed and interpreted correctly as CNV information within IGV Viewer. This script \
											  also converts -Z and -H states into the parent state')
parser.add_argument('-d', default='ploidamatic_scratch/', type=str, help='Directory containing seq.gff', metavar='Directory')
args = parser.parse_args()
indir = args.d
i_path = os.path.join(indir, "seq.gff")
o_path_1 = os.path.join(indir, "seq.seg")
o_path_2 = os.path.join(indir, "seq.bed")

prev_chrom = ""
prev_state = ""
real_start = 0
real_end = 0
with open(i_path) as infile:
	with open(o_path_1, 'w') as outfile:
		outline = "CNV_Track\tchrom\tstart\tend\tCNV_Value\n"
		outfile.write(outline)
		with open(o_path_2, 'w') as outfile2:
			for line in infile:
				line = line.split("\t", 5)
				if line != ["\n"]:
					chrom, stochhmm, state, start, end, extra = line
					if len(state) > 2:
						if state[-1:] == "Z":
							outline = "\t".join(["CNV_Abnormal", chrom, start, end, "Low"]) + "\n"
						else:
							outline = "\t".join(["CNV_Abnormal", chrom, start, end, "High"]) + "\n"
						outfile.write(outline)
					if prev_chrom == "":
						prev_chrom = chrom
						prev_state = state[:1]
						real_start = start
					if chrom != prev_chrom:
						outline = "\t".join(["CNV_Value", prev_chrom, real_start, real_end, prev_state]) + "\n"
						outfile.write(outline)
						if prev_state != "2":
							outline2 = "\t".join([prev_chrom, str(int(real_start) - 1), str(int(real_end) - 1), prev_state]) + "\n" # BED uses 0-based position
							outfile2.write(outline2)
						prev_chrom = chrom
						prev_state = state[:1]
						real_start = start
					if state[:1] != prev_state:
						outline = "\t".join(["CNV_Value", chrom, real_start, real_end, prev_state]) + "\n"
						outfile.write(outline)
						if prev_state != "2":
							outline2 = "\t".join([prev_chrom, str(int(real_start) - 1), str(int(real_end) - 1), prev_state]) + "\n" # BED uses 0-based position
							outfile2.write(outline2)
						prev_state = state[:1]
						real_start = start
					real_end = end