#!/usr/bin/env python3
import argparse
from collections import Counter
from os.path import join

##########################
# Command line arguments #
##########################

parser = argparse.ArgumentParser(description='evaluator_p.py is a program which will take in the path from a StochHMM run on \
								 randomly generated Poisson-based windowed read coverage and compare StochHMM\'s results to \
								 the known locations of duplications.')
parser.add_argument('-p', default='Output/16_win_len/', type=str, help='Directory containing the read simulation files and StochHMM output ("dups", "path")', metavar='CNVs')
parser.add_argument('-l', default=16, type=int, help='Lambda value whose results are being examined', metavar='Lambda')
args = parser.parse_args()
path = args.p
mean = args.l
dup_fn = join(path, "dups_p" + str(mean) + ".txt")
path_fn = dup_fn = join(path, "path_p" + str(mean) + ".txt")

#######################################
# Store window location of CNV states #
#######################################

states = []
dup_win = []
with open(dup_fn) as infile:
	for line in infile:
		line = line.strip()
		line = line.split("\t")
		for i in range(line[0],line[1]+1):
			dup_win.append(line[2])
			states.append(line[2])
states = sorted(list(set(states)))

#################################################
# Store path from StochHMM and DupHMM's results #
#################################################

cnv_results = {s: Counter() for s in states}
with open(path_fn) as infile:
	for line in infile:
		line = line.strip()
		line = line.split("\t")
		tmp_state = line[2][1:]
		s_pos = int(line[3])
		e_pos = int(line[4])
		for i in range(s_pos-1, e_pos):	# Path starts from 1, not 0
			cnv_results[dup_win[i]][path[i]] += 1

##########################################
# Output effectiveness of DupHMM results #
##########################################

o_fn = join(path, "results_p" + str(mean) + ".txt")
with open(o_fn, 'w') as outfile:
	for real_s in sorted(cnv_results):
		first_other_s = sorted(cnv_results[real_s])[0]
		outline = real_s + "X\t" + first_other_s + "X: " + str(cnv_results[real_s][first_other_s]) + "(" + str(round(cnv_results[real_s][first_other_s] / sum(cnv_results[real_s]) * 100,2)) + "%)\n"
		outfile.write(outline)
		for other_s in sorted(cnv_results[real_s]):
			outline = "\t" + other_s + "X: " + str(cnv_results[real_s][other_s]) + "(" + str(round(cnv_results[real_s][other_s] / sum(cnv_results[real_s]) * 100,2)) + "%)\n"
			outfile.write(outline)