#!/usr/bin/env python3
import argparse, copy, random
from scipy.stats import poisson
from math import exp, log, lgamma
from os.path import join

##########################
# Command line arguments #
##########################

parser = argparse.ArgumentParser(description='synth_read_p.py generates a sequence of read counts based on Poisson distributions for CNV states.')
parser.add_argument('-s', default='1,2,3,4', type=str, help='Poisson lambda-based states to be generated', metavar='HMMStates')
parser.add_argument('-a', default='0.5,5', type=str, help='Additional Poisson lambda-based states to generate emission probabilities for', metavar='MoreStates')
parser.add_argument('-l', default='2,4,8,16,32,64,128,256', type=str, help='Comma-separated list of Poisson lambda values to be used in synthetic read creation.', metavar='Lambdas')
parser.add_argument('-d', default=1000, type=int, help='Number of duplications to be incorporated', metavar='DuplicationCount')
parser.add_argument('-sc', default=64, type=int, help='Length of a single copy region', metavar='SingleCopyLength')
parser.add_argument('-nsc', default=64, type=int, help='Length of a non-single copy region', metavar='NonSingleCopyLength')
parser.add_argument('-o', default="Output/", type=str, help='Output directory for results', metavar='OutputDir')
args = parser.parse_args()
states = args.s
states = [float(s) for s in states.split(',')]
add_states = args.a
add_states = [float(s) for s in add_states.split(',')]
lambdas = args.l
lambdas = [int(l) for l in lambdas.split(',')]
dup_count = args.d
sc_len = args.sc
nsc_len = args.nsc
output_dir = args.o

#########################################################################################
# Generate random sequence of coverage values and DupHMM Poisson emission probabilities #
#########################################################################################

for l in lambdas:
	max_val = 0			# Used as Poisson emission probability upper boundary
	seq = []			# Stores entire sequence of hits
	rand_order = []		# Stores randomized order of CNV states to be incorporated
	rand_order_nsc = []
	for s in states:
		if s != 1:
			for i in range(0, round(dup_count / len(states))):
				rand_order_nsc.append(s)
	random.shuffle(rand_order_nsc)
	for s in rand_order_nsc:
		rand_order.append(float(1.0))
		rand_order.append(s)
	
	# Generate randomized sequences and simultaneously output locations of duplications
	cnv_pos = 1
	o_fn = join(output_dir, "dups_p" + str(l) + ".txt")
	with open(o_fn, 'w') as outfile:
		for s in rand_order:
			temp_seq = []
			rand_len = sc_len
			if s != 1:
				rand_len = nsc_len
			temp_seq = list(poisson.rvs(s*l, size=rand_len))
			temp_seq = [int(x) for x in temp_seq]	# Convert numpy.int to int
			temp_max = max(temp_seq)
			if temp_max > max_val:
				max_val = temp_max
			seq += temp_seq
			
			outline = "\t".join([str(cnv_pos), str(cnv_pos+rand_len-1), str(s)]) + "\n"
			outfile.write(outline)
			cnv_pos += rand_len
	
	# Output sequence
	o_fn = join(output_dir,"seq_p" + str(l) + ".txt")
	with open(o_fn, 'w') as outfile:
		outline = ">SEQ" + str(l) + "\n"
		outfile.write(outline)
		outline = ",".join([str(c) for c in seq])
		outfile.write(outline)
	
	# Generate and output Poisson probabilities
	p_dict = {i: [] for i in range(0,max_val)}
	emit_states = copy.copy(states)
	for s in add_states:
		emit_states.append(s)
	emit_states = sorted(emit_states)
	
	for i in range(0, max_val):
		for s in sorted(emit_states):
			prob = exp(i*log(s*l)-lgamma(i+1) - s*l)	# Poisson equation in log-space to handle cases where k is large
			p_dict[i].append(str(prob))
	
	o_fn = join(output_dir,"eprob_p" + str(l) + ".txt")
	with open(o_fn, 'w') as outfile:
		outline = "\t".join([str(s) + "X" for s in emit_states]) + "\n"
		outfile.write(outline)
		for i in range(0, max_val):
			p_list = p_dict[i]
			outline = "\t".join(p_list) + "\n"
			#outline = "\t".join(p_dict[i]) + "\n"
			outfile.write(outline)