#!/usr/bin/env python3
import argparse, subprocess
from os.path import join

##########################
# Command line arguments #
##########################

parser = argparse.ArgumentParser(description='synth_read_c.py generates a sequence of read counts based on real data histograms for CNV states.')
parser.add_argument('-s', default='1,2,3,4', type=str, help='Poisson lambda-based states to be generated', metavar='HMMStates')
parser.add_argument('-a', default='0.5,5', type=str, help='Additional Poisson lambda-based states to generate emission probabilities for', metavar='MoreStates')
parser.add_argument('-d', default=1000, type=int, help='Number of duplications to be incorporated', metavar='DuplicationCount')
parser.add_argument('-sc', default=64, type=int, help='Length of a single copy region', metavar='SingleCopyLength')
parser.add_argument('-nsc', default=64, type=int, help='Length of a non-single copy region', metavar='NonSingleCopyLength')
args = parser.parse_args()
states = args.s
states = [float(s) for s in states.split(',')]
add_states = args.a
add_states = [float(s) for s in add_states.split(',')]
dup_count = args.d
sc_len = args.sc
nsc_len = args.nsc

o_fn = "HMMInputfiles.txt"
i_fn_1 = os.path.join(outdir,str(win) + "bp_window/" + str(chromosome) + "_hist_distfits.txt")
if not os.path.exists(i_fn_1):
	i_fn_2 = os.path.join(outdir, str(win) + "bp_window/" + str(chromosome) + "_hist.txt")
	line = str(i_fn_2) + "\n"
	with open(o_fn, 'w') as outfile:
		outfile.write(line)
	# Now call R script "dist_fit.R" to perform Poisson regression
	Rpath = "Rscript"
	params = ' '.join([str(Rpath) + " dist_fit.R", "--no-save"])
	simulation = subprocess.Popen(params, shell=True)
	simulation.wait()
	os.remove(o_fn)	# Remove temporary file "HMMInputfiles.txt"
	# Read in "dist_fit.R" results file and grab Poisson regression lambda value
	with open(i_fn_1) as infile:
		inlines = infile.readlines()
		inlines = [line.split() for line in inlines]
		singlecopy_lambda[win][chromosome] = float(inlines[1][1])
		# Extra line