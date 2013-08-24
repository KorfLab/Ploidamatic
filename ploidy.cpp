/*****************************************************************************\
 ploidy.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/

#include <StochHMMlib.h>

#define MAX_COUNT 1000

static std::vector<std::vector<double> > EmitTable;
static size_t Last_row(0);

static void initialize_emission_vectors(std::string& filename){
	
	std::ifstream efh;
	efh.open(filename.c_str());
	if (!efh.is_open()){
		std::cerr << "Could not open Emission Probability File\n";
	}
	
	StochHMM::stringList header_row;
	std::string header;
	getline(efh,header,'\n');
	header_row.splitString(header, "\t");
	EmitTable.resize(header_row.size(),
		std::vector<double>(MAX_COUNT,-INFINITY));
	
	std::string line;
	size_t n(0);
	while (getline(efh, line, '\n')){
		StochHMM::stringList row;
		row.splitString(line, "\t");
		
		for (size_t model = 0; model < row.size(); model++){
			double value;
			if (!StochHMM::stringToDouble(row[model], value)) value = 0;		 
			EmitTable[model][n] = log(value);
		}
		n++;
	}
	
	for(size_t i = 0 ;i < EmitTable.size(); i++) {
		EmitTable[i].resize(n);
	}
	
	Last_row = n - 1;
	
}

static double emission (const double val, const std::vector<double>* param) {
	size_t model((*param)[0]);
	size_t value(val);
	
	if (value > Last_row) return EmitTable[model][Last_row];
	else                  return EmitTable[model][value];
}

int main (int argc, const char * argv[]) {

	/* command line */
	std::string usage = "usage: ploidy <hmm file> <emission file> <seq file>";
	if (argc != 4) {
		std::cout << usage << std::endl;
		exit(2);
	}
	std::string hmm_file  =	  argv[1];
	std::string emissions =	  argv[2];
	std::string seq_file  =	  argv[3];
	
	/* init */
	initialize_emission_vectors(emissions);
	StochHMM::StateFuncs hmm_functions;
	hmm_functions.assignPDFFunction("P05x", *emission);
	hmm_functions.assignPDFFunction("P1x",	*emission);
	hmm_functions.assignPDFFunction("P2x",	*emission);
	hmm_functions.assignPDFFunction("P3x",	*emission);
	hmm_functions.assignPDFFunction("P4x",	*emission);
	hmm_functions.assignPDFFunction("P5x",	*emission);
	StochHMM::model hmm;
	hmm.import(hmm_file, &hmm_functions);
	
	/* decode & output */
	StochHMM::seqTracks jobs;
	jobs.loadSeqs(hmm, seq_file);
	StochHMM::seqJob *job=jobs.getJob();
	while (job != NULL) {
		StochHMM::trellis trell(&hmm, job->getSeqs());
		trell.simple_viterbi();
		StochHMM::traceback_path tb(&hmm);
		trell.traceback(tb);
		tb.print_gff(job->getHeader());
		job = jobs.getJob();
	}
	
	return 0;
}
