#ifndef ALIGNCOMMAND_H
#define ALIGNCOMMAND_H

/*
 *  aligncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "alignmentdb.h"

class AlignCommand : public Command {
	
public:
	AlignCommand(string);	
	AlignCommand();
	~AlignCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "align.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "DeSantis TZ, Jr., Hugenholtz P, Keller K, Brodie EL, Larsen N, Piceno YM, Phan R, Andersen GL (2006). NAST: a multiple sequence alignment server for comparative analysis of 16S rRNA genes. Nucleic Acids Res 34: W394-9.\nSchloss PD (2009). A high-throughput DNA sequence aligner for microbial ecology studies. PLoS ONE 4: e8230.\nSchloss PD (2010). The effects of alignment quality, distance calculation method, sequence filtering, and region on the analysis of 16S rRNA gene-based studies. PLoS Comput Biol 6: e1000844.\nhttp://www.mothur.org/wiki/Align.seqs http://www.mothur.org/wiki/Align.seqs"; }
	string getDescription()		{ return "align sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	bool MPIWroteAccnos;
	
	AlignmentDB* templateDB;
	Alignment* alignment;
	
	int driver(linePair*, string, string, string, string);
	int createProcesses(string, string, string, string);
	void appendAlignFiles(string, string); 
	void appendReportFiles(string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
	#endif
	
	string candidateFileName, templateFileName, distanceFileName, search, align, outputDir;
	float match, misMatch, gapOpen, gapExtend, threshold;
	int processors, kmerSize;
	vector<string> candidateFileNames;
	vector<string> outputNames;
	
	bool abort, flip, calledHelp;

};

#endif
