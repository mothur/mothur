#ifndef CHIMERACCODECOMMAND_H
#define CHIMERACCODECOMMAND_H

/*
 *  chimeraccodecommand.h
 *  Mothur
 *
 *  Created by westcott on 3/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "chimera.h"


/***********************************************************/

class ChimeraCcodeCommand : public Command {
public:
	ChimeraCcodeCommand(string);
	ChimeraCcodeCommand();
	~ChimeraCcodeCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.ccode";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Gonzalez JM, Zimmermann J, Saiz-Jimenez C (2005). Evaluating putative chimeric sequences from PCR-amplified products. Bioinformatics 21: 333-7. \nhttp://www.mothur.org/wiki/Chimera.ccode"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
		
private:
	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	
	int driver(linePair*, string, string, string);
	int createProcesses(string, string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long long>&);
	#endif

	bool abort, filter, save;
	string fastafile, templatefile, outputDir, maskfile;
	int processors, window, numwanted, numSeqs, templateSeqsLength;
	Chimera* chimera;
	vector<string> fastaFileNames;
	vector<string> outputNames;
};

/***********************************************************/

#endif

