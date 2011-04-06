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
	
	int driver(linePair*, string, string, string);
	int createProcesses(string, string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
	#endif

	bool abort, filter;
	string fastafile, templatefile, outputDir, maskfile;
	int processors, window, numwanted, numSeqs, templateSeqsLength;
	Chimera* chimera;
	vector<string> fastaFileNames;
	vector<string> outputNames;
};

/***********************************************************/

#endif

