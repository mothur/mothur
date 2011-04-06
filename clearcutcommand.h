#ifndef CLEARCUTCOMMAND_H
#define CLEARCUTCOMMAND_H

/*
 *  clearcutcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/11/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "globaldata.hpp"

/* 
  Evans, J., L. Sheneman, and J.A. Foster (2006) Relaxed Neighbor-Joining: 
  A Fast Distance-Based Phylogenetic Tree Construction Method, 
  J. Mol. Evol., 62, 785-792
 */ 

/****************************************************************************/

class ClearcutCommand : public Command {

public:
	ClearcutCommand(string);
	ClearcutCommand();
	~ClearcutCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "clearcut";			}
	string getCommandCategory()		{ return "Hypothesis Testing";	}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string outputDir, phylipfile, fastafile, matrixout, inputFile, seed, ntrees;
	bool version, verbose, quiet, norandom, shuffle, neighbor, expblen, expdist, stdoutWanted, kimura, jukes, protein, DNA;
	bool abort;
	vector<string> outputNames;
	
};

/****************************************************************************/

#endif

