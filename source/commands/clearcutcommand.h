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

/* 
  Evans, J., L. Sheneman, and J.A. Foster (2006) Relaxed Neighbor-Joining: 
  A Fast Distance-Based Phylogenetic Tree Construction Method, 
  J. Mol. Evol., 62, 785-792
 */ 

/****************************************************************************/

class ClearcutCommand : public Command {

public:
	ClearcutCommand(string);
	~ClearcutCommand() = default;
	
	vector<string> setParameters();
	string getCommandName()			{ return "clearcut";			}
	string getCommandCategory()		{ return "Hypothesis Testing";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Sheneman L, Evans J, Foster JA (2006). Clearcut: a fast implementation of relaxed neighbor joining. Bioinformatics 22: 2823-4. \nhttp://www.mothur.org/wiki/Clearcut"; }
	string getDescription()		{ return "create a tree from a fasta or phylip file"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string  phylipfile, fastafile, matrixout, inputFile, seed, ntrees;
	bool version, verbose, quiet, norandom, shuffle, neighbor, expblen, expdist, stdoutWanted, kimura, jukes, protein, DNA;
	bool abort;
	vector<string> outputNames;
	
};

/****************************************************************************/

#endif

