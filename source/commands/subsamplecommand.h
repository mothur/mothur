#ifndef SUBSAMPLECOMMAND_H
#define SUBSAMPLECOMMAND_H

/*
 *  subsamplecommand.h
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"

#include "listvector.hpp"
#include "rabundvector.hpp"
#include "inputdata.h"
#include "sequence.hpp"
#include "counttable.h"


class SubSampleCommand : public Command {

public:
	SubSampleCommand(string);
	SubSampleCommand();
	~SubSampleCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sub.sample";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Sub.sample"; }
	string getDescription()		{ return "get a sampling of sequences from a list, shared, rabund, sabund or fasta file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:	
	bool abort, pickedGroups, allLines, persample, withReplacement;
	string listfile, groupfile, countfile, sharedfile, rabundfile, sabundfile, fastafile, namefile, taxonomyfile;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir;
	vector<string> Groups, outputNames;
	int size;
	vector<string> names;
	map<string, vector<string> > nameMap;
    CountTable ct;
	
	int getSubSampleShared();
	int getSubSampleList();
	int getSubSampleRabund();
	int getSubSampleSabund();
	int getSubSampleFasta();
	int processShared(SharedRAbundVectors*&, bool&);
	int processRabund(RAbundVector*&, ofstream&);
	int processSabund(SAbundVector*&, ofstream&);
	int processList(ListVector*&, set<string>&);
	int getNames();
	int readNames();
    int getTax(set<string>&);
	
};

#endif

