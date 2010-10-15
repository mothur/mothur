#ifndef GETSHAREDOTUCOMMAND_H
#define GETSHAREDOTUCOMMAND_H

/*
 *  getsharedotucommand.h
 *  Mothur
 *
 *  Created by westcott on 9/22/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "listvector.hpp"
#include "sequence.hpp"
#include "groupmap.h"
#include "globaldata.hpp"

//**********************************************************************************************************************
class GetSharedOTUCommand : public Command {
	
	public:
	
		GetSharedOTUCommand(string);	
		GetSharedOTUCommand();	
		~GetSharedOTUCommand();
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();
		void help();	
		
	private:
		
		GlobalData* globaldata;
		ListVector* list;
		GroupMap* groupMap;
		
		set<string> labels;
		string fastafile, label, groups, listfile, groupfile, output, userGroups, outputDir;
		bool abort, allLines, unique;
		vector<string> Groups;
		map<string, string> groupFinder;
		map<string, string>::iterator it;
		vector<Sequence> seqs;
		vector<string> outputNames;
		map<string, vector<string> > outputTypes;
		
		int process(ListVector*);
		
};
//**********************************************************************************************************************

#endif

