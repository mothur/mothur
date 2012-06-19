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

//**********************************************************************************************************************
class GetSharedOTUCommand : public Command {
	
	public:
	
		GetSharedOTUCommand(string);	
		GetSharedOTUCommand();	
		~GetSharedOTUCommand() {}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.sharedseqs";			}
		string getCommandCategory()		{ return "OTU-Based Approaches";	}
		string getRequiredCommand()		{ return "none";					}
		string getOutputFileNameTag(string, string);
	string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Get.sharedseqs"; }
		string getDescription()		{ return "identifies sequences that are either unique or shared by specific groups"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
		
	private:
		ListVector* list;
		GroupMap* groupMap;
		
		set<string> labels;
		string fastafile, label, groups, listfile, groupfile, output, userGroups, outputDir, format;
		bool abort, allLines, unique;
		vector<string> Groups;
		map<string, string> groupFinder;
		map<string, string>::iterator it;
		vector<Sequence> seqs;
		vector<string> outputNames;
		
		int process(ListVector*);
		
};
//**********************************************************************************************************************

#endif

