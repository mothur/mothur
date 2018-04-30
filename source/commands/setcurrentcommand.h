#ifndef SETCURRENTCOMMAND_H
#define SETCURRENTCOMMAND_H

/*
 *  setcurrentcommand.h
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"

class SetCurrentCommand : public Command {
	
public:
	SetCurrentCommand(string);
	SetCurrentCommand();
	~SetCurrentCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "set.current";	}
	string getCommandCategory()		{ return "General";		}
	
	string getHelpString();	
    string getOutputPattern(string);
	string getCitation() { return "http://www.mothur.org/wiki/Set.current"; }
	string getDescription()		{ return "set current files for mothur"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	
	vector<string> outputNames;
	bool abort;
	
	string clearTypes;
	vector<string> types;
	
	string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile, countfile, summaryfile, currentFile, outputDir;
	string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile, filefile, constaxonomyfile, contigsreportfile;

    int readCurrentFiles();
	
};

#endif


