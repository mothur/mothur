#ifndef PIPELINEPDSCOMMAND_H
#define PIPELINEPDSCOMMAND_H

/*
 *  pipelinepdscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/5/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "commandfactory.hpp"

/****************************************************/

class PipelineCommand : public Command {
	
public:
	PipelineCommand(string);
	PipelineCommand() { abort = true; calledHelp = true; setParameters(); }
	~PipelineCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pipeline.pds";	}
	string getCommandCategory()		{ return "Hidden";			}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Pipeline.pds"; }
	string getDescription()		{ return "pat's pipeline"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;
	CommandFactory* cFactory;
	vector<string> outputNames;
	vector<string> commands;
	string outputDir, sffFile, alignFile, oligosFile, taxonomyFile, pipeFilename, classifyFile, chimeraFile;
	int processors;
	
	bool readUsersPipeline();
	int runUsersPipeline();
	void createPatsPipeline();
	bool parseCommand(string, string&, string&);
	bool checkForValidAndRequiredParameters(string, string, map<string, vector<string> >&);
	bool fillInMothurMade(string&, map<string, vector<string> >&);
};

/****************************************************/

#endif

