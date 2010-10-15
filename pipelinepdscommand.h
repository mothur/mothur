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
	PipelineCommand() {}
	~PipelineCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; } //empty
	int execute();
	void help();
	
private:
	bool abort;
	CommandFactory* cFactory;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	vector<string> commands;
	string outputDir, sffFile, alignFile, oligosFile, taxonomyFile, pipeFilename, classifyFile, chimeraFile;
	int processors;
	
	bool readUsersPipeline();
	int runUsersPipeline();
	void createPatsPipeline();
	bool parseCommand(string, string&, string&);
	bool checkForValidAndRequiredParameters(string, string, map<string, vector<string> >&);
	bool fillInMothurMade(string&, map<string, vector<string> >);
};

/****************************************************/

#endif

