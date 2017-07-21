#ifndef OTUASSOCIATIONCOMMAND_H
#define OTUASSOCIATIONCOMMAND_H

/*
 *  otuassociationcommand.h
 *  Mothur
 *
 *  Created by westcott on 1/19/12.
 *  Copyright 2012 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"


class OTUAssociationCommand : public Command {
public:
	OTUAssociationCommand(string);
	OTUAssociationCommand();
	~OTUAssociationCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "otu.association";			}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Otu.association"; }
	string getDescription()		{ return "calculate the correlation coefficient for the otus in a shared/relabund file"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }	
private:
	
	string sharedfile, relabundfile, metadatafile, groups, label, inputFileName, outputDir, method;
	bool abort, pickedGroups, allLines;
    double cutoff;
	set<string> labels;
    vector< vector< double> > metadata;
	
	vector<string> outputNames, Groups, metadataLabels;
	int processShared();
	int process(SharedRAbundVectors*);
	int processRelabund();
	int process(SharedRAbundFloatVectors*);
    int readMetadata();
    //int getMetadata();
	
};


#endif



