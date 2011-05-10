#ifndef READDISTCOMMAND_H
#define READDISTCOMMAND_H
/*
 *  readdistcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "readmatrix.hpp"
#include "fullmatrix.h"
#include "groupmap.h"

/* The read.dist command is used to read a distance matrix file.  
The read.dist command parameter options are phylipfile, columnfile, namefile, cutoff and precision. 
The read.dist command should be in the following format: read.dist(phylipfile=yourDistFile, 
namefile=yourNameFile, cutoff=yourCutoff, precision=yourPrecision). The phylipfile or columnfile are required and if you use a columnfile the namefile is required.  
If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.  */

class NameAssignment;

class ReadDistCommand : public Command {
public:
	ReadDistCommand(string);
	ReadDistCommand() { abort = true; calledHelp = true; }
	~ReadDistCommand() {}
	
	vector<string> setParameters() {  return outputNames; } //dummy doesn't really do anything
	string getCommandName()			{ return "read.dist";	}
	string getCommandCategory()		{ return "Hidden";	}
	string getHelpString() { return "This command is no longer available. You can provide your distance files directly to the downstream commands like cluster."; }	
	string getCitation() { return "http://www.mothur.org/wiki/Read.dist"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	double cutoff;
	int precision;
	ReadMatrix* read;
	FullMatrix* matrix;
	GroupMap* groupMap;
	string distFileName, format, method;
	string phylipfile, columnfile, namefile, groupfile, outputDir;
	NameAssignment* nameMap;
	vector<string> outputNames;

	bool abort, sim;

};

#endif
