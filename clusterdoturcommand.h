#ifndef CLUSTERDOTURCOMMAND_H
#define CLUSTERDOTURCOMMAND_H

/*
 *  clusterdoturcommand.h
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "nameassignment.hpp"
#include "globaldata.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"


class ClusterDoturCommand : public Command {
	
public:
	ClusterDoturCommand(string);
	ClusterDoturCommand();
	~ClusterDoturCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help();
	
private:
	bool abort, hard;
	string method, fileroot, tag, outputDir, phylipfile, namefile;
	double cutoff;
	int precision, length;
	ofstream sabundFile, rabundFile, listFile;
	NameAssignment* nameMap;
	ListVector* list;
	RAbundVector* rabund;
	RAbundVector oldRAbund;
	ListVector oldList;
	
	void printData(string label);
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
};

#endif

