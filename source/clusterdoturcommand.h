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
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"


class ClusterDoturCommand : public Command {
	
public:
	ClusterDoturCommand(string);
	ClusterDoturCommand();
	~ClusterDoturCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "cluster.classic";		}
	string getCommandCategory()		{ return "Clustering";			}
	string getHelpString();	
	string getCitation() { return "Schloss PD, Handelsman J (2005). Introducing DOTUR, a computer program for defining operational taxonomic units and estimating species richness. Appl Environ Microbiol 71: 1501-6. \nhttp://www.mothur.org/wiki/Cluster.classic"; }
	string getDescription()		{ return "cluster your sequences into OTUs using DOTUR’s method"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort, hard, sim;
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
};

#endif

