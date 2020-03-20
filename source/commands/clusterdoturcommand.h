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
	~ClusterDoturCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "cluster.classic";		}
	string getCommandCategory()		{ return "Clustering";			}
    
	string getHelpString();	
    string getOutputPattern(string);		
	string getCitation() { return "Schloss PD, Westcott SL (2011). Assessing and improving methods used in OTU-based approaches for 16S rRNA gene sequence analysis. Appl Environ Microbiol 77:3219.\nSchloss PD, Handelsman J (2005). Introducing DOTUR, a computer program for defining operational taxonomic units and estimating species richness. Appl Environ Microbiol 71: 1501-6.\nhttp://www.mothur.org/wiki/Cluster.classic\n";}
	string getDescription()		{ return "cluster your sequences into OTUs using DOTUR’s method"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort, sim;
	string method, fileroot, tag, outputDir, phylipfile, namefile, countfile;
	double cutoff;
	int precision, length;
	ofstream sabundFile, rabundFile, listFile;
	NameAssignment* nameMap;
	ListVector* list;
	RAbundVector* rabund;
	RAbundVector oldRAbund;
	ListVector oldList;
	
	void printData(string label, map<string, int>&, bool&);
	vector<string> outputNames;
};

#endif

