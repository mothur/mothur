#ifndef MGCLUSTERCOMMAND_H
#define MGCLUSTERCOMMAND_H

/*
 *  mgclustercommand.h
 *  Mothur
 *
 *  Created by westcott on 12/11/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "readblast.h"
#include "nameassignment.hpp"
#include "cluster.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "counttable.h"
#include "opticluster.h"
#include "optimatrix.h"

/**********************************************************************/

class MGClusterCommand : public Command {

public:
	MGClusterCommand(string);
	MGClusterCommand();
	~MGClusterCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "mgcluster";	}
	string getCommandCategory()		{ return "Clustering";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Handelsman J (2008). A statistical toolbox for metagenomics. BMC Bioinformatics 9: 34. \nhttp://www.mothur.org/wiki/Mgcluster"; }
	string getDescription()		{ return "cluster your sequences into OTUs using a blast file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	ReadBlast* read;
	NameAssignment* nameMap;
	Cluster* cluster;
	ListVector* list;
    ListVector oldList;
    CountTable* ct;
    RAbundVector rav;
	vector<seqDist> overlapMatrix;
	vector<string> outputNames;
	
	string blastfile, method, namefile, countfile, overlapFile, distFile, outputDir, sabundFileName, rabundFileName, listFileName, metric, initialize, tag, fileroot;
	ofstream sabundFile, rabundFile, listFile;
	double cutoff;
	float penalty, adjust, stableMetric;
	int precision, length, precisionLength, maxIters;
	bool abort, minWanted, hclusterWanted, merge, cutoffSet;
	
	void printData(ListVector*, map<string, int>&, bool&);
	ListVector* mergeOPFs(map<string, int>, float);
    vector<seqDist> getSeqs(ifstream&);
    void createRabund(CountTable*&, ListVector*&, RAbundVector*&);
    int runOptiCluster();
    int runMothurCluster();
};

/**********************************************************************/

#endif



