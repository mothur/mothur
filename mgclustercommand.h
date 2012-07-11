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
#include "sparsematrix.hpp"
#include "nameassignment.hpp"
#include "cluster.hpp"
#include "hcluster.h"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "counttable.h"

/**********************************************************************/

class MGClusterCommand : public Command {

public:
	MGClusterCommand(string);
	MGClusterCommand();
	~MGClusterCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "mgcluster";	}
	string getCommandCategory()		{ return "Clustering";	}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "Schloss PD, Handelsman J (2008). A statistical toolbox for metagenomics. BMC Bioinformatics 9: 34. \nhttp://www.mothur.org/wiki/Mgcluster"; }
	string getDescription()		{ return "cluster your sequences into OTUs using a blast file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	ReadBlast* read;
	NameAssignment* nameMap;
	Cluster* cluster;
	HCluster* hcluster;
	ListVector* list;
    CountTable* ct;
	ListVector oldList;
    RAbundVector rav;
	vector<seqDist> overlapMatrix;
	vector<string> outputNames;
	
	string blastfile, method, namefile, countfile, overlapFile, distFile, outputDir;
	ofstream sabundFile, rabundFile, listFile;
	double cutoff;
	float penalty;
	int precision, length, precisionLength;
	bool abort, minWanted, hclusterWanted, merge, hard, large;
	
	void printData(ListVector*);
	ListVector* mergeOPFs(map<string, int>, float);
	void sortHclusterFiles(string, string);
	vector<seqDist> getSeqs(ifstream&);
    void createRabund(CountTable*&, ListVector*&, RAbundVector*&);

};

/**********************************************************************/

#endif



