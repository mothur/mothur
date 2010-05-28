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
#include "globaldata.hpp"
#include "cluster.hpp"
#include "hcluster.h"

/**********************************************************************/

class MGClusterCommand : public Command {

public:
	MGClusterCommand(string);
	~MGClusterCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadBlast* read;
	NameAssignment* nameMap;
	Cluster* cluster;
	HCluster* hcluster;
	ListVector* list;
	ListVector oldList;
	vector<seqDist> overlapMatrix;
	
	string blastfile, method, namefile, overlapFile, distFile, outputDir;
	ofstream sabundFile, rabundFile, listFile;
	double cutoff;
	float penalty;
	int precision, length, precisionLength;
	bool abort, minWanted, hclusterWanted, merge, hard;
	
	void printData(ListVector*);
	ListVector* mergeOPFs(map<string, int>, float);
	void sortHclusterFiles(string, string);
	vector<seqDist> getSeqs(ifstream&);

};

/**********************************************************************/

#endif



