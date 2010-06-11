#ifndef CLUSTERSPLITCOMMAND_H
#define CLUSTERSPLITCOMMAND_H

/*
 *  clustersplitcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "sparsematrix.hpp"
#include "globaldata.hpp"


class ClusterSplitCommand : public Command {
	
public:
	ClusterSplitCommand(string);
	~ClusterSplitCommand();
	int execute();	
	void help();
	
private:
	GlobalData* globaldata;
	vector<int> processIDS;   //processid
	vector<string> outputNames;

	string method, fileroot, tag, outputDir, phylipfile, columnfile, namefile, distfile, format, showabund, timing, splitmethod, taxFile, fastafile;
	double cutoff, splitcutoff;
	int precision, length, processors, taxLevelCutoff;
	bool print_start, abort, hard, large;
	time_t start;
	ofstream outList, outRabund, outSabund;
	
	void printData(ListVector*);
	int createProcesses(vector < vector < map<string, string> > >);
	vector<string> cluster(vector< map<string, string> >, set<string>&);
	int mergeLists(vector<string>, map<float, int>, ListVector*);
	map<float, int> completeListFile(vector<string>, string, set<string>&, ListVector*&);
};

#endif

