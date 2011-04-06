#ifndef HCLUSTERCOMMAND_H
#define HCLUSTERCOMMAND_H

/*
 *  hclustercommand.h
 *  Mothur
 *
 *  Created by westcott on 10/13/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "hcluster.h"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "readcluster.h"

/******************************************************************/
//This command is an implementation of the HCluster algorythmn described in 
//ESPRIT: estimating species richness using large collections of 16S rRNA pyrosequences by
//Yijun Sun1,2,*, Yunpeng Cai2, Li Liu1, Fahong Yu1, Michael L. Farrell3, William McKendree3 
//and William Farmerie1 1 

//Interdisciplinary Center for Biotechnology Research, 2Department of Electrical and Computer Engineering, 
//University of Florida, Gainesville, FL 32610-3622 and 3Materials Technology Directorate, Air Force Technical 
//Applications Center, 1030 S. Highway A1A, Patrick AFB, FL 32925-3002, USA 
//Received January 28, 2009; Revised April 14, 2009; Accepted April 15, 2009 
/************************************************************/
class HClusterCommand : public Command {
	
public:
	HClusterCommand(string);
	HClusterCommand();	
	~HClusterCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "hcluster";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	HCluster* cluster;
	ListVector* list;
	RAbundVector* rabund;
	RAbundVector oldRAbund;
	ListVector oldList;
	ReadCluster* read;
	
	bool abort, sorted, print_start, hard;
	string method, fileroot, tag, distfile, format, phylipfile, columnfile, namefile, sort, showabund, timing, outputDir;
	double cutoff;
	int precision, length;
	ofstream sabundFile, rabundFile, listFile;
	time_t start;
	unsigned long loops;
	vector<string> outputNames;
	
	void printData(string label);
};

/************************************************************/

#endif
