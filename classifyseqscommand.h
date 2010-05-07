#ifndef CLASSIFYSEQSCOMMAND_H
#define CLASSIFYSEQSCOMMAND_H

/*
 *  classifyseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/2/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
//#include "alignment.hpp"
#include "classify.h"

//KNN and Bayesian methods modeled from algorithms in
//Naı¨ve Bayesian Classiﬁer for Rapid Assignment of rRNA Sequences 
//into the New Bacterial Taxonomy􏰎† 
//Qiong Wang,1 George M. Garrity,1,2 James M. Tiedje,1,2 and James R. Cole1* 
//Center for Microbial Ecology1 and Department of Microbiology and Molecular Genetics,2 Michigan State University, 
//East Lansing, Michigan 48824 
//Received 10 January 2007/Accepted 18 June 2007 



class ClassifySeqsCommand : public Command {
	
public:
	ClassifySeqsCommand(string);
	~ClassifySeqsCommand();
	int execute(); 
	void help();	
	
private:
	struct linePair {
		int start;
		int numSeqs;
		linePair(int i, int j) : start(i), numSeqs(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<string> fastaFileNames;
	vector<string> namefileNames;
	vector<string> groupfileNames;
	map<string, int> nameMap;
	map<string, int>::iterator itNames;
	
	Classify* classify;
	
	string fastaFileName, templateFileName, distanceFileName, namefile, search, method, taxonomyFileName, outputDir, groupfile;
	int processors, kmerSize, numWanted, cutoff, iters;
	float match, misMatch, gapOpen, gapExtend;
	bool abort, probs;
	
	int driver(linePair*, string, string, string);
	void appendTaxFiles(string, string);
	void createProcesses(string, string, string); 
	string addUnclassifieds(string, int);
	
	int MPIReadNamesFile(string);
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<long>&);
	#endif
};

#endif

