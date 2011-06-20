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
	ClassifySeqsCommand();
	~ClassifySeqsCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "classify.seqs";		}
	string getCommandCategory()		{ return "Phylotype Analysis";	}
	string getHelpString();	
	string getCitation() { return "Wang Q, Garrity GM, Tiedje JM, Cole JR (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol 73: 5261-7. [ for Bayesian classifier ] \nAltschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25: 3389-402. [ for BLAST ] \nDeSantis TZ, Hugenholtz P, Larsen N, Rojas M, Brodie EL, Keller K, Huber T, Dalevi D, Hu P, Andersen GL (2006). Greengenes, a chimera-checked 16S rRNA gene database and workbench compatible with ARB. Appl Environ Microbiol 72: 5069-72. [ for kmer ] \nhttp://www.mothur.org/wiki/Classify.seqs"; }
	string getDescription()		{ return "classify sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
	
private:
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};

	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<string> fastaFileNames;
	vector<string> namefileNames;
	vector<string> groupfileNames;
	vector<string> outputNames;
	map<string, vector<string> > nameMap;
	map<string,  vector<string> >::iterator itNames;
	
	Classify* classify;
	
	string fastaFileName, templateFileName, distanceFileName, namefile, search, method, taxonomyFileName, outputDir, groupfile;
	int processors, kmerSize, numWanted, cutoff, iters;
	float match, misMatch, gapOpen, gapExtend;
	bool abort, probs;
	
	int driver(linePair*, string, string, string);
	void appendTaxFiles(string, string);
	int createProcesses(string, string, string); 
	string addUnclassifieds(string, int);
	
	int MPIReadNamesFile(string);
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
	#endif
};

#endif

