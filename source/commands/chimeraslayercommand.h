#ifndef CHIMERASLAYERCOMMAND_H
#define CHIMERASLAYERCOMMAND_H

/*
 *  chimeraslayercommand.h
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "mothurchimera.h"
#include "chimeraslayer.h"
#include "sequenceparser.h"
#include "sequencecountparser.h"

/***********************************************************/

class ChimeraSlayerCommand : public Command {
public:
	ChimeraSlayerCommand(string);
	ChimeraSlayerCommand();
	~ChimeraSlayerCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.slayer";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Haas BJ, Gevers D, Earl A, Feldgarden M, Ward DV, Giannokous G, Ciulla D, Tabbaa D, Highlander SK, Sodergren E, Methe B, Desantis TZ, Petrosino JF, Knight R, Birren BW (2011). Chimeric 16S rRNA sequence formation and detection in Sanger and 454-pyrosequenced PCR amplicons. Genome Res  21:494.\nhttp://www.mothur.org/wiki/Chimera.slayer"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:

	vector<int> processIDS;   //processid
	vector<linePair> lines;
	
	int driver(string, string, string, string, map<string, int>&);
	int divideInHalf(Sequence, string&, string&);
	map<string, int> sortFastaFile(string, string);
	map<string, int> sortFastaFile(vector<Sequence>&, map<string, string>&, string newFile);
    int sortFastaFile(vector<Sequence>&, map<string, int>&, string newFile);
	string getNamesFile(string&);
	//int setupChimera(string,);
    int deconvoluteResults(map<string, string>&, string, string, string);
	map<string, int> priority;
	int setUpForSelfReference(SequenceParser*&, map<string, string>&, map<string, map<string, int> >&, int);
    int setUpForSelfReference(SequenceCountParser*&, map<string, string>&, map<string, map<string, int> >&, int);
	int driverGroups(string, string, string, map<string, map<string, int> >&, map<string, string>&, string);

	bool abort, realign, trim, trimera, save, hasName, hasCount, dups;
	string fastafile, groupfile, templatefile, outputDir, search, namefile, countfile, blastlocation;
	int window, iters, increment, numwanted, ksize, match, mismatch, parents, minSimilarity, minCoverage, minBS, minSNP, templateSeqsLength;
    long long numSeqs;
	float divR;
	
    map<string, map<string, string> > group2NameMap;
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	vector<string> groupFileNames;
	
};

/***********************************************************/

#endif


