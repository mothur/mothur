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
#include "chimera.h"

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
	string getCitation() { return "Haas BJ, Gevers D, Earl A, Feldgarden M, Ward DV, Giannokous G, Ciulla D, Tabbaa D, Highlander SK, Sodergren E, Methe B, Desantis TZ, Petrosino JF, Knight R, Birren BW (2011). Chimeric 16S rRNA sequence formation and detection in Sanger and 454-pyrosequenced PCR amplicons. Genome Res. \nhttp://www.mothur.org/wiki/Chimera.slayer"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
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
	
	int driver(linePair*, string, string, string, string);
	int createProcesses(string, string, string, string);
	int divideInHalf(Sequence, string&, string&);
	map<string, int> sortFastaFile(string, string);
		
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
	#endif

	bool abort, realign, trim, trimera, save;
	string fastafile, templatefile, outputDir, search, namefile, blastlocation;
	int processors, window, iters, increment, numwanted, ksize, match, mismatch, parents, minSimilarity, minCoverage, minBS, minSNP, numSeqs, templateSeqsLength;
	float divR;
	Chimera* chimera;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	
};

/***********************************************************/

#endif


