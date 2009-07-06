#ifndef CHIMERACOMMAND_H
#define CHIMERACOMMAND_H

/*
 *  chimeraseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "filterseqscommand.h"
#include "sequence.hpp"
#include "sparsematrix.hpp"
#include "dist.h"

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;  //maps sequence to all distance for that seqeunce

/***********************************************************/

class ChimeraSeqsCommand : public Command {
public:
	ChimeraSeqsCommand(string);
	~ChimeraSeqsCommand();
	int execute();
	void help();
	
private:
	//Dist* distCalculator;
	
	struct Preference {
		string leftParent;
		string rightParent;
		float score;

	};

	Dist* distCalculator;
	bool abort;
	string method, fastafile;
	bool filter, correction;
	int processors, midpoint;
	FilterSeqsCommand* filterSeqs;
	vector<Sequence> seqs;
	vector<Preference> pref;
	
	int findAverageMidPoint();
	void readSeqs();
	void generatePreferences(SparseMatrix*, SparseMatrix*);
	int createSparseMatrix(int, int, SparseMatrix*, vector<Sequence>);
	

};

/***********************************************************/

#endif

