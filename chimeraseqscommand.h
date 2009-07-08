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

struct Preference {
		string name;
		vector<string> leftParent; //keep the name of closest left associated with the two scores
		vector<string> rightParent; //keep the name of closest right associated with the two scores
		vector<float> score;  //so you can keep last score and calc this score and keep whichever is bigger.
		vector<float> closestLeft;  //keep the closest left associated with the two scores
		vector<float> closestRight; //keep the closest right associated with the two scores
		int midpoint;

};



/***********************************************************/

class ChimeraSeqsCommand : public Command {
public:
	ChimeraSeqsCommand(string);
	~ChimeraSeqsCommand();
	int execute();
	void help();
	
		
private:
	
	Dist* distCalculator;
	bool abort;
	string method, fastafile;
	bool filter, correction;
	int processors, midpoint, averageLeft, averageRight, window, iters, increment;
	FilterSeqsCommand* filterSeqs;
	ListVector* list;
	vector<Sequence> seqs;
	vector<Preference> pref;
	
	void readSeqs();
	void generatePreferences(vector<SeqMap>, vector<SeqMap>, int);
	int createSparseMatrix(int, int, SparseMatrix*, vector<Sequence>);
	

};

/***********************************************************/

#endif

