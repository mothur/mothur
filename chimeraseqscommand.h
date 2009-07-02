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


	bool abort;
	string method, fastafile;
	bool filter, correction;
	int processors, midpoint;
	FilterSeqsCommand* filterSeqs;
	vector<Sequence> seqs;
	vector<Preference> pref;
	
	int findAverageMidPoint();
	void readSeqs();
	

};

/***********************************************************/

#endif

