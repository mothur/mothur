#ifndef FLOWDATA_H
#define FLOWDATA_H

/*
 *  flowdata.h
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"

class FlowData {

public:
	FlowData();
	FlowData(ifstream&, float, float, int);
	~FlowData(){};
	void capFlows(int);
	bool hasMinFlows(int);
	Sequence getSequence();
	
	int getSeqLength();
	void printFlows(ofstream&);
	void printFlows(ofstream&, string);
	void printFASTA(ofstream&);
private:
	MothurOut* m;

	void findDeadSpot(float, float, int);
	void translateFlow();
	
	string seqName, locationString, sequence, baseFlow;
	int numFlows, seqLength, deadSpot;
	vector<float> flowData;
};

#endif
