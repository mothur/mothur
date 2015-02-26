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
	FlowData(int, float, float, int, string);
	~FlowData();
	bool getNext(ifstream&);
	string getName();
	void capFlows(int);
	bool hasMinFlows(int);
    bool hasGoodHomoP();

	Sequence getSequence();

	void printFlows(ofstream&);
	void printFlows(ofstream&, string);
private:
	MothurOut* m;
	
	void updateEndFlow();
	void translateFlow();
	float signalIntensity, noiseIntensity;
	int maxHomoP;
	string seqName, locationString, sequence, baseFlow;
	int numFlows, maxFlows, endFlow;
	vector<float> flowData;
    string getSequenceName(ifstream&);
};

#endif
