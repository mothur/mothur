#ifndef RAREFACT_H
#define RAREFACT_H

using namespace std;

#include "rarefactioncurvedata.h"
#include "raredisplay.h"
#include "ordervector.hpp"
#include "mothur.h"


class Rarefact {
	
public:
	Rarefact(OrderVector* o, vector<Display*> disp) :
			numSeqs(o->getNumSeqs()), order(o), displays(disp), label(o->getLabel())  {};
	Rarefact(SharedOrderVector* sharedorder, vector<Display*> disp) :
					numSeqs(sharedorder->getNumSeqs()), sharedorder(sharedorder), displays(disp), label(sharedorder->getLabel())  {};

	~Rarefact(){};
	void getCurve(int, int);
	void getSharedCurve(int, int);
	
private:
	SharedOrderVector* sharedorder;
	GlobalData* globaldata;
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb;
	string label;
	void mergeVectors(SharedRAbundVector*, SharedRAbundVector*);
	vector<SharedRAbundVector*> lookup; 

};


#endif

