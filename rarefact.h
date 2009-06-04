#ifndef RAREFACT_H
#define RAREFACT_H

#include "rarefactioncurvedata.h"
#include "raredisplay.h"
#include "ordervector.hpp"
#include "mothur.h"


class Rarefact {
	
public:
	Rarefact(OrderVector* o, vector<Display*> disp) :
			numSeqs(o->getNumSeqs()), order(o), displays(disp), label(o->getLabel())  {};
	Rarefact(vector<SharedRAbundVector*> shared, vector<Display*> disp) :
					 lookup(shared), displays(disp) {};

	~Rarefact(){};
	void getCurve(int, int);
	void getSharedCurve(int, int);
	
private:
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb;
	string label;
	void mergeVectors(SharedRAbundVector*, SharedRAbundVector*);
	vector<SharedRAbundVector*> lookup; 

};


#endif

