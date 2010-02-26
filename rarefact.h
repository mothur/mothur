#ifndef RAREFACT_H
#define RAREFACT_H

#include "rarefactioncurvedata.h"
#include "raredisplay.h"
#include "ordervector.hpp"
#include "mothur.h"
#include "globaldata.hpp"


class Rarefact {
	
public:
	Rarefact(OrderVector* o, vector<Display*> disp) :
			numSeqs(o->getNumSeqs()), order(o), displays(disp), label(o->getLabel())  { m = MothurOut::getInstance(); }
	Rarefact(vector<SharedRAbundVector*> shared, vector<Display*> disp) :
					 lookup(shared), displays(disp) {  globaldata = GlobalData::getInstance(); m = MothurOut::getInstance(); }

	~Rarefact(){};
	void getCurve(int, int);
	void getSharedCurve(int, int);
	
private:
	GlobalData* globaldata;
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb;
	string label;
	void mergeVectors(SharedRAbundVector*, SharedRAbundVector*);
	vector<SharedRAbundVector*> lookup; 
	MothurOut* m;

};


#endif

