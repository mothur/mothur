#ifndef RAREFACT_H
#define RAREFACT_H

#include "rarefactioncurvedata.h"
#include "raredisplay.h"
#include "ordervector.hpp"
#include "mothur.h"
#include "globaldata.hpp"


class Rarefact {
	
public:
	Rarefact(OrderVector* o, vector<Display*> disp, int p) :
			numSeqs(o->getNumSeqs()), order(o), displays(disp), label(o->getLabel()), processors(p)  { m = MothurOut::getInstance(); }
	Rarefact(vector<SharedRAbundVector*> shared, vector<Display*> disp) :
					 lookup(shared), displays(disp) {  globaldata = GlobalData::getInstance(); m = MothurOut::getInstance(); }

	~Rarefact(){};
	int getCurve(float, int);
	int getSharedCurve(float, int);
	
private:
	GlobalData* globaldata;
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb, processors;
	string label;
	void mergeVectors(SharedRAbundVector*, SharedRAbundVector*);
	vector<SharedRAbundVector*> lookup; 
	MothurOut* m;
	
	int createProcesses(vector<int>&, RarefactionCurveData*, int);
	int driver(RarefactionCurveData*, int, int);

};


#endif

