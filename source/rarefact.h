#ifndef RAREFACT_H
#define RAREFACT_H

#include "rarefactioncurvedata.h"
#include "raredisplay.h"
#include "ordervector.hpp"
#include "mothur.h"


class Rarefact {
	
public:
	Rarefact(OrderVector& o, vector<Display*> disp, int p, set<int> en) :
			numSeqs(o.getNumSeqs()), order(o), displays(disp), label(o.getLabel()), processors(p), ends(en)  { m = MothurOut::getInstance(); }
	Rarefact(vector<RAbundVector*> shared, vector<Display*> disp) :
					 lookup(shared), displays(disp) {  m = MothurOut::getInstance(); }

	~Rarefact(){};
	int getCurve(float, int);
	int getSharedCurve(float, int);
	
private:
	
	OrderVector order;
	vector<Display*> displays;
	int numSeqs, numGroupComb, processors;
	string label;
    set<int> ends;
	void mergeVectors(RAbundVector*, RAbundVector*);
	vector<RAbundVector*> lookup;
	MothurOut* m;
	
	int createProcesses(vector<int>&, RarefactionCurveData*, int, int);
	int driver(RarefactionCurveData*, int, int);

};


#endif

