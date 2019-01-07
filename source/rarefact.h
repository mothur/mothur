#ifndef RAREFACT_H
#define RAREFACT_H

#include "rarefactioncurvedata.h"
#include "raredisplay.h"
#include "ordervector.hpp"
#include "mothur.h"


class Rarefact {
	
public:
	Rarefact(OrderVector& o, vector<Display*> disp, set<int> en, int proc) :
    numSeqs(o.getNumSeqs()), order(o), displays(disp), label(o.getLabel()),  ends(en)  { m = MothurOut::getInstance(); jumble = false; processors = proc; }
	Rarefact(vector<SharedRAbundVector*> shared, vector<Display*> disp, bool j) :
					 lookup(shared), displays(disp), jumble(j) {  m = MothurOut::getInstance(); }

	~Rarefact(){};
	int getCurve(float, int);
	int getSharedCurve(float, int);
	
private:
	
	OrderVector order;
	vector<Display*> displays;
	int numSeqs, numGroupComb, processors;
	string label;
    set<int> ends;
	void mergeVectors(SharedRAbundVector*, SharedRAbundVector*);
	vector<SharedRAbundVector*> lookup;
	MothurOut* m;
    bool jumble;
    Utils util;
	
	int driver(vector<Display*>&, int, int);

};


#endif

