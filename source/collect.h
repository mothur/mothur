#ifndef COLLECT_H
#define COLLECT_H

#include "collectorscurvedata.h"
#include "display.h"
#include "ordervector.hpp"
#include "sharedordervector.h"

/***********************************************************************/

class Collect {
	
public:
	Collect(OrderVector* order, vector<Display*> disp) :
					numSeqs(order->getNumSeqs()), order(order), displays(disp), label(order->getLabel())  { m = MothurOut::getInstance(); };
	Collect(SharedOrderVector* sharedorder, vector<Display*> disp) :
					numSeqs(sharedorder->getNumSeqs()), sharedorder(sharedorder), displays(disp), label(sharedorder->getLabel())  { m = MothurOut::getInstance(); }

	~Collect(){		};
	int getCurve(float);
	int getSharedCurve(float);
	
private:
    MothurOut* m;
	SharedOrderVector* sharedorder;
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb;
	string label, groupLabel;
	vector<string> groupComb;
    
	bool validGroup(vector<string>, string);
    map<string, int> getGroupComb(vector<string>);
	
};


#endif

