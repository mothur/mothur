#ifndef COLLECT_H
#define COLLECT_H

#include "collectorscurvedata.h"
#include "display.h"
#include "ordervector.hpp"
#include "sabundvector.hpp"
#include "rabundvector.hpp"
#include "sharedordervector.h"
#include "datavector.hpp"
#include "mothurout.h"

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
	SharedOrderVector* sharedorder;
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb, totalNumSeq;
	string label, groupLabel;
	void getGroupComb();
	vector<string> groupComb;
	bool validGroup(vector<string>, string);
	MothurOut* m;
};


#endif

