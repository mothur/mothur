#ifndef COLLECT_H
#define COLLECT_H

using namespace std;


#include "collectorscurvedata.h"
#include "display.h"
#include "ordervector.hpp"
#include "sabundvector.hpp"
#include "rabundvector.hpp"
#include "sharedordervector.h"
#include "datavector.hpp"
#include "globaldata.hpp"

/***********************************************************************/

class Collect {
	
public:
	Collect(OrderVector* order, vector<Display*> disp) :
					numSeqs(order->getNumSeqs()), order(order), displays(disp), label(order->getLabel())  {};
	Collect(SharedOrderVector* sharedorder, vector<Display*> disp) :
					numSeqs(sharedorder->getNumSeqs()), sharedorder(sharedorder), displays(disp), label(sharedorder->getLabel())  {};

	~Collect(){		};
	void getCurve(int);
	void getSharedCurve(int);
	
private:
	SharedOrderVector* sharedorder;
	GlobalData* globaldata;
	OrderVector* order;
	vector<Display*> displays;
	int numSeqs, numGroupComb, totalNumSeq;
	string label, groupLabel;
	void getGroupComb();
	vector<string> groupComb;
};


#endif

