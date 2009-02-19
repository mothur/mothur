#ifndef COLLECTORSCURVEDATA_H
#define COLLECTORSCURVEDATA_H

#include "mothur.h"
#include "sabundvector.hpp"
#include "sharedrabundvector.h"
#include "display.h"
#include "observable.h"

using namespace std;

/***********************************************************************/

class CollectorsCurveData : public Observable {
	
public:
	CollectorsCurveData() : rank(0) {};
	
	void registerDisplay(Display* o)		{	displays.insert(o);				};
	void removeDisplay(Display* o)			{	displays.erase(o);	delete o;	};
	SAbundVector* getRankData()				{	return rank;					};
	void rankDataChanged()					{	notifyDisplays();				};
	void updateRankData(SAbundVector* rv)	{	rank = rv; rankDataChanged();	};

	void notifyDisplays(){	
		for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
			(*pos)->update(rank);
		}	
	};
	
private:
	set<Display*> displays;
	SAbundVector* rank;
	
};

/***********************************************************************/


class SharedCollectorsCurveData : public Observable {
	
public:
	SharedCollectorsCurveData() : shared1(0), shared2(0) {};
	
	void registerDisplay(Display* o)		{	displays.insert(o);				};
	void removeDisplay(Display* o)			{	displays.erase(o);	delete o;	};
	void SharedDataChanged()					{	notifyDisplays();				};
	void updateSharedData(SharedRAbundVector* rv, SharedRAbundVector* rv2, int numSeqs, int numGroupComb)	{	shared1 = rv; shared2 = rv2; NumSeqs = numSeqs; NumGroupComb = numGroupComb; SharedDataChanged();	};

	void notifyDisplays(){	
		for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
			(*pos)->update(shared1, shared2, NumSeqs, NumGroupComb);
		}	
	};
	
private:
	set<Display*> displays;
	SharedRAbundVector* shared1;
	SharedRAbundVector* shared2;
	int NumSeqs, NumGroupComb;
};

/***********************************************************************/

#endif

