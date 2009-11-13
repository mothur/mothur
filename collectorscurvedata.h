#ifndef COLLECTORSCURVEDATA_H
#define COLLECTORSCURVEDATA_H

#include "sabundvector.hpp"
#include "sharedrabundvector.h"
#include "display.h"
#include "observable.h"


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
	SharedCollectorsCurveData() { }; //: shared1(0), shared2(0)
	
	void registerDisplay(Display* o)		{	displays.insert(o);				};
	void removeDisplay(Display* o)			{	displays.erase(o);	delete o;	};
	void SharedDataChanged()				{	notifyDisplays();				};
	void updateSharedData(vector<SharedRAbundVector*> s, int numSeqs, int numGroupComb)	{	shared = s; NumSeqs = numSeqs; NumGroupComb = numGroupComb; SharedDataChanged();	};
	
	void notifyDisplays(){	
		for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
//cout << (*pos)->getName() << endl;
			if ( ( ((*pos)->isCalcMultiple() == true) && ((*pos)->getAll() == true) ) || (shared.size() == 2)  ) {
				(*pos)->update(shared, NumSeqs, NumGroupComb);
			}
		}	
	};
	
private:
	set<Display*> displays;
	vector<Display*> multiDisplays;
	vector<SharedRAbundVector*> shared;
	int NumSeqs, NumGroupComb;
};

/***********************************************************************/

#endif

