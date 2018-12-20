#ifndef RAREFACTIONCURVEDATA_H
#define RAREFACTIONCURVEDATA_H

#include "mothur.h"
#include "sabundvector.hpp"
#include "display.h"
#include "observable.h"

/***********************************************************************/
//Has a display for each calculator
class RarefactionCurveData : public Observable {
	
public:
	RarefactionCurveData() : rank(0) {};
	
	void registerDisplay(Display* o)	{	displays.insert(o);				};
    SAbundVector* getRankData()			{	return rank;					};
	void rankDataChanged()				{	notifyDisplays();				};
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

class SharedRarefactionCurveData : public Observable {
	
public:
	SharedRarefactionCurveData() {}; //: shared1(0), shared2(0) 
	
	void registerDisplay(Display* o)	{	displays.insert(o);				};
	//void removeDisplay(Display* o)		{	displays.erase(o);	delete o;	};
	void SharedDataChanged()			{	notifyDisplays();				};
	void updateSharedData(vector<SharedRAbundVector*> r, int numSeqs, int numGroupComb)	{
        shared = r; NumSeqs = numSeqs; NumGroupComb = numGroupComb;
        groups.clear(); for (int i = 0; i < r.size(); i++) { groups.push_back(r[i]->getGroup()); }
        SharedDataChanged();
    };

	void notifyDisplays(){	
		for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
				(*pos)->update(shared, NumSeqs, NumGroupComb, groups);
		}	
	};
	
private:
	set<Display*> displays;
	vector<SharedRAbundVector*> shared;
	int NumSeqs, NumGroupComb;
    vector<string> groups;
	
};

/***********************************************************************/


#endif

