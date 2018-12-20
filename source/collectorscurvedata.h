#ifndef COLLECTORSCURVEDATA_H
#define COLLECTORSCURVEDATA_H

#include "sabundvector.hpp"
#include "sharedrabundvectors.hpp"
#include "display.h"
#include "observable.h"


/***********************************************************************/

class CollectorsCurveData : public Observable {
	
public:
	CollectorsCurveData() : rank(0) {};
	
	void registerDisplay(Display* o)		{	displays.insert(o);				};
	//void removeDisplay(Display* o)			{	displays.erase(o);	delete o;	};
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
	
	void registerDisplay(Display* o)		{	displays.insert(o);			};
	//void removeDisplay(Display* o)			{	displays.erase(o);	delete o;	};
	void SharedDataChanged()				{	notifyDisplays();				};
    void updateSharedData(vector<SharedRAbundVector*> s, int numSeqs, int numGroupComb, bool p, vector<string> g)	{	pairs = p; shared = s; NumSeqs = numSeqs; NumGroupComb = numGroupComb; groups = g; SharedDataChanged();	};
	
	void notifyDisplays(){	
		for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
			
			if ((*pos)->calcNeedsAll() ) {
				(*pos)->update(shared, NumSeqs, NumGroupComb, groups);
			}else{
				
				if ( ((*pos)->isCalcMultiple() ) && ((*pos)->getAll() ) && (!pairs) ) {
					(*pos)->update(shared, NumSeqs, NumGroupComb, groups);
				}else {
					vector<SharedRAbundVector*> temp; temp.push_back(shared[0]); temp.push_back(shared[1]);
					shared = temp; 
					
					(*pos)->update(shared, NumSeqs, NumGroupComb, groups);
				}
			}
		}	
	};
	
private:
	set<Display*> displays;
	vector<Display*> multiDisplays;
	vector<SharedRAbundVector*> shared;
    vector<string> groups;
	int NumSeqs, NumGroupComb;
	bool pairs;
};

/***********************************************************************/

#endif

