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
	
	void registerDisplay(Display* o)            {	displays.insert(o);				}
    void registerDisplays(vector<Display*> o)	{	for(int i=0;i<o.size();i++){ registerDisplay(o[i]); 	} }
	void updateRankData(SAbundVector& rv)       {	rank = rv; for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){ (*pos)->update(rank); }	}
	
private:
	set<Display*> displays;
	SAbundVector rank;
	
};

/***********************************************************************/


class SharedCollectorsCurveData : public Observable {
	
public:
	SharedCollectorsCurveData() {}
	
	void registerDisplay(Display* o)		    {	displays.insert(o);			                       }
    void registerDisplays(vector<Display*> o)	{	for(int i=0;i<o.size();i++){ displays.insert(o[i]); } }
    
    void updateSharedData(vector<SharedRAbundVector*>& shared, int numSeqs, map<string, int>& groupComboToColumn)	{
        
       
        for (int k = 0; k < (shared.size() - 1); k++) { // pass cdd each set of groups to commpare
            for (int l = k+1; l < shared.size(); l++) {
                
                for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
                    vector<SharedRAbundVector*> subset;
                    //add new pair of sharedrabund vectors
                    subset.push_back(shared[k]); subset.push_back(shared[l]);
                    if ((*pos)->calcNeedsAll()) {
                        //load subset with rest of lookup for those calcs that need everyone to calc for a pair
                        for (int w = 0; w < shared.size(); w++) {
                            if ((w != k) && (w != l)) { subset.push_back(shared[w]); }
                        }
                        (*pos)->update(subset, numSeqs, true, groupComboToColumn);
                    }else {
                        (*pos)->update(subset, numSeqs, true, groupComboToColumn);
                    }
                }
            }
        }
        
        //if this is a calculator that can do multiples then do them
        for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
            if ((*pos)->isCalcMultiple() && (*pos)->getAll()) {
                (*pos)->update(shared, numSeqs, false, groupComboToColumn);
            }
        }
    }
	
private:
	set<Display*> displays;
};

/***********************************************************************/

#endif

