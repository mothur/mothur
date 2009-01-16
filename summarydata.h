#ifndef SUMMARYDATA_H
#define SUMMARYDATA_H

#include <set>
#include "sabundvector.hpp"
#include "display.h"
#include "observable.h"

using namespace std;

/***********************************************************************/

class SummaryData : public Observable {
	
public:
	SummaryData() : sabund(0) {};
	
	void registerDisplay(Display* o)	{	displays.insert(o);				};
	void removeDisplay(Display* o)		{	displays.erase(o);	delete o;	};
	SAbundVector* getSabundData()			{	return sabund;					};
	void sabundDataChanged()				{	notifyDisplays();				};
	void updatesabundData(SAbundVector* rv)	{	sabund = rv; sabundDataChanged();	};
	
	void notifyDisplays(){	
		for(set<Display*>::iterator pos=displays.begin();pos!=displays.end();pos++){
			(*pos)->update(sabund);
		}	
	};
	
private:
	set<Display*> displays;
	SAbundVector* sabund;
	
};

/***********************************************************************/

#endif
