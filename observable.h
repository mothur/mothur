#ifndef OBSERVABLE_H
#define OBSERVABLE_H


#include "collectdisplay.h"

using namespace std;

/***********************************************************************/

class Observable {
	
public:
	virtual void registerDisplay(Display*) = 0;
	virtual void removeDisplay(Display*) = 0;
	virtual void notifyDisplays() = 0;	
};

/***********************************************************************/

#endif
