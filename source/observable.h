#ifndef OBSERVABLE_H
#define OBSERVABLE_H


#include "collectdisplay.h"


/***********************************************************************/

class Observable {
	
public:
	virtual void registerDisplay(Display*) = 0;
    virtual void registerDisplays(vector<Display*>) = 0;
	virtual ~Observable() {}
};

/***********************************************************************/

#endif
