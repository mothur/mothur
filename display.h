#ifndef DISPLAY_H
#define DISPLAY_H

#include "sabundvector.hpp"
#include "sharedsabundvector.h"
#include "calculator.h"
#include "fileoutput.h"


/***********************************************************************/

class Display {
	
public:
	virtual void update(SAbundVector* rank) = 0;
	virtual void update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroupComb) = 0;
	virtual void init(string) = 0;
	virtual void reset() = 0;
	virtual void close() = 0;
	virtual bool isCalcMultiple() = 0;
	
private:

};

/***********************************************************************/

#endif
