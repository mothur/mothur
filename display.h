#ifndef DISPLAY_H
#define DISPLAY_H

#include "sabundvector.hpp"
#include "sharedsabundvector.h"
#include "calculator.h"
#include "fileoutput.h"


using namespace std;

/***********************************************************************/

class Display {
	
public:
	virtual void update(SAbundVector* rank) = 0;
	virtual void update(SharedRAbundVector* shared1, SharedRAbundVector* shared2, int numSeqs, int numGroupComb) = 0;
	virtual void init(string) = 0;
	virtual void reset() = 0;
	virtual void close() = 0;
	
private:

};

/***********************************************************************/

#endif
