#ifndef DISPLAY_H
#define DISPLAY_H

#include "calculator.h"
#include "fileoutput.h"

/***********************************************************************/

class Display {
	
public:
	virtual void update(SAbundVector* rank) = 0;
	virtual void update(vector<RAbundVector*> shared, int numSeqs, int numGroupComb);
	virtual void init(string) = 0;
	virtual void reset() = 0;
	virtual void close() = 0;
	virtual void outputTempFiles(string) {}
	virtual void inputTempFiles(string) {}
	virtual bool isCalcMultiple() = 0;
	virtual void setAll(bool){}
	virtual bool hasLciHci(){ return false; }
	virtual bool getAll()	{	bool a; return a;	}
	virtual bool calcNeedsAll()    { bool a; return a;	}
	virtual string getName() { return ""; };
	virtual ~Display() {}
	Display() {  m = MothurOut::getInstance();  }
	
protected:
	MothurOut* m;
	
};

/***********************************************************************/

#endif
