#ifndef RAREDISPLAY_H
#define RAREDISPLAY_H

#include "sabundvector.hpp"
#include "calculator.h"
#include "fileoutput.h"
#include "display.h"

/***********************************************************************/
//Each display is responsible for one calculator. The FileOutput class handles creating the outputfile for the calc.
//This class uses mutex and lock_guard to prevent thread errors.

class RareDisplay : public Display {
	
public:
	RareDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file), nIters(1) {};
	~RareDisplay()					{	delete estimate; delete output;		}
	void init(string);
	void reset();
	void update(SAbundVector&);
	void update(vector<SharedRAbundVector*> shared, int numSeqs);
	void close();
	bool isCalcMultiple() { return estimate->getMultiple(); }
	
private:
	Calculator* estimate;
	FileOutput* output;
	string label;
	map<int, vector<double> > results; //maps seqCount to results for that number of sequences
	int nIters;
    Utils util;
    std::mutex mutex;
    
};

#endif

