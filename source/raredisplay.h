#ifndef RAREDISPLAY_H
#define RAREDISPLAY_H

#include "sabundvector.hpp"
#include "calculator.h"
#include "fileoutput.h"
#include "display.h"

/***********************************************************************/

class RareDisplay : public Display {
	
public:
	RareDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file), nIters(1) {};
	~RareDisplay()					{	delete estimate; delete output;		};
	void init(string);
	void reset();
	void update(SAbundVector*);
	void update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroupComb, vector<string>);
	void close();
	bool isCalcMultiple() { return estimate->getMultiple(); }
	
	void outputTempFiles(string);
	void inputTempFiles(string);
	
private:
	Calculator* estimate;
	FileOutput* output;
	string label;
	map<int, vector<double> > results; //maps seqCount to results for that number of sequences
	int nIters;
    vector<string> Groups;
};

#endif

