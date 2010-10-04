#ifndef RAREDISPLAY_H
#define RAREDISPLAY_H

#include "sabundvector.hpp"
#include "calculator.h"
#include "fileoutput.h"
#include "display.h"

/***********************************************************************/

class RareDisplay : public Display {
	
public:
	RareDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file), nIters(1), index(0) {};
	~RareDisplay()					{	delete estimate; delete output;		};
	void init(string);
	void reset();
	void update(SAbundVector*);
	void update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroupComb);
	void close();
	bool isCalcMultiple() { return estimate->getMultiple(); }
	
	void outputTempFiles(string);
	void inputTempFiles(string);
	
private:
	Calculator* estimate;
	FileOutput* output;
	string label;
	vector<int> seqs;  
	vector<double> results;
	vector<double> var;
	int index, nIters;
};

#endif

