#ifndef RAREDISPLAY_H
#define RAREDISPLAY_H

#include "sabundvector.hpp"
#include "calculator.h"
#include "fileoutput.h"
#include "display.h"


using namespace std;

/***********************************************************************/

class RareDisplay : public Display {
	
public:
	RareDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file), nIters(1),
							tempInName(getPathName(output->getFileName()) + ".tempin"), tempOutName(getPathName(output->getFileName()) + ".tempout") {};
	~RareDisplay()					{	delete estimate; delete output;		};
	void init(string);
	void reset();
	void update(SAbundVector*);
	void update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroupComb);
	void close();
	bool isCalcMultiple() { return estimate->getMultiple(); }
	
private:
	Calculator* estimate;
	FileOutput* output;
	string label;
	int nIters;
	string tempInName, tempOutName;
	ifstream tempInFile;
	ofstream tempOutFile;
	int renameOk;

};

#endif

