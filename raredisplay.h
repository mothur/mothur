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
							tempInName(getPathName(output->getFileName()) + ".tempin."+ getSimpleName(output->getFileName())), tempOutName(getPathName(output->getFileName()) + ".tempout."+ getSimpleName(output->getFileName())) {};
	~RareDisplay()					{	delete estimate; delete output;		};
	void init(string);
	void reset();
	void update(SAbundVector*);
	void update(SharedRAbundVector* shared1, SharedRAbundVector* shared2, int numSeqs, int numGroupComb);
	void close();
	
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

