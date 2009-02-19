#ifndef SUMMARY_H
#define SUMMARY_H

#include "mothur.h"
#include "display.h"



class SummaryDisplay : public Display {

public:
	SummaryDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file){	};
	
	void update(SAbundVector* sabund)	{	output->output(sabund->getLabel(), estimate->getValues(rank));	};
	void init(string s)				{	output->initFile(s);	};
	void reset()					{	output->resetFile();	};
	void close(){};
	string getLabel(){	return estimate->getName();	}
	
private:
	Calculator* estimate;
	FileOutput* output;
	
	int nIters;
	string tempInName, tempOutName;
	ifstream tempInFile;
	ofstream tempOutFile;
};

/***********************************************************************/

#endif
