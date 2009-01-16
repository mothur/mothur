#ifndef COLLECTDISPLAY_H
#define COLLECTDISPLAY_H

#include "sabundvector.hpp"
#include "sharedsabundvector.h"
#include "calculator.h"
#include "fileoutput.h"
#include "display.h"
#include <vector>

using namespace std;

/***********************************************************************/

class CollectDisplay : public Display {
	
public:
	CollectDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file) {timesCalled = 0;};
	~CollectDisplay()	{	delete estimate; delete output;		}
	void update(SAbundVector* rank){
		nSeqs=rank->getNumSeqs();
		data = estimate->getValues(rank);
		output->output(nSeqs, data);	
	};
	
	void update(SharedRAbundVector* shared1, SharedRAbundVector* shared2, int numSeqs, int numGroupComb){
		timesCalled++;
		data = estimate->getValues(shared1, shared2);  //passes estimators a shared vector from each group to be compared
		//fills groupdata with datas info
		for (int i = 0; i < data.size(); i++) {
			groupData.push_back(data[i]);
		}
		//when you get all your groups info then output
		if ((timesCalled % numGroupComb) == 0) {
			output->output(numSeqs, groupData);	
			groupData.clear();
		}
	};
	
	void init(string s)		{	output->initFile(s);	};
	void reset()			{	output->resetFile();	};
	void close()			{	output->resetFile();	};
	
private:
	Calculator* estimate;
	FileOutput* output;
	int nSeqs, timesCalled;
	vector<double> data;
	vector<double> groupData;
};

/***********************************************************************/

#endif