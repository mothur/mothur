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
	
	void update(SharedRAbundVector* shared1, SharedRAbundVector* shared2, int numSeqs, int numGroups){
		timesCalled++;
		data = estimate->getValues(shared1, shared2);  //passes estimators a shared vector from each group to be compared
		
		//figure out what groups are being compared in getValues
		//because the jumble parameter randomizes the order we need to put the results in the correct column in the output file
		int group1Index, group2Index, pos;
		group1Index = shared1->getGroupIndex();
		group2Index = shared2->getGroupIndex();
		
		numGroupComb = 0;
		int n = 1;
		for (int i = 0; i < (numGroups - 1); i++) {
			for (int l = n; l < numGroups; l++) {
				if ((group1Index == i) && (group2Index == l)) {
					pos = numGroupComb;  //pos tells you which column in the output file you are in
				}else if ((group1Index == l) && (group2Index == i)) {
					pos = numGroupComb;
				}
				numGroupComb++;
			}
			n++;
		}

		groupData.resize((numGroupComb*data.size()), 0);
		
		//fills groupdata with datas info
		for (int i = 0; i < data.size(); i++) {
			groupData[pos+i] = data[i];
		}
		
		//when you get all your groups info then output
		if ((timesCalled % numGroupComb) == 0) {
			output->output(numSeqs, groupData);	
		}
	};
	
	void init(string s)		{	output->initFile(s);	};
	void reset()			{	output->resetFile();	};
	void close()			{	output->resetFile();	};
	
private:
	Calculator* estimate;
	FileOutput* output;
	int nSeqs, timesCalled, numGroupComb;
	vector<double> data;
	vector<double> groupData;
};

/***********************************************************************/

#endif
