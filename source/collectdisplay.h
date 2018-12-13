#ifndef COLLECTDISPLAY_H
#define COLLECTDISPLAY_H

#include "calculator.h"
#include "fileoutput.h"
#include "display.h"


/***********************************************************************/

class CollectDisplay : public Display {
	
public:
	CollectDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file) {timesCalled = 0;};
	~CollectDisplay()	{	delete estimate; delete output;		}
	
	
	void update(SAbundVector* rank){
		nSeqs=rank->getNumSeqs();
		data = estimate->getValues(rank);
		output->updateOutput(nSeqs, data);
	};
	
    void update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroups, vector<string> mGroups){
        timesCalled++;
        data = estimate->getValues(shared);  //passes estimators a shared vector from each group to be compared
        
        //figure out what groups are being compared in getValues
        //because we randomizes the order we need to put the results in the correct column in the output file
        int group1Index, group2Index, pos;
        
        for (int i = 0; i < mGroups.size(); i++) {
            if (shared[0]->getGroup() == mGroups[i]) { group1Index = i; }
            if (shared[1]->getGroup() == mGroups[i]) { group2Index = i; }
        }
        
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
			
		if ((estimate->getMultiple() ) && all) { 
			numGroupComb++; 
			groupData.resize((numGroupComb*data.size()), 0);
			//is this the time its called with all values
			if  ((timesCalled % numGroupComb) == 0) { 
				//last spot
				pos = ((groupData.size()-1) * data.size());
			}
			//fills groupdata with datas info
			for (int i = 0; i < data.size(); i++) {
				groupData[pos+i] = data[i];
			}
		}else {
			groupData.resize((numGroupComb*data.size()), 0);
			//fills groupdata with datas info
			for (int i = 0; i < data.size(); i++) {
				groupData[pos+i] = data[i];
			}
		}
		
		//when you get all your groups info then output
		if ((timesCalled % numGroupComb) == 0) {
			output->updateOutput(numSeqs, groupData);
		}
	};
									
	void init(string s)		{	output->setLabelName(s);	}
	void reset()			{	output->resetFile();        }
	void close()			{	output->resetFile();        }
	void setAll(bool a)		{	all = a;                    }
	bool getAll()			{	return all;                 }
	
	
	bool isCalcMultiple()	{ return estimate->getMultiple(); }
	bool calcNeedsAll()     { return estimate->getNeedsAll(); }
	bool hasLciHci()	{
		if (estimate->getCols() == 3) { return true; } 
		else{ return false; } 
	}
	
	string getName()	{  return estimate->getName();  }
	
	
private:
	
	Calculator* estimate;
	FileOutput* output;
	int nSeqs, timesCalled, numGroupComb;
	vector<double> data;
	vector<double> groupData;
	bool all;
	
};

/***********************************************************************/

#endif
