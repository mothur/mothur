#ifndef COLLECTDISPLAY_H
#define COLLECTDISPLAY_H

#include "calculator.h"
#include "fileoutput.h"
#include "display.h"

/*  There is a display for each calculator. The CollectorsCurveData class manages
    the displays. It sends each one either a sabundvector (collect.single) or a
    set of samples for a shared file (collect.shared) to find the results for. The
    display arranges the data and sends it to the FileOutput class for writing.
 */
/***********************************************************************/

class CollectDisplay : public Display {
	
public:
	CollectDisplay(Calculator* calc, FileOutput* file) : estimate(calc), output(file) { timesCalled = 0; }
	~CollectDisplay()	{	delete estimate; delete output;		}
	
    //used by collect.single
	void update(SAbundVector& rank){
		nSeqs=rank.getNumSeqs();
		data = estimate->getValues(&rank);
		output->updateOutput(nSeqs, data);
	}
	
    /* This function is called by the collect class. The collect class is passing pairs of samples,
        as well as the all samples if a multi calc is used. This function assembles a row of output data.
        It makes sure the output is assembled in the same order as the labels in the header column.
        It then sends the entire row to the file output class to handle the file writing. */
    void update(vector<SharedRAbundVector*>& shared, int numSeqs, bool pairs, map<string, int> groupComboToColumn){
        timesCalled++;
        data = estimate->getValues(shared);  //passes estimators a shared vector from each group to be compared
        
        //figure out what groups are being compared in getValues
        //because we randomizes the order we need to put the results in the correct column in the output file
        //pos tells you which column in the output file you are in
        string groupComboName = shared[0]->getGroup() +"_"+ shared[1]->getGroup();
        numGroupComb = groupComboToColumn.size();
        if (!pairs && all) { groupComboName = "all";  groupComboToColumn[groupComboName] = numGroupComb; numGroupComb++; }
        map<string, int>::iterator it = groupComboToColumn.find(groupComboName);
        int pos = 0;
        if (it != groupComboToColumn.end()) {
            pos = it->second * data.size(); //combo location * 1, or comboLocation * 3 for lci/hci
        }else {
            cout << groupComboName << " shouldn't get here\n";
        }
        
        //fills groupdata with datas info
        groupData.resize((numGroupComb*data.size()), 0);
        for (int i = 0; i < data.size(); i++) { groupData[pos+i] = data[i]; }
        
		//when you get all your groups info then output
		if ((timesCalled % numGroupComb) == 0) {
			output->updateOutput(numSeqs, groupData);
		}
	}
									
	void init(string s)		{	output->setLabelName(s);	}
	void reset()			{	output->resetFile();        }
	void close()			{	output->resetFile();        }
	void setAll(bool a)		{	all = a;                    }
	bool getAll()			{	return all;                 }
	string getName()        {  return estimate->getName();  }
	
	bool isCalcMultiple()	{ return estimate->getMultiple(); }
	bool calcNeedsAll()     { return estimate->getNeedsAll(); }
	bool hasLciHci()	{
		if (estimate->getCols() == 3) { return true; } 
		else{ return false; } 
	}
	
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
