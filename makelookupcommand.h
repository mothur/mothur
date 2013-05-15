//
//  makelookupcommand.h
//  Mothur
//
//  Created by SarahsWork on 5/14/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_makelookupcommand_h
#define Mothur_makelookupcommand_h

#include "command.hpp"
#include "sequence.hpp"

/**************************************************************************************************/

class MakeLookupCommand : public Command {
public:
    MakeLookupCommand(string);
    MakeLookupCommand();
    ~MakeLookupCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.lookup";			}
    string getCommandCategory()		{ return "Sequence Processing";		}
    
    string getOutputPattern(string);
	string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/Make.lookup"; }
    string getDescription()		{ return "create custom lookup files for use with shhh.flows"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort;
    string outputDir, flowFileName, errorFileName, flowOrder, refFastaFileName, barcodeSequence, keySequence;
    vector<string> outputNames;
    int thresholdCount;
    
    vector<double> convertSeqToFlow(string sequence, string order);
    int alignFlowGrams(vector<double>& flowgram, vector<double>& refFlow, double gapOpening, vector<vector<double> > penaltyMatrix, string flowOrder);
    int regress(vector<double>& data, int N);
};

/**************************************************************************************************/




#endif
