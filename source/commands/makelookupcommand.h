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
    ~MakeLookupCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.lookup";			}
    string getCommandCategory()		{ return "Sequence Processing";		}
    
    string getOutputPattern(string);
	string getHelpString();
    string getCitation() { return "Quince, C., A. Lanzén, T. P. Curtis, R. J. Davenport, N. Hall, I. M. Head, L. F. Read, and W. T. Sloan. 2009. Accurate determination of microbial diversity from 454 pyrosequencing data. Nat Methods 6:639-41. http://www.mothur.org/wiki/Make.lookup"; }
    string getDescription()		{ return "Creates a lookup file for use with shhh.flows using user-supplied mock community data and flow grams"; }
    
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
