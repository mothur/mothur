#ifndef Mothur_getcoremicrobiomcommand_h
#define Mothur_getcoremicrobiomcommand_h


//
//  getcoremicrobiomcommand.h
//  Mothur
//
//  Created by John Westcott on 5/8/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "inputdata.h"

/**************************************************************************************************/

class GetCoreMicroBiomCommand : public Command {
public:
    GetCoreMicroBiomCommand(string);
    GetCoreMicroBiomCommand();
    ~GetCoreMicroBiomCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.coremicrobiom";			}
    string getCommandCategory()		{ return "OTU-Based Approaches";		} 
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
    string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/get.coremicrobio"; }
    string getDescription()		{ return "description"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    string relabundfile, sharedfile, inputFileName, format;
    bool allLines;
    vector<string> Groups;
    set<string> labels;
    bool abort;
    string outputDir;
    vector<string> outputNames;
    
    int createTable(vector<SharedRAbundFloatVector*>&);

};

/**************************************************************************************************/




#endif
