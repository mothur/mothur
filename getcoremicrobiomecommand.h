#ifndef Mothur_getcoremicrobiomcommand_h
#define Mothur_getcoremicrobiomcommand_h


//
//  GetCoreMicroBiomeCommand.h
//  Mothur
//
//  Created by John Westcott on 5/8/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "inputdata.h"

/**************************************************************************************************/

class GetCoreMicroBiomeCommand : public Command {
public:
    GetCoreMicroBiomeCommand(string);
    GetCoreMicroBiomeCommand();
    ~GetCoreMicroBiomeCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.coremicrobiome";			}
    string getCommandCategory()		{ return "OTU-Based Approaches";		} 
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
    string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Get.coremicrobiome"; }
    string getDescription()		{ return "determines the fraction of OTUs that are found in varying numbers of samples for different minimum relative abundances"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    string relabundfile, sharedfile, inputFileName, format, output;
    bool allLines;
    vector<string> Groups;
    set<string> labels;
    bool abort;
    string outputDir;
    vector<string> outputNames;
    int samples, abund;
    
    int createTable(vector<SharedRAbundFloatVector*>&);

};

/**************************************************************************************************/




#endif
