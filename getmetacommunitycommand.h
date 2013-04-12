//
//  getmetacommunitycommand.h
//  Mothur
//
//  Created by SarahsWork on 4/9/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_getmetacommunitycommand_h
#define Mothur_getmetacommunitycommand_h

#include "command.hpp"
#include "inputdata.h"

/**************************************************************************************************/

class GetMetaCommunityCommand : public Command {
public:
    GetMetaCommunityCommand(string);
    GetMetaCommunityCommand();
    ~GetMetaCommunityCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.metacommunity";		}
    string getCommandCategory()		{ return "OTU-Based Approaches";         }
    
    string getOutputPattern(string);
    
	string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/Get.metacommunity"; }
    string getDescription()		{ return "brief description"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines;
    string outputDir;
    vector<string> outputNames;
    string sharedfile;
    int minpartitions, maxpartitions, optimizegap;
    vector<string> Groups;
    set<string> labels;
    
    int process(vector<SharedRAbundVector*>&);
    vector<double> generateDesignFile(int, string);
    int generateSummaryFile(int, string, string, string, string);

};

/**************************************************************************************************/
struct summaryData {
    
    string name;
    double refMean, difference;
    vector<double> partMean, partLCI, partUCI;
    
};
/**************************************************************************************************/




#endif
