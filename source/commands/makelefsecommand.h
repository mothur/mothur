//
//  makelefse.h
//  Mothur
//
//  Created by SarahsWork on 6/3/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__makelefse__
#define __Mothur__makelefse__

#include "mothurout.h"
#include "command.hpp"
#include "inputdata.h"
#include "sharedutilities.h"
#include "phylosummary.h"

/**************************************************************************************************/

class MakeLefseCommand : public Command {
public:
    MakeLefseCommand(string);
    MakeLefseCommand();
    ~MakeLefseCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.lefse";			}
    string getCommandCategory()		{ return "General";		}
    
    string getOutputPattern(string);
	string getHelpString();
    string getCitation() { return "http://huttenhower.sph.harvard.edu/galaxy/root?tool_id=lefse_upload http://www.mothur.org/wiki/Make.lefse"; }
    string getDescription()		{ return "creates LEfSe input file"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines, otulabel, hasGroupInfo;
    string outputDir;
    vector<string> outputNames, Groups;
    string sharedfile, designfile, constaxonomyfile, relabundfile, scale, label, inputFile;
    
    int runRelabund(map<int, consTax2>&, SharedRAbundFloatVectors*&);
    
    SharedRAbundFloatVectors* getRelabund();
    SharedRAbundFloatVectors* getSharedRelabund();
};

/**************************************************************************************************/




#endif /* defined(__Mothur__makelefse__) */
