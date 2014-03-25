//
//  getmimarkspackagecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 3/25/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#ifndef Mothur_getmimarkspackagecommand_h
#define Mothur_getmimarkspackagecommand_h

#include "command.hpp"

/**************************************************************************************************/

class GetMIMarksPackageCommand : public Command {
public:
    GetMIMarksPackageCommand(string);
    GetMIMarksPackageCommand();
    ~GetMIMarksPackageCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.mimarkspackage";			}
    string getCommandCategory()		{ return "Sequence Processing";         }
    
    string getOutputPattern(string);
	string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/get.mimarkspackage"; }
    string getDescription()		{ return "create blank mimarks package form for sra command"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort;
    string oligosfile, groupfile, package;
    string outputDir;
    vector<string> outputNames;
};

/**************************************************************************************************/




#endif
