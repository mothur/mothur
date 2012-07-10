#ifndef Mothur_loadlogfilecommand_h
#define Mothur_loadlogfilecommand_h

//
//  loadlogfilecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 6/13/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"

/**************************************************************************************************/

class LoadLogfileCommand : public Command {
public:
    LoadLogfileCommand(string);
    LoadLogfileCommand();
    ~LoadLogfileCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "load.logfile";		}
    string getCommandCategory()		{ return "General";             } 
    string getOutputFileNameTag(string, string) { return ""; }
	string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Load.logfile"; }
    string getDescription()		{ return "extracts current files from a logfile"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort;
    string outputDir, logfile;
    vector<string> outputNames;
};

/**************************************************************************************************/




#endif
