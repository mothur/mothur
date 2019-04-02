//
//  setseedcommand.h
//  Mothur
//
//  Created by Sarah Westcott on 3/24/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__setseedcommand__
#define __Mothur__setseedcommand__

#include "command.hpp"
#include "commandfactory.hpp"

/**********************************************************/

class SetSeedCommand : public Command {
    
public:
    SetSeedCommand(string);
    SetSeedCommand() { abort = true; calledHelp = true; setParameters(); }
    ~SetSeedCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "set.seed";		}
    string getCommandCategory()		{ return "General";		}
    
    string getHelpString();
    string getOutputPattern(string){ return ""; }
    string getCitation() { return "http://www.mothur.org/wiki/Set.seed"; }
    string getDescription()		{ return "set random seed"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort;
    int random;
    vector<string> outputNames;
    
    
};

/**********************************************************/

#endif /* defined(__Mothur__setseedcommand__) */
