//
//  makeclrcommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/20/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef makeclrcommand_hpp
#define makeclrcommand_hpp

#include "command.hpp"
#include "inputdata.h"

class MakeCLRCommand : public Command {

public:
    MakeCLRCommand(string);
    ~MakeCLRCommand() {}
    
    vector<string> setParameters();
    string getCommandName()            { return "make.clr";                 }
    string getCommandCategory()        { return "OTU-Based Approaches";     }
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation()            { return "http://www.mothur.org/wiki/Make.clr";                 }
    string getDescription()         { return "create a log centered ratio file from a shared file"; }

    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
    
private:
    bool abort, allLines;
    set<string> labels; //holds labels to be used
    string groups, label,  sharedfile, zeroReplacement;
    vector<string> Groups, outputNames;
    double zeroReplacementValue;
        
    void process(SharedRAbundVectors*&, ofstream&, bool&);
    
};

#endif /* makeclrcommand_hpp */
