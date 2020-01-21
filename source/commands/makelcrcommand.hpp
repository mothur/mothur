//
//  makelcrcommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/20/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef makelcrcommand_hpp
#define makelcrcommand_hpp

#include "command.hpp"
#include "inputdata.h"

class MakeLCRCommand : public Command {

public:
    MakeLCRCommand(string);
    MakeLCRCommand();
    ~MakeLCRCommand() {}
    
    vector<string> setParameters();
    string getCommandName()            { return "make.lcr";                 }
    string getCommandCategory()        { return "OTU-Based Approaches";     }
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation()            { return "http://www.mothur.org/wiki/Make.lcr";                 }
    string getDescription()         { return "create a log centered ratio file from a shared file"; }

    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
    
private:
    bool abort, allLines;
    set<string> labels; //holds labels to be used
    string groups, label, outputDir, sharedfile, zeroReplacement;
    vector<string> Groups, outputNames;
    double zeroReplacementValue;
        
    void process(SharedRAbundVectors*&);
    
};

#endif /* makelcrcommand_hpp */
