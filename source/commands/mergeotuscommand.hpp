//
//  mergeotuscommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 12/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef mergeotuscommand_hpp
#define mergeotuscommand_hpp

#include "command.hpp"

class MergeOTUsCommand : public Command {
    
public:
    MergeOTUsCommand(string);
    MergeOTUsCommand();
    ~MergeOTUsCommand();
    
    vector<string> setParameters();
    string getCommandName()			{ return "merge.otus";				}
    string getCommandCategory()		{ return "OTU-Based Approaches";	}
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Merge.otus"; }
    string getDescription()		{ return "combine otus based on inputs"; }
    
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines;
    string label, outputDir, constaxfile, sharedfile, listfile, relabundfile;
    vector<string> Groups, outputNames;
    set<string> labels;
    
        
};


#endif /* mergeotuscommand_hpp */
