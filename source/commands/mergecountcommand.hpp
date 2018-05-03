//
//  mergecountcommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/3/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//

#ifndef mergecountcommand_hpp
#define mergecountcommand_hpp

#include "command.hpp"

class MergeCountCommand : public Command {
    
#ifdef UNIT_TEST
    //friend class TestMergeCountCommand;
#endif
    
public:
    MergeCountCommand(string);
    MergeCountCommand();
    ~MergeCountCommand() {}
    
    vector<string> setParameters();
    string getCommandName()			{ return "merge.count";	}
    string getCommandCategory()		{ return "General";			}
    
    string getHelpString();
    string getOutputPattern(string) { return "";  }
    string getCitation() { return "http://www.mothur.org/wiki/Merge.count"; }
    string getDescription()		{ return "reads count files and combines them into a single count file"; }
    
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    
    bool abort;
    string outputDir, inputDir, countfile, output, outputFileName;
    vector<string> outputNames, fileNames;
    int numInputFiles;
    
};


#endif /* mergecountcommand_hpp */
