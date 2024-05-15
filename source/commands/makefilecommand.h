//
//  makefilecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 6/24/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__makefilecommand__
#define __Mothur__makefilecommand__

#include "command.hpp"

class MakeFileCommand : public Command {
    
public:
    MakeFileCommand(string);
    ~MakeFileCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.file";	}
    string getCommandCategory()		{ return "General";		}
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Make.file"; }
    string getDescription()		{ return "creates a file file containing fastq filenames"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }

private:
    
    string inputDir,  typeFile, prefix, delim;
    vector<string> outputNames;
    int numCols;
    bool abort;
    
    vector< vector<string> > prepareOutputs(vector< vector<string> > paired);
    vector< vector<string> > prepareOutputs(vector<string> files);
    vector<string> denoiseGroupNames(vector<string> files);
    vector<string> createDefaultGroupNames(vector<string> files);

    int fillAccnosFile(string tempFile);
    vector< vector<string> > pairFiles(vector<string> files, vector<string>& singles);
    void createFastqGzFile(vector<string>);
    void createFastaFile(vector<string>);
};


#endif /* defined(__Mothur__makefilecommand__) */
