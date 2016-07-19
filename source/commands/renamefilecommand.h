//
//  renamefilecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 4/18/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__renamefilecommand__
#define __Mothur__renamefilecommand__

#include "command.hpp"

class RenameFileCommand : public Command {
    
#ifdef UNIT_TEST
    friend class TestRenameFileCommand;
#endif
    
public:
    
    RenameFileCommand(string);
    RenameFileCommand();
    ~RenameFileCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "rename.file";				}
    string getCommandCategory()		{ return "General";		}
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/rename.file"; }
    string getDescription()		{ return "renames file and updates current"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
    
private:
    string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile, countfile, summaryfile, inputfile, outputDir;
    string treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, flowfile, filefile, outputfile, constaxonomyfile, prefix;
    bool mothurGenerated, abort, deleteOld;
    
    vector<string> outputNames;
    
    string getNewName(string name, string type);
    string renameOrCopy(string oldName, string newName);
};

#endif /* defined(__Mothur__renamefilecommand__) */
