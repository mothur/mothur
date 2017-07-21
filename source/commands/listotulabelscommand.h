#ifndef Mothur_listotucommand_h
#define Mothur_listotucommand_h

//
//  listotucommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "listvector.hpp"
#include "sharedrabundvectors.hpp"
#include "sharedrabundfloatvectors.hpp"

/**************************************************************************************************/

class ListOtuLabelsCommand : public Command {
public:
    ListOtuLabelsCommand(string);
    ListOtuLabelsCommand();
    ~ListOtuLabelsCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "list.otulabels";          }
    string getCommandCategory()		{ return "OTU-Based Approaches";	} 
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
    
	string getHelpString();	
    string getOutputPattern(string);	
    string getCitation() { return "http://www.mothur.org/wiki/List.otulabels"; }
    string getDescription()		{ return "lists otu labels from shared or relabund file. Can be used by get.otulabels with output from classify.otu, otu.association, or corr.axes to select specific otus."; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort, allLines;
    string outputDir, sharedfile, relabundfile, label, inputFileName, format, listfile;
    vector<string> outputNames;
    vector<string> Groups;
    set<string> labels;
    
    int createList(SharedRAbundFloatVectors*);
    int createList(SharedRAbundVectors*);
    int createList(ListVector*&);

};

/**************************************************************************************************/






#endif
