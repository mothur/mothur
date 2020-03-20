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
#include "phylotree.h"
#include "inputdata.h"

class MergeOTUsCommand : public Command {
    
public:
    MergeOTUsCommand(string);
    ~MergeOTUsCommand();
    
    vector<string> setParameters();
    string getCommandName()			{ return "merge.otus";				}
    string getCommandCategory()		{ return "OTU-Based Approaches";	}
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Merge.otus"; }
    string getDescription()		{ return "combine otus based on consensus taxonomy"; }
    
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines;
    string label, outputDir, constaxfile, sharedfile, listfile, relabundfile;
    int taxLevelCutoff;
    vector<string> Groups, outputNames;
    set<string> labels;
    map<string, string> otuLabel2ConsTax;
    map<string, int> otuLabel2ConsSize; //for use with list file, since list file only contains uniques
    
    int mergeListOTUs(vector<TaxNode>&);
    int mergeSharedOTUs(vector<TaxNode>&);
    int mergeRelabundOTUs(vector<TaxNode>&);
    int process(SharedRAbundVectors*&, ofstream&, bool&, vector<TaxNode>& nodes);
    int process(ListVector*&, ofstream&, bool&, vector<TaxNode>& nodes);
    int process(SharedRAbundFloatVectors*&, ofstream&, bool&, vector<TaxNode>& nodes);

};


#endif /* mergeotuscommand_hpp */
