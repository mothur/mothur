//
//  biominfocommand.h
//  Mothur
//
//  Created by Sarah Westcott on 8/5/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__biominfocommand__
#define __Mothur__biominfocommand__

#include "command.hpp"
#include "inputdata.h"
#include "phylosummary.h"

#define MAX_NAME 1024

class BiomInfoCommand : public Command {
    
#ifdef UNIT_TEST
    friend class TestBiomInfoCommand;
#endif
    
public:
    BiomInfoCommand(string);
    ~BiomInfoCommand() {}
    
    vector<string> setParameters();
    string getCommandName()			{ return "biom.info";				}
    string getCommandCategory()		{ return "OTU-Based Approaches";	}
    
    string getCommonQuestions();
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Biom.info"; }
    string getDescription()		{ return "create 'mothur' files from a biom file. ie: shared, taxonomy, constaxonomy"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
protected:
    void printSharedData(SharedRAbundVectors*, ofstream&);
    int createFilesFromBiom();
    int extractFilesFromHDF5();
    string getTag(string&);
    string getName(string);
    string getTaxonomy(string, string);
    vector< vector<string> > readRows(string, int&, bool&);
    int getDims(string, int&, int&);
    SharedRAbundVectors* readData(string, string, string, vector<string>&, int);
    vector<string> getNamesAndTaxonomies(string);
    
    vector<string> outputNames, otuNames, sampleNames, taxonomy;
    vector<int> indices, indptr, otudata;
    string fileroot, outputDir, biomfile, label, basis, output, format;
    bool firsttime, abort, relabund;
    int maxLevel, printlevel, nnz;

    #ifdef USE_HDF5
    void processAttributes(H5::Group&, set<string>&);
    void checkGroups(H5::H5File&, map<string, vector<string> >&);
    #endif
};


#endif /* defined(__Mothur__biominfocommand__) */
