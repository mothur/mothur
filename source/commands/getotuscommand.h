#ifndef Mothur_getotulabelscommand_h
#define Mothur_getotulabelscommand_h

//
//  getotuscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/21/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"


/**************************************************************************************************/

class GetOtusCommand : public Command {
public:
    GetOtusCommand(string);
    ~GetOtusCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.otus";          }
    string getCommandCategory()		{ return "OTU-Based Approaches";	} 
    
	string getHelpString();	
    string getOutputPattern(string);	
    string getCitation() { return "http://www.mothur.org/wiki/Get.otus"; }
    string getDescription()		{ return "Can be used with output from classify.otu, otu.association, or corr.axes to select specific otus."; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort;
    string  accnosfile, constaxonomyfile, otucorrfile, corraxesfile, listfile, sharedfile, label;
    vector<string> outputNames;
    unordered_set<string> labels;
    ListVector* list;
    
    int readClassifyOtu();
    int readOtuAssociation();
    int readCorrAxes();
    int readList();
    int readShared();
    int getListVector();
    SharedRAbundVectors* getShared();
};

/**************************************************************************************************/






#endif
