#ifndef Mothur_createdatabasecommand_h
#define Mothur_createdatabasecommand_h

//
//  createdatabasecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 3/28/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "command.hpp"
#include "listvector.hpp"
#include "sequence.hpp"

class CreateDatabaseCommand : public Command {
public:
	CreateDatabaseCommand(string);
	CreateDatabaseCommand();
	~CreateDatabaseCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "create.database";		}
	string getCommandCategory()		{ return "OTU-Based Approaches"; }
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Create.database"; }
	string getDescription()		{ return "creates database file that includes, abundances across groups, representative sequences, and taxonomy for each OTU"; }
    
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	
	bool abort;
	string listfile, groupfile, repfastafile, repnamesfile, contaxonomyfile, label, outputDir;
	
	vector<string> outputNames;
		
	vector<int> readFasta(vector<Sequence>&);
    vector<int> readTax(vector<string>&);
    int readNames(map<string, string>&); 
	ListVector* getList();
	
};




#endif
