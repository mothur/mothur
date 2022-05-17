//
//  removedistscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 1/29/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_removedistscommand_h
#define Mothur_removedistscommand_h

#include "command.hpp"

class RemoveDistsCommand : public Command {
	
public:
	
	RemoveDistsCommand(string);	
	~RemoveDistsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "remove.dists";			}
	string getCommandCategory()		{ return "General";                 }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Remove.dists"; }
	string getDescription()		{ return "removes distances from a phylip or column file related to groups or sequences listed in an accnos file"; }
    
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
    unordered_set<string> names;
	string accnosfile, phylipfile, columnfile;
	bool abort;
	vector<string> outputNames;
	
	int readPhylip();
	int readColumn();
	
};


#endif
