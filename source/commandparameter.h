#ifndef COMMANDPARAMETER_H
#define COMMANDPARAMETER_H


/*
 *  commandparameter.h
 *  Mothur
 *
 *  Created by westcott on 3/23/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */



#include "mothur.h"

//**********************************************************************************************************************

class CommandParameter {
	
	public:
    CommandParameter() { name = ""; type = ""; options = ""; optionsDefault = ""; chooseOnlyOneGroup = ""; chooseAtLeastOneGroup = ""; linkedGroup = ""; multipleSelectionAllowed = false; required = false; important = false; outputTypes = ""; }
		CommandParameter(string n, string t, string o, string d, string only, string atLeast, string linked, string opt, bool m, bool r, bool i) : name(n), type(t), options(o), optionsDefault(d), 
				chooseOnlyOneGroup(only), chooseAtLeastOneGroup(atLeast), linkedGroup(linked), outputTypes(opt),multipleSelectionAllowed(m), required(r), important(i) {}
        CommandParameter(string n, string t, string o, string d, string only, string atLeast, string linked, string opt, bool m, bool r) : name(n), type(t), options(o), optionsDefault(d), 
    chooseOnlyOneGroup(only), chooseAtLeastOneGroup(atLeast), linkedGroup(linked), outputTypes(opt), multipleSelectionAllowed(m), required(r)  { important = false; }
		~CommandParameter() {}
	
		string name;		//something like fasta, processors, method
		string type;  //must be set to "Boolean", "Multiple", "Number", "String", "InputTypes" - InputTypes is for file inputs
		string options; //if the parameter has specific options allowed, used for parameters of type "Multiple", something like "furthest-nearest-average", or "sobs-chao...", leave blank for command that do not required specific options
		string optionsDefault;   //the default for this parameter, could be something like "F" for a boolean or "100" for a number or "sobs-chao" for multiple
		
	
		//for chooseOnlyOneGroup, chooseAtLeastOneGroup and linkedGroup if no group is needed set to "none".
		string chooseOnlyOneGroup; //for file inputs: if a command has several options for input files but you can only choose one then put them in a group
									//for instance in the read.dist command you can use a phylip or column file but not both so set chooseOnlyOneGroup for both parameters to something like "DistanceFileGroup"
		string chooseAtLeastOneGroup; //for file inputs: if a command has several options for input files and you want to make sure one is choosen then put them in a group
									//for instance in the read.dist command you must provide a phylip or column file so set chooseAtLeastOneGroup for both parameters to something like "DistanceFileGroup"
		string linkedGroup; //for file inputs: if a command has a file option were if you provide one you must provide another you can put them in a group
										//for instance in the cluster command if you provide a column file you must provide a name file so set linkedGroup for both parameters to something like "ColumnNameGroup"

		bool multipleSelectionAllowed; //for "Multiple" type to say whether you can select multiple options, for instance for calc parameter set to true, but for method set to false
		bool required; //is this parameter required
    
        bool important; //is this parameter important.  The gui will put "important" parameters first in the option panel.
        
        string outputTypes; //types on files created by the command if this parameter is given.  ie. get.seqs command fasta parameter makes a fasta file. can be multiple values split by dashes.
		
	private:
};

//**********************************************************************************************************************

#endif
