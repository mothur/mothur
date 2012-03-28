//
//  createdatabasecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/28/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "createdatabasecommand.h"
#include "sequence.hpp"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> CreateDatabaseCommand::setParameters(){	
	try {
		CommandParameter pfasta("repfasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("repname", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pname);
		CommandParameter pcontaxonomy("contaxonomy", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pcontaxonomy);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(plist);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CreateDatabaseCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The create.database command reads a listfile, *.cons.taxonomy, *.rep.fasta, *.rep.names and optional groupfile, and creates a database file.\n";
		helpString += "The create.database command parameters are repfasta, list, repname, contaxonomy, group and label. List, repfasta, repnames, and contaxonomy are required.  \n";
        
        helpString += "The create.database command should be in the following format: \n";
		helpString += "create.database(repfasta=yourFastaFileFromGetOTURep, repname=yourNameFileFromGetOTURep, contaxonomy=yourConTaxFileFromClassifyOTU, list=yourListFile) \n";	
		helpString += "Example: create.database(repfasta=final.an.0.03.rep.fasta, name=final.an.0.03.rep.names, list=fina.an.list, label=0.03, contaxonomy=final.an.0.03.cons.taxonomy) \n";
		helpString += "Note: No spaces between parameter labels (i.e. repfasta), '=' and parameters (i.e.yourFastaFileFromGetOTURep).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************


