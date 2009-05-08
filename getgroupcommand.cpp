/*
 *  getgroupcommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getgroupcommand.h"


GetgroupCommand::GetgroupCommand(){
	try {
		globaldata = GlobalData::getInstance();
		groupMap = globaldata->gGroupmap;
		outputFile = globaldata->inputFileName + ".bootGroups";
		openOutputFile(outputFile, out);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetgroupCommand class Function GetgroupCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetgroupCommand class function GetgroupCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

GetgroupCommand::~GetgroupCommand(){
}

//**********************************************************************************************************************

int GetgroupCommand::execute(){
	try {
		vector<string> groupNames = groupMap->namesOfGroups;	
		for(int i = 0; i < groupNames.size(); i++) {
			cout << groupNames[i] << "\n";
			out << groupNames[i] << "\t" << groupNames[i] << "\n";
		}
		out.close();
		return 0;	
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetgroupCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetgroupCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


