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
		
		//open shared file
		sharedfile = globaldata->getSharedFile();
		openInputFile(sharedfile, in);
		
		//open output file
		outputFile = getRootName(globaldata->inputFileName) + "bootGroups";
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
		int num, inputData, count;
		count = 0;  
		string holdLabel, nextLabel, groupN, label;
		
		//read in first row since you know there is at least 1 group.
		in >> label >> groupN >> num;
		holdLabel = label;
		
		//output first group
		cout << groupN << endl;
		out << groupN << '\t' << groupN << endl;	
			
		//get rest of line
		for(int i=0;i<num;i++){
			in >> inputData;
		}
		
		if (in.eof() != true) { in >> nextLabel; }
		
		//read the rest of the groups info in
		while ((nextLabel == holdLabel) && (in.eof() != true)) {
			in >> groupN >> num;
			count++;
			
			//output next group
			cout << groupN << endl;
			out << groupN << '\t' << groupN << endl;				
			
			//fill vector.  
			for(int i=0;i<num;i++){
				in >> inputData;
			}
			
			if (in.eof() != true) { in >> nextLabel; }
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


