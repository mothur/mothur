/*
 *  GetlineCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlinecommand.h"


GetlineCommand::GetlineCommand(){
	try {
		globaldata = GlobalData::getInstance();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlineCommand class Function GetlineCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlineCommand class function GetlineCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

GetlineCommand::~GetlineCommand(){
}

//**********************************************************************************************************************

int GetlineCommand::execute(){
	try {
		filename = globaldata->inputFileName;
		ifstream in;
		openInputFile(filename, in);
		string label;
		int numBins = 0;
		int count = -1;
		int line = 1;
		while(in.good())
		{
			if(count > numBins)
				count = 0;
			if(count == 0)
			{
				cout << line << "\n";
				in >> numBins;
				line++;
			}
			in >> label;
			count++;
		}
		return 0;		
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlineCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlineCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}



