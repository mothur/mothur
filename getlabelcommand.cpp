/*
 *  GetlabelCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlabelcommand.h"



GetlabelCommand::GetlabelCommand(){
	try {
		globaldata = GlobalData::getInstance();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlabelCommand class Function GetlabelCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlabelCommand class function GetlabelCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

GetlabelCommand::~GetlabelCommand(){
}

//**********************************************************************************************************************

int GetlabelCommand::execute(){
	try {
		filename = globaldata->inputFileName;
		ifstream in;
		openInputFile(filename, in);
		string label;
		int numBins = 0;
		int count = -1;
		while(in.good())
		{
			if(count > numBins)
				count = 0;
			if(count == 0)
			{
				cout << label << "\n";
				in >> numBins;
			}
			in >> label;
			count++;
		}	
		return 0;	
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlabelCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlabelCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

