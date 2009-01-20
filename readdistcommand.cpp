/*
 *  readdistcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readdistcommand.h"

ReadDistCommand::ReadDistCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		filename = globaldata->inputFileName;
		format = globaldata->getFormat();	
		
		if (format == "column") { read = new ReadColumnMatrix(filename); }	
		else if (format == "phylip") { read = new ReadPhylipMatrix(filename); }	
			
		if(globaldata->getPrecision() != ""){
			convert(globaldata->getPrecision(), precision);	
		}
		
		if(globaldata->getCutOff() != ""){
			convert(globaldata->getCutOff(), cutoff);	
			cutoff += (5 / (precision * 10.0));
		}
		read->setCutoff(cutoff);
	
		if(globaldata->getNameFile() != ""){	
			nameMap = new NameAssignment(globaldata->getNameFile());
			nameMap->readMap(1,2);
		}
		else{
			nameMap = NULL;
		}
			
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistCommand class Function ReadDistCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistCommand class function ReadDistCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************
ReadDistCommand::~ReadDistCommand(){
	delete read;
	delete nameMap;
}

//**********************************************************************************************************************
int ReadDistCommand::execute(){
	try {
		read->read(nameMap);
		globaldata->setListVector(read->getListVector());
		globaldata->setSparseMatrix(read->getMatrix());
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}