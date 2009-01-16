/*
 *  readdistphylipfilecommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readdistphylipfilecommand.h"

//**********************************************************************************************************************

ReadDistPhylipFileCommand::ReadDistPhylipFileCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		filename = globaldata->inputFileName;
				
		format = globaldata->getFormat();	
		read = new ReadPhylipMatrix(filename);	
		
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
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistPhylipFileCommand class Function ReadDistPhylipFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistPhylipFileCommand class function ReadDistPhylipFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ReadDistPhylipFileCommand::~ReadDistPhylipFileCommand(){
	delete read;
	delete nameMap;
}

//**********************************************************************************************************************

int ReadDistPhylipFileCommand::execute(){
	try {
		read->read(nameMap);
		globaldata->setListVector(read->getListVector());
		globaldata->setSparseMatrix(read->getMatrix());
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistPhylipFileCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistPhylipFileCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
