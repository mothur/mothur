/*
 *  mergefilecommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/14/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mergefilecommand.h"

//**********************************************************************************************************************

MergeFileCommand::MergeFileCommand(string option){
	try {
		abort = false;
		
		if(option == "help") {
			help();
			abort = true; 
		}
		else {
			//valid paramters for this command
			string Array[] =  {"input", "output"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			string fileList = validParameter.validFile(parameters, "input", false);			
			if(fileList == "not found") { cout << "you must enter two or more file names" << endl;  abort=true;  }
			else{ 	splitAtDash(fileList, fileNames);	}
			
			numInputFiles = fileNames.size();
			ifstream testFile;
			if(numInputFiles == 0){
				cout << "you must enter two or more file names and you entered " << fileNames.size() <<  " file names" << endl;
				abort=true;  
			}
			else{
				for(int i=0;i<numInputFiles;i++){
					if(openInputFile(fileNames[i], testFile)){	abort = true;	}
					testFile.close();
				}
			}   
			
			outputFileName = validParameter.validFile(parameters, "output", false);			
			if (outputFileName == "not found") { cout << "you must enter an output file name" << endl;  abort=true;  }
		}
			
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the MergeFileCommand class Function MergeFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the MergeFileCommand class function MergeFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

MergeFileCommand::~MergeFileCommand()	{	/*	do nothing	*/	}

//**********************************************************************************************************************

int MergeFileCommand::execute(){
	try {
		if (abort == true) {	return 0;	}
		
		ofstream outputFile;
		openOutputFile(outputFileName, outputFile);
		
		ifstream inputFile;
		char c;
		for(int i=0;i<numInputFiles;i++){
			openInputFile(fileNames[i], inputFile);
			
			while(!inputFile.eof()){	c = inputFile.get(); outputFile << c;	}
			
			inputFile.close();
		}
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function BinSeqCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function BinSeqCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
	
}

//**********************************************************************************************************************

void MergeFileCommand::help(){
	try {
		cout << "The merge.file command..." << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the MergeFileCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the MergeFileCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
	
}

//**********************************************************************************************************************
