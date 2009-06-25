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
			if(fileList == "not found") { mothurOut("you must enter two or more file names"); mothurOutEndLine();  abort=true;  }
			else{ 	splitAtDash(fileList, fileNames);	}
			
			numInputFiles = fileNames.size();
			ifstream testFile;
			if(numInputFiles == 0){
				mothurOut("you must enter two or more file names and you entered " + toString(fileNames.size()) +  " file names"); mothurOutEndLine();
				abort=true;  
			}
			else{
				for(int i=0;i<numInputFiles;i++){
					if(openInputFile(fileNames[i], testFile)){	abort = true;	}
					testFile.close();
				}
			}   
			
			outputFileName = validParameter.validFile(parameters, "output", false);			
			if (outputFileName == "not found") { mothurOut("you must enter an output file name"); mothurOutEndLine();  abort=true;  }
		}
			
	}
	catch(exception& e) {
		errorOut(e, "MergeFileCommand", "MergeFileCommand");
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
		
		outputFile.close();
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MergeFileCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void MergeFileCommand::help(){
	try {
		mothurOut("The merge.file command..."); mothurOutEndLine();
	}
	catch(exception& e) {
		errorOut(e, "MergeFileCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************
