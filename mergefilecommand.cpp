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

MergeFileCommand::MergeFileCommand(string option)  {
	try {
		abort = false;
		
		if(option == "help") {
			help();
			abort = true; 
		}
		else {
			//valid paramters for this command
			string Array[] =  {"input", "output","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			
			string fileList = validParameter.validFile(parameters, "input", false);			
			if(fileList == "not found") { m->mothurOut("you must enter two or more file names"); m->mothurOutEndLine();  abort=true;  }
			else{ 	m->splitAtDash(fileList, fileNames);	}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			string outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found")	{	outputDir = "";		}
			
			
			numInputFiles = fileNames.size();
			ifstream testFile;
			if(numInputFiles == 0){
				m->mothurOut("you must enter two or more file names and you entered " + toString(fileNames.size()) +  " file names"); m->mothurOutEndLine();
				abort=true;  
			}
			else{
				for(int i=0;i<numInputFiles;i++){
					if (inputDir != "") {
						string path = m->hasPath(fileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fileNames[i] = inputDir + fileNames[i];		}
					}
					
					if(m->openInputFile(fileNames[i], testFile)){	abort = true;	}
					testFile.close();
				}
			}   
			
			outputFileName = validParameter.validFile(parameters, "output", false);			
			if (outputFileName == "not found") { m->mothurOut("you must enter an output file name"); m->mothurOutEndLine();  abort=true;  }
			else if (outputDir != "") { outputFileName = outputDir + m->getSimpleName(outputFileName); }
		}
			
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "MergeFileCommand");
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
		m->openOutputFile(outputFileName, outputFile);
		
		char c;
		for(int i=0;i<numInputFiles;i++){
			ifstream inputFile; //declaration must be inside for loop of windows throws an error
			
			m->openInputFile(fileNames[i], inputFile);
			
			while(!inputFile.eof()){	
				if (m->control_pressed) { inputFile.close(); outputFile.close(); remove(outputFileName.c_str()); return 0;  }
			
				c = inputFile.get(); 
				//-1 is eof char
				if (int(c) != -1) { outputFile << c; }   
			}
			
			inputFile.close();
		}
		
		outputFile.close();
		
		if (m->control_pressed) { remove(outputFileName.c_str()); return 0;  }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void MergeFileCommand::help(){
	try {
		m->mothurOut("The merge.file command..."); m->mothurOutEndLine();
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************
