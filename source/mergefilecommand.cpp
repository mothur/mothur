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
vector<string> MergeFileCommand::setParameters(){	
	try {
		CommandParameter pinput("input", "String", "", "", "", "", "",false,true); parameters.push_back(pinput);
		CommandParameter poutput("output", "String", "", "", "", "", "",false,true); parameters.push_back(poutput);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MergeFileCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The merge.file command takes a list of files separated by dashes and merges them into one file."; 
		helpString += "The merge.file command parameters are input and output."; 
		helpString += "Example merge.file(input=small.fasta-large.fasta, output=all.fasta).";
		helpString += "Note: No spaces between parameter labels (i.e. output), '=' and parameters (i.e.yourOutputFileName).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
MergeFileCommand::MergeFileCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["merge"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "MergeFileCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MergeFileCommand::MergeFileCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		if(option == "help") {
			help();
			abort = true; calledHelp = true;
		}else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["merge"] = tempOutNames;
			
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
			else if (outputDir != "") { outputFileName = outputDir + m->getSimpleName(outputFileName);  }
		}
			
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "MergeFileCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeFileCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		ofstream outputFile;
		m->openOutputFile(outputFileName, outputFile);
		
		char c;
		for(int i=0;i<numInputFiles;i++){
			ifstream inputFile; //declaration must be inside for loop of windows throws an error
			
			m->openInputFile(fileNames[i], inputFile);
			
			while(!inputFile.eof()){	
				if (m->control_pressed) { outputTypes.clear(); inputFile.close(); outputFile.close(); m->mothurRemove(outputFileName); return 0;  }
			
				c = inputFile.get(); 
				//-1 is eof char
				if (int(c) != -1) { outputFile << c; }   
			}
			
			inputFile.close();
		}
		
		outputFile.close();
		
		if (m->control_pressed) { outputTypes.clear();  m->mothurRemove(outputFileName); return 0;  }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	outputNames.push_back(outputFileName); outputTypes["merge"].push_back(outputFileName);
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeFileCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
