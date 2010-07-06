/*
 *  sensspeccommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/6/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "sensspeccommand.h"

//***************************************************************************************************************

SensSpecCommand::SensSpecCommand(string option)  {
	try {
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			string temp;

			//valid paramters for this command
			string AlignArray[] =  {"list", "phylip", "column", "name", "hard", "label", "cutoff", "outputdir", "inputdir"};
			
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
			}
			//check for required parameters
			listFile = validParameter.validFile(parameters, "list", true);
			if (listFile == "not found") { m->mothurOut("list is a required parameter for the sens.spec command."); m->mothurOutEndLine(); abort = true; }
			else if (listFile == "not open") { abort = true; }	
			
			distFile = validParameter.validFile(parameters, "column", true);
			format = "column";
			if(distFile == "not found")	{
				distFile = validParameter.validFile(parameters, "phylip", true);
				format = "phylip";	
			}
			if(distFile == "not found") { m->mothurOut("either column or phylip are required for the sens.spec command."); m->mothurOutEndLine(); abort = true; }
			else if (distFile == "not open") { abort = true; }	
		
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);
			if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(listFile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			temp = validParameter.validFile(parameters, "hard", false);
			if (temp == "not found"){	hard = 0;	}
			else if(!isTrue(temp))	{	hard = 0;	}
			else if(isTrue(temp))	{	hard = 1;	}
			
//			temp = validParameter.validFile(parameters, "name", true);
//			if (temp == "not found")	{	nameFile = "";		}
//			else if(temp == "not open")	{	abort = true;		}
//			else						{	nameFile = temp;	}
//			cout << "name:\t" << nameFile << endl;

			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found") { temp = ""; }
			convert(temp, cutoff);  
			cout << "cutoff:\t" << cutoff << endl;
			
//			cutoff = 0.0349;
			
			lineLabel = validParameter.validFile(parameters, "label", false);	if (lineLabel == "not found") { lineLabel = ""; }
			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "SensSpecCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SensSpecCommand::help(){
	try {
		m->mothurOut("The sens.spec command reads a fastaFile and creates .....\n");

		
		
		m->mothurOut("Example sens.spec(...).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n");
		m->mothurOut("For more details please check out the wiki http://www.mothur.org/wiki/Trim.seqs .\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "help");
		exit(1);
	}
}


//***************************************************************************************************************

SensSpecCommand::~SensSpecCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int SensSpecCommand::execute(){
	try{
		if (abort == true) { return 0; }

		if(format == "phylip")		{	processPhylip();	}
//		else if(format == "column")	{	processColumn(seqMap);	}
		
		
//		string seqList;
//		map<string, int> seqMap;
		
		
		
		return 0;	
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************

void SensSpecCommand::processPhylip(){
	try{
		ifstream inputListFile;
		openInputFile(listFile, inputListFile);

		string label;
		int numOTUs;

		map<string, int> seqMap;
		string seqList;
		
		while(inputListFile){
			inputListFile >> label >> numOTUs;
			for(int i=0;i<numOTUs;i++){
				inputListFile >> seqList;
				int seqListLength = seqList.length();
				string seqName = "";
				for(int j=0;j<seqListLength;j++){
					
					if(seqList[j] == ','){
						seqMap[seqName] = i;
						seqName = "";
					}
					else{
						seqName += seqList[j];
					}
					
				}
				seqMap[seqName] = i;
			}
		}
		inputListFile.close();
		
		int lNumSeqs = seqMap.size();
		int pNumSeqs = 0;

		ifstream phylipFile;
		openInputFile(distFile, phylipFile);
		phylipFile >> pNumSeqs;
		if(pNumSeqs != lNumSeqs){	cout << "numSeq mismatch!" << endl;	}
		
		string seqName;
		double distance;
		vector<int> otuIndices(lNumSeqs, -1);
			
		truePositives = 0;
		falsePositives = 0;
		trueNegatives = 0;
		falseNegatives = 0;
		
		
		for(int i=0;i<lNumSeqs;i++){
			phylipFile >> seqName;
			otuIndices[i] = seqMap[seqName];
			
			for(int j=0;j<i;j++){
				phylipFile >> distance;
				if(distance <= cutoff){
					if(otuIndices[i] == otuIndices[j]){
						truePositives++;
					}
					else{
						falseNegatives++;
					}
				}
				else{
					if(otuIndices[i] == otuIndices[j]){
						falsePositives++;
					}
					else{
						trueNegatives++;
					}
				}
			}
		}
		phylipFile.close();
		
		cout << "truePositives:\t"	<< truePositives	<< endl;
		cout << "trueNegatives:\t"	<< trueNegatives	<< endl;
		cout << "falsePositives:\t"	<< falsePositives	<< endl;
		cout << "falseNegatives:\t"	<< falseNegatives	<< endl;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "processPhylip");
		exit(1);
	}
}

//***************************************************************************************************************

//void SensSpecCommand::processColumn(map<string, int> seqMap){
//	
//	truePositives = 0;
//	falsePositives = 0;
//	trueNegatives = numDists;
//	falseNegatives = 0;
//	
//	ifstream columnFile;
//	openInputFile(distFile, columnFile);
//	
//	string seqNameA, seqNameB, oldSeqNameA;
//	int otuA, otuB, oldOTUA;
//	double distance;
//
//	while(columnFile){
//		columnFile >> seqNameA >> seqNameB >> distance;
//
//		if(seqNameA == oldSeqNameA)	{	otuA = oldOTUA;	}
//		else						{	otuA = seqMap[seqNameA];	oldOTUA = otuA;	}
//		
//		otuB = seqMap[seqNameB];
//		
//		if(distance <= cutoff){
//			if(otuA == otuB){
//				truePositives++;
//			}
//			else{
//				falseNegatives++;
//			}
//			trueNegatives--;
//		}
//		else{
//			if(otuA == otuB){
//				falsePositives++;
//				trueNegatives--;
//			}
//		}
//	
//		gobble(columnFile);
//	}
//	columnFile.close();
//	
//	cout << "truePositives:\t"	<< truePositives	<< endl;
//	cout << "trueNegatives:\t"	<< trueNegatives	<< endl;
//	cout << "falsePositives:\t"	<< falsePositives	<< endl;
//	cout << "falseNegatives:\t"	<< falseNegatives	<< endl;
//}

//***************************************************************************************************************
