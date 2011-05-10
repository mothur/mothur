/*
 *  sensspeccommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/6/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "sensspeccommand.h"

//**********************************************************************************************************************
vector<string> SensSpecCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(plist);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none",false,false); parameters.push_back(pphylip);
		//CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none",false,false); parameters.push_back(pcolumn);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pcutoff("cutoff", "Number", "", "-1.00", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);
		CommandParameter phard("hard", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(phard);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SensSpecCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sens.spec command....\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
SensSpecCommand::SensSpecCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sensspec"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "SensSpecCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SensSpecCommand::SensSpecCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			string temp;
			
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["sensspec"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				//it = parameters.find("name");
				//user has given a template file
				//if(it != parameters.end()){ 
					//path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					//if (path == "") {	parameters["name"] = inputDir + it->second;		}
				//}
				
			}
			//check for required parameters
			listFile = validParameter.validFile(parameters, "list", true);
			if (listFile == "not found") { 		
				listFile = m->getListFile(); 
				if (listFile != "") { m->mothurOut("Using " + listFile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listFile == "not open") { abort = true; }	
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not found") { phylipfile = "";  }
			else if (phylipfile == "not open") { abort = true; }	
			else { distFile = phylipfile; format = "phylip";   }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not found") { columnfile = ""; }
			else if (columnfile == "not open") { abort = true; }	
			else { distFile = columnfile; format = "column";   }
			
			if ((phylipfile == "") && (columnfile == "")) { //is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  distFile = columnfile; format = "column";  m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") {  distFile = phylipfile; format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a sens.spec command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);
			if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(listFile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			temp = validParameter.validFile(parameters, "hard", false);
			if (temp == "not found"){	hard = 0;	}
			else if(!m->isTrue(temp))	{	hard = 0;	}
			else if(m->isTrue(temp))	{	hard = 1;	}
			
//			temp = validParameter.validFile(parameters, "name", true);
//			if (temp == "not found")	{	nameFile = "";		}
//			else if(temp == "not open")	{	abort = true;		}
//			else						{	nameFile = temp;	}
//			cout << "name:\t" << nameFile << endl;

			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found") { temp = "-1.00"; }
			convert(temp, cutoff);  
//			cout << cutoff << endl;
			
			temp = validParameter.validFile(parameters, "precision", false);	if (temp == "not found") { temp = "100"; }
			convert(temp, precision);  
//			cout << precision << endl;
			
			lineLabel = validParameter.validFile(parameters, "label", false);	if (lineLabel == "not found") { lineLabel = ""; }
			
			sensSpecFileName = outputDir + m->getRootName(m->getSimpleName(listFile)) + ".sensspec";
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "SensSpecCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SensSpecCommand::execute(){
	try{
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		setUpOutput();
		outputNames.push_back(sensSpecFileName); outputTypes["sensspec"].push_back(sensSpecFileName);
		if(format == "phylip")		{	processPhylip();	}
		else if(format == "column")	{	processColumn();	}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(sensSpecFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		
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
		//probably need some checking to confirm that the names in the distance matrix are the same as those in the list file
		
		ifstream inputListFile;
		m->openInputFile(listFile, inputListFile);
		
		string origCutoff = "";
		bool getCutoff = 0;
		if(cutoff == -1.00)	{	getCutoff = 1;															}
		else				{	origCutoff = toString(cutoff);	cutoff += (0.49 / double(precision));	}
		
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
			m->gobble(inputListFile);
		
			int lNumSeqs = seqMap.size();
			int pNumSeqs = 0;

			ifstream phylipFile;
			m->openInputFile(distFile, phylipFile);
			phylipFile >> pNumSeqs;
			if(pNumSeqs != lNumSeqs){	cout << "numSeq mismatch!" << endl;	}
			
			string seqName;
			double distance;
			vector<int> otuIndices(lNumSeqs, -1);
				
			truePositives = 0;
			falsePositives = 0;
			trueNegatives = 0;
			falseNegatives = 0;
			
			if(getCutoff == 1){
				if(label != "unique"){
					origCutoff = label;
					convert(label, cutoff);
					if(hard == 0){	cutoff += (0.49 / double(precision));	}		
				}
				else{
					origCutoff = "unique";
					cutoff = 0.0000;
				}
			}
				   
			cout << label << endl;
			
			for(int i=0;i<lNumSeqs;i++){
				phylipFile >> seqName;
				otuIndices[i] = seqMap[seqName];
				
				for(int j=0;j<i;j++){
					phylipFile >> distance;
					
					if(distance <= cutoff){
						if(otuIndices[i] == otuIndices[j])	{	truePositives++;	}
						else								{	falseNegatives++;	}
					}
					else{
						if(otuIndices[i] == otuIndices[j])	{	falsePositives++;	}
						else								{	trueNegatives++;	}
					}
				}
			}
			phylipFile.close();
			
			outputStatistics(label, origCutoff);
		}
		inputListFile.close();

	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "processPhylip");
		exit(1);
	}
}

//***************************************************************************************************************

void SensSpecCommand::processColumn(){
	try{		
		ifstream inputListFile;
		m->openInputFile(listFile, inputListFile);
		
		string origCutoff = "";
		bool getCutoff = 0;
		if(cutoff == -1.00)	{	getCutoff = 1;															}
		else				{	origCutoff = toString(cutoff);	cutoff += (0.49 / double(precision));	}
		
		set<string> seqPairSet;
		
		string label, seqList;
		int numOTUs;
		int numSeqs;
		
		while(inputListFile){
			numSeqs = 0;
			
			inputListFile >> label >> numOTUs;
			for(int i=0;i<numOTUs;i++){
				
				vector<string> seqNameVector;
				
				inputListFile >> seqList;
				int seqListLength = seqList.length();
				string seqName = "";
				for(int j=0;j<seqListLength;j++){
					
					if(seqList[j] == ','){
						seqNameVector.push_back(seqName);
						seqName = "";
					}
					else{
						seqName += seqList[j];
					}
				}
				seqNameVector.push_back(seqName);

				numSeqs += seqNameVector.size();
				
				int numSeqsInOTU = seqNameVector.size();
				for(int j=0;j<numSeqsInOTU;j++){
					string seqPairString = "";				
					for(int k=0;k<j;k++){
						if(seqNameVector[j] < seqNameVector[k])	{	seqPairString = seqNameVector[j] + '\t' + seqNameVector[k];	}
						else									{	seqPairString = seqNameVector[k] + '\t' + seqNameVector[j];	}
						seqPairSet.insert(seqPairString);
					}
				}
			}
			m->gobble(inputListFile);
			
			int numDists = (numSeqs * (numSeqs-1) / 2);

			ifstream columnFile;
			m->openInputFile(distFile, columnFile);
			string seqNameA, seqNameB, seqPairString;
			double distance;
			
			truePositives = 0;
			falsePositives = 0;
			trueNegatives = numDists;
			falseNegatives = 0;
			
			if(getCutoff == 1){
				if(label != "unique"){
					origCutoff = label;
					convert(label, cutoff);
					if(hard == 0){	cutoff += (0.49 / double(precision));	}		
				}
				else{
					origCutoff = "unique";
					cutoff = 0.0000;
				}
			}
			
			cout << label << endl;
			
			while(columnFile){
				columnFile >> seqNameA >> seqNameB >> distance;
				if(seqNameA < seqNameB)	{	seqPairString = seqNameA + '\t' + seqNameB;	}
				else					{	seqPairString = seqNameB + '\t' + seqNameA;	}

				set<string>::iterator it = seqPairSet.find(seqPairString);
			
				if(distance <= cutoff){
					if(it != seqPairSet.end()){
						truePositives++;
						seqPairSet.erase(it);	
					}
					else{
						falseNegatives++;
					}
					trueNegatives--;
				}
				else if(it != seqPairSet.end()){	
					falsePositives++;
					trueNegatives--;
					seqPairSet.erase(it);	
				}
				
				m->gobble(columnFile);
			}
			falsePositives += seqPairSet.size();
			
			outputStatistics(label, origCutoff);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "processColumn");
		exit(1);
	}
}

//***************************************************************************************************************

void SensSpecCommand::setUpOutput(){
	try{		
		ofstream sensSpecFile;
		m->openOutputFile(sensSpecFileName, sensSpecFile);
		
		sensSpecFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

		sensSpecFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "setUpOutput");
		exit(1);
	}
}

//***************************************************************************************************************

void SensSpecCommand::outputStatistics(string label, string cutoff){
	try{		
		double tp = (double) truePositives;
		double fp = (double) falsePositives;
		double tn = (double) trueNegatives;
		double fn = (double) falseNegatives;
		
		double p = tp + fn;
		double n = fp + tn;
		double pPrime = tp + fp;
		double nPrime = tn + fn;
		
		double sensitivity = tp / p;	
		double specificity = tn / n;
		double positivePredictiveValue = tp / pPrime;
		double negativePredictiveValue = tn / nPrime;
		double falseDiscoveryRate = fp / pPrime;
		
		double accuracy = (tp + tn) / (p + n);
		double matthewsCorrCoef = (tp * tn - fp * fn) / sqrt(p * n * pPrime * nPrime);	if(p == 0 || n == 0){	matthewsCorrCoef = 0;	}
		double f1Score = 2.0 * tp / (p + pPrime);
		
		
		if(p == 0)			{	sensitivity = 0;	matthewsCorrCoef = 0;	}
		if(n == 0)			{	specificity = 0;	matthewsCorrCoef = 0;	}
		if(p + n == 0)		{	accuracy = 0;								}
		if(p + pPrime == 0)	{	f1Score = 0;								}
		if(pPrime == 0)		{	positivePredictiveValue = 0;	falseDiscoveryRate = 0;	matthewsCorrCoef = 0;	}
		if(nPrime == 0)		{	negativePredictiveValue = 0;	matthewsCorrCoef = 0;							}
		
		ofstream sensSpecFile;
		m->openOutputFileAppend(sensSpecFileName, sensSpecFile);
		
		sensSpecFile << label << '\t' << cutoff << '\t';
		sensSpecFile << truePositives << '\t' << trueNegatives << '\t' << falsePositives << '\t' << falseNegatives << '\t';
		sensSpecFile << setprecision(4);
		sensSpecFile << sensitivity << '\t' << specificity << '\t' << positivePredictiveValue << '\t' << negativePredictiveValue << '\t';
		sensSpecFile << falseDiscoveryRate << '\t' << accuracy << '\t' << matthewsCorrCoef << '\t' << f1Score << endl;
		
		sensSpecFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "outputStatistics");
		exit(1);
	}
}

//***************************************************************************************************************



