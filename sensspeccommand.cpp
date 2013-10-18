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
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","sensspec",false,true,true); parameters.push_back(plist);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none","",false,false); parameters.push_back(pphylip);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none","",false,false); parameters.push_back(pcolumn);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pcutoff("cutoff", "Number", "", "-1.00", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(phard);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
string SensSpecCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "sensspec") {  pattern = "[filename],sensspec"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SensSpecCommand", "getOutputPattern");
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
		allLines = 1;
		
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
			}
			//check for required parameters
			listFile = validParameter.validFile(parameters, "list", true);
			if (listFile == "not found") { 		
				listFile = m->getListFile(); 
				if (listFile != "") { m->mothurOut("Using " + listFile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listFile == "not open") { abort = true; }	
			else { m->setListFile(listFile); }
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not found") { phylipfile = "";  }
			else if (phylipfile == "not open") { abort = true; }	
			else { distFile = phylipfile; format = "phylip"; m->setPhylipFile(phylipfile);  }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not found") { columnfile = ""; }
			else if (columnfile == "not open") { abort = true; }	
			else { distFile = columnfile; format = "column";   m->setColumnFile(columnfile); }
			
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
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found") { temp = "-1.00"; }
			m->mothurConvert(temp, cutoff);  
//			cout << cutoff << endl;
			
			temp = validParameter.validFile(parameters, "precision", false);	if (temp == "not found") { temp = "100"; }
			m->mothurConvert(temp, precision);  
//			cout << precision << endl;
			
			string label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
            map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listFile));
			sensSpecFileName = getOutputFileName("sensspec",variables);
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
        
        int startTime = time(NULL);
        
        //create list file with only unique names, saves time and memory by removing redundant names from list file that are not in the distance file.
        string newListFile = preProcessList();
        if (newListFile != "") { listFile = newListFile; }
        
		setUpOutput();
		outputNames.push_back(sensSpecFileName); outputTypes["sensspec"].push_back(sensSpecFileName);
		if(format == "phylip")		{	processPhylip();	}
		else if(format == "column")	{	processColumn();	}
		
        //remove temp file if created
        if (newListFile != "") { m->mothurRemove(newListFile); }
        
		if (m->control_pressed) { m->mothurRemove(sensSpecFileName); return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - startTime) + " to run sens.spec."); m->mothurOutEndLine();
        
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
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
bool SensSpecCommand::testFile(){
	try{
        ifstream fileHandle;
        m->openInputFile(phylipfile, fileHandle);
        
        bool square = false;
        string numTest, name;
        fileHandle >> numTest >> name;

        if (!m->isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }

        char d;
        while((d=fileHandle.get()) != EOF){
            if(isalnum(d)){
                square = true;
                break;
            }
            if(d == '\n'){
                square = false;
                break;
            }
        }
        fileHandle.close();
        
        return square;
    }
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "testFile");
		exit(1);
	}
}

//***************************************************************************************************************

int SensSpecCommand::processPhylip(){
	try{
		//probably need some checking to confirm that the names in the distance matrix are the same as those in the list file
        square = testFile();
 		string origCutoff = "";
		bool getCutoff = 0;
		if(cutoff == -1.00)	{	getCutoff = 1;															}
		else				{	origCutoff = toString(cutoff);	cutoff += (0.49 / double(precision));	}
		
		map<string, int> seqMap;
		string seqList;
		
		InputData input(listFile, "list");
		ListVector* list = input.getListVector();
		string lastLabel = list->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(m->control_pressed){
                for (int i = 0; i < outputNames.size(); i++){	m->mothurRemove(outputNames[i]);  }  delete list;  return 0;
            }
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){			
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//process
				fillSeqMap(seqMap, list);
				process(seqMap, list->getLabel(), getCutoff, origCutoff);
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input.getListVector(lastLabel);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//process
				fillSeqMap(seqMap, list);
				process(seqMap, list->getLabel(), getCutoff, origCutoff);
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}		
			
			lastLabel = list->getLabel();			
			
			delete list;
			list = input.getListVector();
		}
		
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input.getListVector(lastLabel);
			
			//process
			fillSeqMap(seqMap, list);
			process(seqMap, list->getLabel(), getCutoff, origCutoff);
			
			delete list;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "processPhylip");
		exit(1);
	}
}

//***************************************************************************************************************

int SensSpecCommand::fillSeqMap(map<string, int>& seqMap, ListVector*& list){
	try {
		//for each otu
		for(int i=0;i<list->getNumBins();i++){
			
			if (m->control_pressed) { return 0; }
			
			string seqList = list->get(i);
			int seqListLength = seqList.length();
			string seqName = "";
			
			//parse bin by name, mapping each name to its otu number
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
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "fillSeqMap");
		exit(1);
	}
}
//***************************************************************************************************************
int SensSpecCommand::fillSeqPairSet(set<string>& seqPairSet, ListVector*& list){
	try {
		int numSeqs = 0;
		
		//for each otu
		for(int i=0;i<list->getNumBins();i++){
			
			if (m->control_pressed) { return 0; }
			
			vector<string> seqNameVector;
			string bin = list->get(i);
			m->splitAtComma(bin, seqNameVector);
			
			numSeqs += seqNameVector.size();
			
			for(int j=0;j<seqNameVector.size();j++){
				string seqPairString = "";				
				for(int k=0;k<j;k++){
					if(seqNameVector[j] < seqNameVector[k])	{	seqPairString = seqNameVector[j] + '\t' + seqNameVector[k];	}
					else									{	seqPairString = seqNameVector[k] + '\t' + seqNameVector[j];	}
					seqPairSet.insert(seqPairString);
				}
			}
		}
		
		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "fillSeqPairSet");
		exit(1);
	}
}
//***************************************************************************************************************
int SensSpecCommand::process(map<string, int>& seqMap, string label, bool& getCutoff, string& origCutoff){
	try {
						
		int lNumSeqs = seqMap.size();
		int pNumSeqs = 0;
		
		ifstream phylipFile;
		m->openInputFile(distFile, phylipFile);
		phylipFile >> pNumSeqs;
		if(pNumSeqs != lNumSeqs){	m->mothurOut("numSeq mismatch!\n"); /*m->control_pressed = true;*/ }
		
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
		
		m->mothurOut(label); m->mothurOutEndLine();
		
		for(int i=0;i<pNumSeqs;i++){
			
			if (m->control_pressed) { return 0; }
			
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
            
            if (square) { m->getline(phylipFile); } //get rest of line - redundant distances
            m->gobble(phylipFile);
		}
		phylipFile.close();
		
		outputStatistics(label, origCutoff);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "process");
		exit(1);
	}
}
//***************************************************************************************************************
int SensSpecCommand::process(set<string>& seqPairSet, string label, bool& getCutoff, string& origCutoff, int numSeqs){
	try {
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
		
		m->mothurOut(label); m->mothurOutEndLine();
		
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
		
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "process");
		exit(1);
	}
}
//***************************************************************************************************************

int SensSpecCommand::processColumn(){
	try{		
		string origCutoff = "";
		bool getCutoff = 0;
		if(cutoff == -1.00)	{	getCutoff = 1;															}
		else				{	origCutoff = toString(cutoff);	cutoff += (0.49 / double(precision));	}
		
		set<string> seqPairSet;
		int numSeqs = 0;
		
		InputData input(listFile, "list");
		ListVector* list = input.getListVector();
		string lastLabel = list->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }  delete list;  return 0;  }
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){			
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//process
				numSeqs = fillSeqPairSet(seqPairSet, list);
				process(seqPairSet, list->getLabel(), getCutoff, origCutoff, numSeqs);
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input.getListVector(lastLabel);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//process
				numSeqs = fillSeqPairSet(seqPairSet, list);
				process(seqPairSet, list->getLabel(), getCutoff, origCutoff, numSeqs);
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}		
			
			lastLabel = list->getLabel();			
			
			delete list;
			list = input.getListVector();
		}
		
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input.getListVector(lastLabel);
			
			//process
			numSeqs = fillSeqPairSet(seqPairSet, list);
			delete list;
			process(seqPairSet, list->getLabel(), getCutoff, origCutoff, numSeqs);
		}
		
		return 0;
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

string SensSpecCommand::preProcessList(){
    try {
        set<string> uniqueNames;
        //get unique names from distance file
        if (format == "phylip") {
            
            ifstream phylipFile;
            m->openInputFile(distFile, phylipFile);
            string numTest;
            int pNumSeqs;
			phylipFile >> numTest; m->gobble(phylipFile);
			
			if (!m->isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
            else {
                m->mothurConvert(numTest, pNumSeqs);
            }
            
            string seqName;
            for(int i=0;i<pNumSeqs;i++){
                if (m->control_pressed) { return ""; }
                phylipFile >> seqName;  m->getline(phylipFile);  m->gobble(phylipFile);
                uniqueNames.insert(seqName);
            }
            phylipFile.close();
        }else {
            ifstream columnFile;
            m->openInputFile(distFile, columnFile);
            string seqNameA, seqNameB;
            double distance;
            
            while(columnFile){
                if (m->control_pressed) { return ""; }
                columnFile >> seqNameA >> seqNameB >> distance;
                uniqueNames.insert(seqNameA); uniqueNames.insert(seqNameB);
                m->gobble(columnFile);
            }
            columnFile.close();
        }
        
        //read list file, if numSeqs > unique names then remove redundant names
        string newListFile = listFile + ".temp";
        ofstream out;
        m->openOutputFile(newListFile, out);
        ifstream in;
		m->openInputFile(listFile, in);
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(newListFile);  return ""; }
            
			//read in list vector
			ListVector list(in);
            
            //listfile is already unique
            if (list.getNumSeqs() == uniqueNames.size()) { in.close(); out.close(); m->mothurRemove(newListFile);  return ""; }
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list.getLabel());
            vector<string> binLabels = list.getLabels();
            vector<string> newLabels;
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
                
				//parse out names that are in accnos file
				string binnames = list.get(i);
                vector<string> bnames;
                m->splitAtComma(binnames, bnames);
				
				string newNames = "";
                for (int j = 0; j < bnames.size(); j++) {
					string name = bnames[j];
					//if that name is in the .accnos file, add it
					if (uniqueNames.count(name) != 0) {  newNames += name + ",";  }
				}
                
				//if there are names in this bin add to new list
				if (newNames != "") { 
					newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
					newList.push_back(newNames);
                    newLabels.push_back(binLabels[i]);
				}
			}
            
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
                newList.setLabels(newLabels);
                if (!m->printedListHeaders) { newList.printHeaders(out); }
				newList.print(out);
			}
			
			m->gobble(in);
		}
		in.close();	
		out.close();

        if (wroteSomething) { return newListFile; }
        else { m->mothurRemove(newListFile); }
        
        return ""; 
    }
    catch(exception& e) {
        m->errorOut(e, "SensSpecCommand", "preProcessList");
        exit(1);
    }
}


//***************************************************************************************************************



