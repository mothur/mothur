/*
 *  parsesffcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 2/6/10.
 *  Copyright 2010 Patrick D. Schloss. All rights reserved.
 *
 */

#include "parsesffcommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> ParseSFFCommand::getValidParameters(){	
	try {
		string Array[] =  {"sff", "oligos", "minlength", "outputdir", "inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ParseSFFCommand::ParseSFFCommand(){	
	try {
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["flow"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "ParseSFFCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ParseSFFCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"sff"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ParseSFFCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

ParseSFFCommand::ParseSFFCommand(string option){
	try {
		abort = false;
		
		if(option == "help") {
			help();
			abort = true; 
		}
		else {
			//valid paramters for this command
			string Array[] =  {"sff", "oligos", "minlength", "outputdir", "inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;

			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["flow"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("sff");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sff"] = inputDir + it->second;		}
				}
				
				it = parameters.find("oligos");
				//user has given an oligos file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			sffFile = validParameter.validFile(parameters, "sff", true);
			if (sffFile == "not found"){
				m->mothurOut("sff is a required parameter for the parse.sff command.");
				m->mothurOutEndLine();
				abort = true;
			}
			else if (sffFile == "not open")		{	abort = true;	}	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);
			if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(sffFile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...			
			oligoFile = validParameter.validFile(parameters, "oligos", true);
			if (oligoFile == "not found")	{	oligoFile = "";		}
			else if(oligoFile == "not open"){	abort = true;		} 
			
			string temp = validParameter.validFile(parameters, "minlength", false);
			if (temp == "not found") { temp = "0"; }
			convert(temp, minLength); 
		}		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "ParseSFFCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

ParseSFFCommand::~ParseSFFCommand()	{	/*	do nothing	*/	}

//**********************************************************************************************************************

int ParseSFFCommand::execute(){
	try {
		if (abort == true) {	return 0;	}

		ifstream inSFF;
		m->openInputFile(sffFile, inSFF);
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);
		cout << setprecision(2);
			
		vector<ofstream*> flowFileNames;
		if(oligoFile != ""){
			getOligos(flowFileNames);
		}
		else{
			flowFileNames.push_back(new ofstream((outputDir + m->getRootName(m->getSimpleName(sffFile)) + "flow").c_str(), ios::ate));
			outputNames.push_back((outputDir + m->getRootName(m->getSimpleName(sffFile)) + "flow")); outputTypes["flow"].push_back((outputDir + m->getRootName(m->getSimpleName(sffFile)) + "flow"));
		}
		
		for(int i=0;i<flowFileNames.size();i++){
			flowFileNames[i]->setf(ios::fixed, ios::floatfield);
			flowFileNames[i]->setf(ios::showpoint);
			*flowFileNames[i] << setprecision(2);
		}			
		
		if (m->control_pressed) { outputTypes.clear(); for(int i=0;i<flowFileNames.size();i++){	flowFileNames[i]->close();  } return 0; }
		
//		ofstream fastaFile;
//		m->openOutputFile(m->getRootName(sffFile) + "fasta", fastaFile);

//		ofstream qualFile;
//		m->openOutputFile(m->getRootName(sffFile) + "qual", qualFile);
		
		string commonHeader = m->getline(inSFF);
		string magicNumber = m->getline(inSFF);		
		string version = m->getline(inSFF);
		string indexOffset = m->getline(inSFF);
		string indexLength = m->getline(inSFF);
		int numReads = parseHeaderLineToInt(inSFF);
		string headerLength = m->getline(inSFF);
		string keyLength = m->getline(inSFF);
		int numFlows = parseHeaderLineToInt(inSFF);
		string flowgramCode = m->getline(inSFF);
		string flowChars = m->getline(inSFF);
		string keySequence = m->getline(inSFF);
		m->gobble(inSFF);

		string seqName;
		bool good = 0;
		
		for(int i=0;i<numReads;i++){
			
			if (m->control_pressed) { outputTypes.clear(); for(int i=0;i<flowFileNames.size();i++){	flowFileNames[i]->close();  } return 0; }
			
			inSFF >> seqName;
			seqName = seqName.substr(1);
			m->gobble(inSFF);
			
			string runPrefix = parseHeaderLineToString(inSFF);
			string regionNumber = parseHeaderLineToString(inSFF);
			string xyLocation = parseHeaderLineToString(inSFF);
			m->gobble(inSFF);
			
			string runName = parseHeaderLineToString(inSFF);
			string analysisName = parseHeaderLineToString(inSFF);
			string fullPath = parseHeaderLineToString(inSFF);
			m->gobble(inSFF);
			
			string readHeaderLen = parseHeaderLineToString(inSFF);
			string nameLength = parseHeaderLineToString(inSFF);
			int numBases = parseHeaderLineToInt(inSFF);
			string clipQualLeft = parseHeaderLineToString(inSFF);
			int clipQualRight = parseHeaderLineToInt(inSFF);
			string clipAdapLeft = parseHeaderLineToString(inSFF);
			string clipAdapRight = parseHeaderLineToString(inSFF);
			m->gobble(inSFF);
			
			vector<float> flowVector = parseHeaderLineToFloatVector(inSFF, numFlows);
			vector<int> flowIndices = parseHeaderLineToIntVector(inSFF, numBases);
			string bases = parseHeaderLineToString(inSFF);
			string qualityScores = parseHeaderLineToString(inSFF);
			m->gobble(inSFF);
			

			
			int flowLength = flowIndices[clipQualRight-1];
						
			screenFlow(flowVector, flowLength);
			string sequence = flow2seq(flowVector, flowLength);
			
			int group = 0;
	
			if(minLength != 0 || numFPrimers != 0  || numBarcodes != 0 || numRPrimers != 0){		
				good = screenSeq(sequence, group);
			}

			if(good){
				*flowFileNames[group] << seqName << ' ' << flowLength;
				for(int i=0;i<numFlows;i++){
					*flowFileNames[group] << ' ' << flowVector[i];
				}
				*flowFileNames[group] << endl;				
			}
			
//			string fastaHeader = '>' + seqName + "\tregion=" + regionNumber + " xy=" + xyLocation;
//			fastaFile << fastaHeader << endl;
//			fastaFile << stripSeqQual(bases, clipQualLeft, clipQualRight) << endl;
//
//			qualFile << fastaHeader << endl;
//			qualFile << stripQualQual(qualityScores, clipQualLeft, clipQualRight) << endl;

		}
		for(int i=0;i<flowFileNames.size();i++){
			flowFileNames[i]->close();
		}

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

//		fastaFile.close();
//		qualFile.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void ParseSFFCommand::help(){
	try {
		m->mothurOut("The parse.sff command...");
		m->mothurOutEndLine();
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

void ParseSFFCommand::getOligos(vector<ofstream*>& outSFFFlowVec){
	try {

		ifstream inOligos;
		m->openInputFile(oligoFile, inOligos);
		
		string type, oligo, group;
		
		int index = 0;
		
		while(!inOligos.eof()){
			inOligos >> type;

			if(type[0] == '#'){	m->getline(inOligos);	} // get rest of line if there's any crap there
			else{
				inOligos >> oligo;
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				if(type == "forward"){
					forPrimer.push_back(oligo);
				}
				else if(type == "reverse"){
					Sequence oligoRC("reverse", oligo);
					oligoRC.reverseComplement();
					revPrimer.push_back(oligoRC.getUnaligned());
				}
				else if(type == "barcode"){
					inOligos >> group;
					barcodes[oligo]=index++;
					groupVector.push_back(group);
					
					outSFFFlowVec.push_back(new ofstream((outputDir + m->getRootName(m->getSimpleName(sffFile)) + group + ".flow").c_str(), ios::ate));
					outputNames.push_back((outputDir + m->getRootName(m->getSimpleName(sffFile)) + group + "flow"));
					outputTypes["flow"].push_back((outputDir + m->getRootName(m->getSimpleName(sffFile)) + group + "flow"));
				}
			}
			m->gobble(inOligos);
		}
		
		inOligos.close();
		
		numFPrimers = forPrimer.size();
		numRPrimers = revPrimer.size();
		numBarcodes = barcodes.size();
	}
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "getOligos");
		exit(1);
	}
	
}

//**********************************************************************************************************************

int ParseSFFCommand::parseHeaderLineToInt(ifstream& file){
	
	int number;

	while (!file.eof())	{

		char c = file.get(); 
		if (c == ':'){
			file >> number;
			break;
		}
		
	}
	m->gobble(file);
	return number;
}

//**********************************************************************************************************************

string ParseSFFCommand::parseHeaderLineToString(ifstream& file){
	
	string text;
	
	while (!file.eof())	{
		char c = file.get(); 
		
		if (c == ':'){
			m->gobble(file);
			text = m->getline(file);			
			break;
		}
	}
	m->gobble(file);

	return text;
}

//**********************************************************************************************************************

vector<float> ParseSFFCommand::parseHeaderLineToFloatVector(ifstream& file, int length){
	
	vector<float> floatVector(length);
	
	while (!file.eof())	{
		char c = file.get(); 
		if (c == ':'){
			for(int i=0;i<length;i++){
				file >> floatVector[i];
			}
			break;
		}
	}
	m->gobble(file);	
	return floatVector;
}

//**********************************************************************************************************************

vector<int> ParseSFFCommand::parseHeaderLineToIntVector(ifstream& file, int length){
	
	vector<int> intVector(length);
	
	while (!file.eof())	{
		char c = file.get(); 
		if (c == ':'){
			for(int i=0;i<length;i++){
				file >> intVector[i];
			}
			break;
		}
	}
	m->gobble(file);	
	return intVector;
}

//**********************************************************************************************************************


void ParseSFFCommand::screenFlow(vector<float> flowgram, int& length){
	try{

		int newLength = 0;

		while(newLength * 4 < length){
			
			int signal = 0;
			int noise = 0;
			for(int i=0;i<4;i++){
				float flow = flowgram[i + 4 * newLength];

				if(flow > 0.50){
					signal++;
					if(flow <= 0.69){ // not sure why, but if i make it <0.70 it doesn't work...
						noise++;
					}
				}
			}
			if(noise > 0 || signal == 0){
				break;
			}			
			newLength++;
		}
		length = newLength * 4;
	}
	
	catch(exception& e) {
		m->errorOut(e, "ParseSFFCommand", "screenFlow");
		exit(1);
	}
}

//**********************************************************************************************************************

string ParseSFFCommand::flow2seq(vector<float> flowgram, int length){

	string flow = "TACG";
	string sequence = "";
	for(int i=8;i<length;i++){
		int signal = int(flowgram[i] + 0.5);
		char base = flow[ i % 4 ];
		for(int j=0;j<signal;j++){
			sequence += base;
		}
	}
	return sequence;
}

//**********************************************************************************************************************

bool ParseSFFCommand::screenSeq(string& sequence, int& group){

	int length = 1;
	group = -1;
	
	if(sequence.length() < minLength){	length = 0;	}
	
	int barcode = 1;
	int barcodeLength = 0;

	for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
		if(compareDNASeq(it->first, sequence.substr(0,(it->first).length()))){
			barcode = 1;
			barcodeLength = (it->first).size();
			group = it->second;
			break;
		}
		else{
			barcode = 0;
		}
	}
	
	int fPrimer = 1;
	for(int i=0;i<numFPrimers;i++){
		if(compareDNASeq(forPrimer[i], sequence.substr(barcodeLength,forPrimer[i].length()))){
			fPrimer = 1;
			break;
		}
		else{
			fPrimer = 0;
		}
	}

	int rPrimer = 1;
	for(int i=0;i<numRPrimers;i++){
		if(compareDNASeq(revPrimer[i], sequence.substr(sequence.length()-revPrimer[i].length(),revPrimer[i].length()))){
			rPrimer = 1;
			break;
		}
		else{
			rPrimer = 0;
		}
	}

	return fPrimer * rPrimer * length * barcode;
		
}

//**********************************************************************************************************************
	   
bool ParseSFFCommand::compareDNASeq(string oligo, string seq){
   try {
	   bool success = 1;
	   int length = oligo.length();
	   
	   for(int i=0;i<length;i++){
		   
		   if(oligo[i] != seq[i]){
			   if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C')		{	success = 0;	}
			   else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))					{	success = 0;	}
			   else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))						{	success = 0;	}
			   else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))						{	success = 0;	}
			   else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))						{	success = 0;	}
			   else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))						{	success = 0;	}
			   else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))						{	success = 0;	}
			   else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))						{	success = 0;	}
			   else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
			   else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
			   else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))	{	success = 0;	}
			   else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))	{	success = 0;	}			
			   
			   if(success == 0)	{	break;	}
		   }
		   else{
			   success = 1;
		   }
	   }
	   
	   return success;
   }
   catch(exception& e) {
	   m->errorOut(e, "TrimSeqsCommand", "compareDNASeq");
	   exit(1);
   }
}
	   
//**********************************************************************************************************************

//string ParseSFFCommand::stripSeqQual(string qScores, int start, int end){
//	
//	
//	return qScores.substr(start-1, end-start+1);
//
//}

//**********************************************************************************************************************

//string ParseSFFCommand::stripQualQual(string qScores, int start, int end){
//	
//	start--;
//	
//	int startCount = 0;
//	int startIndex = 0;
//	
//	while(startCount < start && startIndex < qScores.length()){
//		if(isspace(qScores[startIndex])){
//			startCount++;
//		}
//	   startIndex++;
//	}
//	
//	int endCount = startCount;
//	int endIndex = startIndex;
//	
//	while(endCount < end && endIndex < qScores.length()){
//		if(isspace(qScores[endIndex])){
//			endCount++;
//		}
//		endIndex++;
//	}
//	
//   return qScores.substr(startIndex, endIndex-startIndex-1);//, endCount-startCount);
//	
//}

//**********************************************************************************************************************


