/*
 *  seqerrorcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "seqerrorcommand.h"
#include "reportfile.h"
#include "qualityscores.h"
#include "refchimeratest.h"

//**********************************************************************************************************************
vector<string> SeqErrorCommand::setParameters(){	
	try {
		CommandParameter pquery("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pquery);
		CommandParameter preference("reference", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(preference);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "QualReport",false,false); parameters.push_back(pqfile);
		CommandParameter preport("report", "InputTypes", "", "", "none", "none", "QualReport",false,false); parameters.push_back(preport);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pignorechimeras("ignorechimeras", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pignorechimeras);
		CommandParameter pthreshold("threshold", "Number", "", "1.0", "", "", "",false,false); parameters.push_back(pthreshold);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SeqErrorCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The seq.error command reads a query alignment file and a reference alignment file and creates .....\n";
		helpString += "Example seq.error(...).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/seq.error .\n\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
SeqErrorCommand::SeqErrorCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["error.summary"] = tempOutNames;
		outputTypes["error.seq"] = tempOutNames;
		outputTypes["error.quality"] = tempOutNames;
		outputTypes["error.qual.forward"] = tempOutNames;
		outputTypes["error.qual.reverse"] = tempOutNames;
		outputTypes["error.forward"] = tempOutNames;
		outputTypes["error.reverse"] = tempOutNames;
		outputTypes["error.count"] = tempOutNames;
		outputTypes["error.matrix"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "SeqErrorCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SeqErrorCommand::SeqErrorCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
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
			outputTypes["error.summary"] = tempOutNames;
			outputTypes["error.seq"] = tempOutNames;
			outputTypes["error.quality"] = tempOutNames;
			outputTypes["error.qual.forward"] = tempOutNames;
			outputTypes["error.qual.reverse"] = tempOutNames;
			outputTypes["error.forward"] = tempOutNames;
			outputTypes["error.reverse"] = tempOutNames;
			outputTypes["error.count"] = tempOutNames;
			outputTypes["error.matrix"] = tempOutNames;

			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a names file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}

				it = parameters.find("qfile");
				//user has given a quality score file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("report");
				//user has given a alignment report file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["report"] = inputDir + it->second;		}
				}
				
			}
			//check for required parameters
			queryFileName = validParameter.validFile(parameters, "fasta", true);
			if (queryFileName == "not found") { 
				queryFileName = m->getFastaFile(); 
				if (queryFileName != "") { m->mothurOut("Using " + queryFileName + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fasta file and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (queryFileName == "not open") { abort = true; }	
			
			referenceFileName = validParameter.validFile(parameters, "reference", true);
			if (referenceFileName == "not found") { m->mothurOut("reference is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (referenceFileName == "not open") { abort = true; }	
			

			//check for optional parameters
			namesFileName = validParameter.validFile(parameters, "name", true);
			if(namesFileName == "not found"){	namesFileName = "";	}
			else if (namesFileName == "not open") { namesFileName = ""; abort = true; }	
			
			qualFileName = validParameter.validFile(parameters, "qfile", true);
			if(qualFileName == "not found"){	qualFileName = "";	}
			else if (qualFileName == "not open") { qualFileName = ""; abort = true; }	

			reportFileName = validParameter.validFile(parameters, "report", true);
			if(reportFileName == "not found"){	reportFileName = "";	}
			else if (reportFileName == "not open") { reportFileName = ""; abort = true; }	
			
			if((reportFileName != "" && qualFileName == "") || (reportFileName == "" && qualFileName != "")){
				m->mothurOut("if you use either a qual file or a report file, you have to have both.");
				m->mothurOutEndLine();
				abort = true; 
			}
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);
			if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(queryFileName); //if user entered a file with a path then preserve it	
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found") { temp = "1.00"; }
			convert(temp, threshold);  
			
			temp = validParameter.validFile(parameters, "ignorechimeras", false);	if (temp == "not found") { temp = "1"; }
			convert(temp, ignoreChimeras);  

			substitutionMatrix.resize(6);
			for(int i=0;i<6;i++){	substitutionMatrix[i].resize(6,0);	}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "SeqErrorCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SeqErrorCommand::execute(){
	try{
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		maxLength = 2000;
		
		string errorSummaryFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.summary";
		m->openOutputFile(errorSummaryFileName, errorSummaryFile);
		outputNames.push_back(errorSummaryFileName); outputTypes["error.summary"].push_back(errorSummaryFileName);
		printErrorHeader();
		
		string errorSeqFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.seq";
		m->openOutputFile(errorSeqFileName, errorSeqFile);
		outputNames.push_back(errorSeqFileName); outputTypes["error.seq"].push_back(errorSeqFileName);

		getReferences();	//read in reference sequences - make sure there's no ambiguous bases

		map<string, int> weights;
		if(namesFileName != ""){	weights = getWeights();	}
		
		ifstream queryFile;
		m->openInputFile(queryFileName, queryFile);
		
		ifstream reportFile;
		ifstream qualFile;

		ReportFile report;
		QualityScores quality;
		vector<vector<int> > qualForwardMap;
		vector<vector<int> > qualReverseMap;
		
		if(qualFileName != "" && reportFileName != ""){
			m->openInputFile(qualFileName, qualFile);
			report = ReportFile(reportFile, reportFileName);
			
			qualForwardMap.resize(maxLength);
			qualReverseMap.resize(maxLength);
			for(int i=0;i<maxLength;i++){
				qualForwardMap[i].assign(41,0);
				qualReverseMap[i].assign(41,0);
			}				
		}
		
		int totalBases = 0;
		int totalMatches = 0;
		
		vector<int> misMatchCounts(11, 0);
		int maxMismatch = 0;
		int numSeqs = 0;
		
		map<string, int>::iterator it;
		map<char, vector<int> > qScoreErrorMap;
		qScoreErrorMap['m'].assign(41, 0);
		qScoreErrorMap['s'].assign(41, 0);
		qScoreErrorMap['i'].assign(41, 0);
		qScoreErrorMap['a'].assign(41, 0);
		
		map<char, vector<int> > errorForward;
		errorForward['m'].assign(maxLength,0);
		errorForward['s'].assign(maxLength,0);
		errorForward['i'].assign(maxLength,0);
		errorForward['d'].assign(maxLength,0);
		errorForward['a'].assign(maxLength,0);
		
		map<char, vector<int> > errorReverse;
		errorReverse['m'].assign(maxLength,0);
		errorReverse['s'].assign(maxLength,0);
		errorReverse['i'].assign(maxLength,0);
		errorReverse['d'].assign(maxLength,0);
		errorReverse['a'].assign(maxLength,0);	
		

		string errorChimeraFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.chimera";
		RefChimeraTest chimeraTest(referenceSeqs, errorChimeraFileName);
		outputNames.push_back(errorChimeraFileName); outputTypes["error.chimera"].push_back(errorChimeraFileName);
		
		vector<string> megaAlignVector(numRefs, "");

		int index = 0;
		bool ignoreSeq = 0;
		
		while(queryFile){

			if (m->control_pressed) { errorSummaryFile.close();	errorSeqFile.close(); for (int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); } return 0; }
		
			Sequence query(queryFile);
			
			int numParentSeqs = chimeraTest.analyzeQuery(query.getName(), query.getAligned());
			int closestRefIndex = chimeraTest.getClosestRefIndex();

			if(numParentSeqs > 1 && ignoreChimeras == 1)	{	ignoreSeq = 1;	}
			else											{	ignoreSeq = 0;	}

			Compare minCompare = getErrors(query, referenceSeqs[closestRefIndex]);
			
			if(namesFileName != ""){
				it = weights.find(query.getName());
				minCompare.weight = it->second;
			}
			else{	minCompare.weight = 1;	}

			printErrorData(minCompare, numParentSeqs);

			if(!ignoreSeq){
				
				for(int i=0;i<minCompare.total;i++){
					char letter = minCompare.sequence[i];

					errorForward[letter][i] += minCompare.weight;
					errorReverse[letter][minCompare.total-i-1] += minCompare.weight;				
				}
			}

			if(qualFileName != "" && reportFileName != ""){
				report = ReportFile(reportFile);
				
//				int origLength = report.getQueryLength();
				int startBase = report.getQueryStart();
				int endBase = report.getQueryEnd();

				quality = QualityScores(qualFile);

				if(!ignoreSeq){
					quality.updateQScoreErrorMap(qScoreErrorMap, minCompare.sequence, startBase, endBase, minCompare.weight);
					quality.updateForwardMap(qualForwardMap, startBase, endBase, minCompare.weight);
					quality.updateReverseMap(qualReverseMap, startBase, endBase, minCompare.weight);
				}
			}			

			if(minCompare.errorRate < threshold && !ignoreSeq){
				totalBases += (minCompare.total * minCompare.weight);
				totalMatches += minCompare.matches * minCompare.weight;
				if(minCompare.mismatches > maxMismatch){
					maxMismatch = minCompare.mismatches;
					misMatchCounts.resize(maxMismatch + 1, 0);
				}				
				misMatchCounts[minCompare.mismatches] += minCompare.weight;
				numSeqs++;
				
				megaAlignVector[closestRefIndex] += query.getInlineSeq() + '\n';
			}

			index++;
			
			if(index % 1000 == 0){	m->mothurOut(toString(index) + '\n');	}
		}
		queryFile.close();
		errorSummaryFile.close();	
		errorSeqFile.close();

		if(qualFileName != "" && reportFileName != ""){		
			printErrorQuality(qScoreErrorMap);
			printQualityFR(qualForwardMap, qualReverseMap);
		}
		
		printErrorFRFile(errorForward, errorReverse);
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); } return 0; }

		string errorCountFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.count";
		ofstream errorCountFile;
		m->openOutputFile(errorCountFileName, errorCountFile);
		outputNames.push_back(errorCountFileName);  outputTypes["error.count"].push_back(errorCountFileName);
		m->mothurOut("Overall error rate:\t" + toString((double)(totalBases - totalMatches) / (double)totalBases) + "\n\n");
		m->mothurOut("Errors\tSequences\n");
		errorCountFile << "Errors\tSequences\n";		
		for(int i=0;i<misMatchCounts.size();i++){
			m->mothurOut(toString(i) + '\t' + toString(misMatchCounts[i]) + '\n');
			errorCountFile << i << '\t' << misMatchCounts[i] << endl;
		}
		errorCountFile.close();
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); } return 0; }

		printSubMatrix();
				
		string megAlignmentFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.ref-query";
		ofstream megAlignmentFile;
		m->openOutputFile(megAlignmentFileName, megAlignmentFile);
		outputNames.push_back(megAlignmentFileName);  outputTypes["error.ref-query"].push_back(megAlignmentFileName);
		
		for(int i=0;i<numRefs;i++){
			megAlignmentFile << referenceSeqs[i].getInlineSeq() << endl;
			megAlignmentFile << megaAlignVector[i] << endl;
		}
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
		
		return 0;	
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::getReferences(){
	try {
		
		ifstream referenceFile;
		m->openInputFile(referenceFileName, referenceFile);
		
		int numAmbigSeqs = 0;
		
		int maxStartPos = 0;
		int minEndPos = 100000;
		
		while(referenceFile){
			Sequence currentSeq(referenceFile);
			int numAmbigs = currentSeq.getAmbigBases();
			if(numAmbigs > 0){	numAmbigSeqs++;	}
			
			int startPos = currentSeq.getStartPos();
			if(startPos > maxStartPos)	{	maxStartPos = startPos;	}

			int endPos = currentSeq.getEndPos();
			if(endPos < minEndPos)		{	minEndPos = endPos;		}
			referenceSeqs.push_back(currentSeq);
			m->gobble(referenceFile);
		}
		referenceFile.close();
		numRefs = referenceSeqs.size();

		
		for(int i=0;i<numRefs;i++){
			referenceSeqs[i].padToPos(maxStartPos);
			referenceSeqs[i].padFromPos(minEndPos);
		}
		
		if(numAmbigSeqs != 0){
			m->mothurOut("Warning: " + toString(numAmbigSeqs) + " reference sequences have ambiguous bases, these bases will be ignored\n");
		}		
		
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getReferences");
		exit(1);
	}
}

//***************************************************************************************************************

Compare SeqErrorCommand::getErrors(Sequence query, Sequence reference){
	try {
		if(query.getAlignLength() != reference.getAlignLength()){
			m->mothurOut("Warning: " + toString(query.getName()) + " and " + toString(reference.getName()) + " are different lengths\n");
		}
		int alignLength = query.getAlignLength();
	
		string q = query.getAligned();
		string r = reference.getAligned();

		int started = 0;
		Compare errors;

		for(int i=0;i<alignLength;i++){
			if(r[i] != 'N' && q[i] != '.' && r[i] != '.' && (q[i] != '-' || r[i] != '-')){			//	no missing data and no double gaps
				started = 1;
				
				if(q[i] == 'A'){
					if(r[i] == 'A'){	errors.AA++;	errors.matches++;	errors.sequence += 'm';	}
					if(r[i] == 'T'){	errors.AT++;	errors.sequence += 's';	}
					if(r[i] == 'G'){	errors.AG++;	errors.sequence += 's';	}
					if(r[i] == 'C'){	errors.AC++;	errors.sequence += 's';	}
					if(r[i] == '-'){	errors.Ai++;	errors.sequence += 'i';	}
				}
				else if(q[i] == 'T'){
					if(r[i] == 'A'){	errors.TA++;	errors.sequence += 's';	}
					if(r[i] == 'T'){	errors.TT++;	errors.matches++;	errors.sequence += 'm';	}
					if(r[i] == 'G'){	errors.TG++;	errors.sequence += 's';	}
					if(r[i] == 'C'){	errors.TC++;	errors.sequence += 's';	}
					if(r[i] == '-'){	errors.Ti++;	errors.sequence += 'i';	}
				}
				else if(q[i] == 'G'){
					if(r[i] == 'A'){	errors.GA++;	errors.sequence += 's';	}
					if(r[i] == 'T'){	errors.GT++;	errors.sequence += 's';	}
					if(r[i] == 'G'){	errors.GG++;	errors.matches++;	errors.sequence += 'm';	}
					if(r[i] == 'C'){	errors.GC++;	errors.sequence += 's';	}
					if(r[i] == '-'){	errors.Gi++;	errors.sequence += 'i';	}
				}
				else if(q[i] == 'C'){
					if(r[i] == 'A'){	errors.CA++;	errors.sequence += 's';	}
					if(r[i] == 'T'){	errors.CT++;	errors.sequence += 's';	}
					if(r[i] == 'G'){	errors.CG++;	errors.sequence += 's';	}
					if(r[i] == 'C'){	errors.CC++;	errors.matches++;	errors.sequence += 'm';	}
					if(r[i] == '-'){	errors.Ci++;	errors.sequence += 'i';	}
				}
				else if(q[i] == 'N'){
					if(r[i] == 'A'){	errors.NA++;	errors.sequence += 'a';	}
					if(r[i] == 'T'){	errors.NT++;	errors.sequence += 'a';	}
					if(r[i] == 'G'){	errors.NG++;	errors.sequence += 'a';	}
					if(r[i] == 'C'){	errors.NC++;	errors.sequence += 'a';	}
					if(r[i] == '-'){	errors.Ni++;	errors.sequence += 'a';	}
				}
				else if(q[i] == '-' && r[i] != '-'){
					if(r[i] == 'A'){	errors.dA++;	errors.sequence += 'd';	}
					if(r[i] == 'T'){	errors.dT++;	errors.sequence += 'd';	}
					if(r[i] == 'G'){	errors.dG++;	errors.sequence += 'd';	}
					if(r[i] == 'C'){	errors.dC++;	errors.sequence += 'd';	}
				}
				errors.total++;	
				
			}
			else if(q[i] == '.' && r[i] != '.'){		//	reference extends beyond query
				if(started == 1){	break;	}
			}
			else if(q[i] != '.' && r[i] == '.'){		//	query extends beyond reference
				if(started == 1){	break;	}
			}
			else if(q[i] == '.' && r[i] == '.'){		//	both are missing data
				if(started == 1){	break;	}			
			}
			
		}
		errors.mismatches = errors.total-errors.matches;
		errors.errorRate = (double)(errors.total-errors.matches) / (double)errors.total;
		errors.queryName = query.getName();
		errors.refName = reference.getName();
		
		return errors;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getErrors");
		exit(1);
	}
}

//***************************************************************************************************************

map<string, int> SeqErrorCommand::getWeights(){
	ifstream nameFile;
	m->openInputFile(namesFileName, nameFile);
	
	string seqName;
	string redundantSeqs;
	map<string, int> nameCountMap;
	
	while(nameFile){
		nameFile >> seqName >> redundantSeqs;
		nameCountMap[seqName] = m->getNumNames(redundantSeqs); 
		m->gobble(nameFile);
	}
	return nameCountMap;
}


//***************************************************************************************************************

void SeqErrorCommand::printErrorHeader(){
	try {
		errorSummaryFile << "query\treference\tweight\t";
		errorSummaryFile << "AA\tAT\tAG\tAC\tTA\tTT\tTG\tTC\tGA\tGT\tGG\tGC\tCA\tCT\tCG\tCC\tNA\tNT\tNG\tNC\tAi\tTi\tGi\tCi\tNi\tdA\tdT\tdG\tdC\t";
		errorSummaryFile << "insertions\tdeletions\tsubstitutions\tambig\tmatches\tmismatches\ttotal\terror\tnumparents\n";
		
		errorSummaryFile << setprecision(6);
		errorSummaryFile.setf(ios::fixed);
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorHeader");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::printErrorData(Compare error, int numParentSeqs){
	try {

		errorSummaryFile << error.queryName << '\t' << error.refName << '\t' << error.weight << '\t';
		errorSummaryFile << error.AA << '\t' << error.AT << '\t' << error.AG << '\t' << error.AC << '\t';
		errorSummaryFile << error.TA << '\t' << error.TT << '\t' << error.TG << '\t' << error.TC << '\t';
		errorSummaryFile << error.GA << '\t' << error.GT << '\t' << error.GG << '\t' << error.GC << '\t';
		errorSummaryFile << error.CA << '\t' << error.CT << '\t' << error.CG << '\t' << error.CC << '\t';
		errorSummaryFile << error.NA << '\t' << error.NT << '\t' << error.NG << '\t' << error.NC << '\t';
		errorSummaryFile << error.Ai << '\t' << error.Ti << '\t' << error.Gi << '\t' << error.Ci << '\t' << error.Ni << '\t';
		errorSummaryFile << error.dA << '\t' << error.dT << '\t' << error.dG << '\t' << error.dC << '\t';
		
		errorSummaryFile << error.Ai + error.Ti + error.Gi + error.Ci << '\t';			//insertions
		errorSummaryFile << error.dA + error.dT + error.dG + error.dC << '\t';			//deletions
		errorSummaryFile << error.mismatches - (error.Ai + error.Ti + error.Gi + error.Ci) - (error.dA + error.dT + error.dG + error.dC) - (error.NA + error.NT + error.NG + error.NC + error.Ni) << '\t';	//substitutions
		errorSummaryFile << error.NA + error.NT + error.NG + error.NC + error.Ni << '\t';	//ambiguities
		errorSummaryFile << error.matches << '\t' << error.mismatches << '\t' << error.total << '\t' << error.errorRate << '\t' << numParentSeqs << endl;

		errorSeqFile << '>' << error.queryName << "\tref:" << error.refName << '\n' << error.sequence << endl;
		
		int a=0;		int t=1;		int g=2;		int c=3;
		int gap=4;		int n=5;

		if(numParentSeqs == 1 || ignoreChimeras == 0){
			substitutionMatrix[a][a] += error.weight * error.AA;
			substitutionMatrix[a][t] += error.weight * error.TA;
			substitutionMatrix[a][g] += error.weight * error.GA;
			substitutionMatrix[a][c] += error.weight * error.CA;
			substitutionMatrix[a][gap] += error.weight * error.dA;
			substitutionMatrix[a][n] += error.weight * error.NA;
			
			substitutionMatrix[t][a] += error.weight * error.AT;
			substitutionMatrix[t][t] += error.weight * error.TT;
			substitutionMatrix[t][g] += error.weight * error.GT;
			substitutionMatrix[t][c] += error.weight * error.CT;
			substitutionMatrix[t][gap] += error.weight * error.dT;
			substitutionMatrix[t][n] += error.weight * error.NT;

			substitutionMatrix[g][a] += error.weight * error.AG;
			substitutionMatrix[g][t] += error.weight * error.TG;
			substitutionMatrix[g][g] += error.weight * error.GG;
			substitutionMatrix[g][c] += error.weight * error.CG;
			substitutionMatrix[g][gap] += error.weight * error.dG;
			substitutionMatrix[g][n] += error.weight * error.NG;

			substitutionMatrix[c][a] += error.weight * error.AC;
			substitutionMatrix[c][t] += error.weight * error.TC;
			substitutionMatrix[c][g] += error.weight * error.GC;
			substitutionMatrix[c][c] += error.weight * error.CC;
			substitutionMatrix[c][gap] += error.weight * error.dC;
			substitutionMatrix[c][n] += error.weight * error.NC;

			substitutionMatrix[gap][a] += error.weight * error.Ai;
			substitutionMatrix[gap][t] += error.weight * error.Ti;
			substitutionMatrix[gap][g] += error.weight * error.Gi;
			substitutionMatrix[gap][c] += error.weight * error.Ci;
			substitutionMatrix[gap][n] += error.weight * error.Ni;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorData");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::printSubMatrix(){
	try {
		string subMatrixFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.matrix";
		ofstream subMatrixFile;
		m->openOutputFile(subMatrixFileName, subMatrixFile);
		outputNames.push_back(subMatrixFileName);  outputTypes["error.matrix"].push_back(subMatrixFileName);
		vector<string> bases(6);
		bases[0] = "A";
		bases[1] = "T";
		bases[2] = "G";
		bases[3] = "C";
		bases[4] = "Gap";
		bases[5] = "N";
		vector<int> refSums(5,1);

		for(int i=0;i<5;i++){
			subMatrixFile << "\tr" << bases[i];
			
			for(int j=0;j<6;j++){
				refSums[i] += substitutionMatrix[i][j];				
			}
		}
		subMatrixFile << endl;
		
		for(int i=0;i<6;i++){
			subMatrixFile << 'q' << bases[i];
			for(int j=0;j<5;j++){
				subMatrixFile << '\t' << substitutionMatrix[j][i];				
			}
			subMatrixFile << endl;
		}

		subMatrixFile << "total";
		for(int i=0;i<5;i++){
			subMatrixFile << '\t' << refSums[i];
		}
		subMatrixFile << endl;
		subMatrixFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printSubMatrix");
		exit(1);
	}
}
//***************************************************************************************************************

void SeqErrorCommand::printErrorFRFile(map<char, vector<int> > errorForward, map<char, vector<int> > errorReverse){
	try{
		string errorForwardFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.seq.forward";
		ofstream errorForwardFile;
		m->openOutputFile(errorForwardFileName, errorForwardFile);
		outputNames.push_back(errorForwardFileName);  outputTypes["error.forward"].push_back(errorForwardFileName);

		errorForwardFile << "position\ttotalseqs\tmatch\tsubstitution\tinsertion\tdeletion\tambiguous" << endl;
		for(int i=0;i<maxLength;i++){
			float match = (float)errorForward['m'][i];
			float subst = (float)errorForward['s'][i];
			float insert = (float)errorForward['i'][i];
			float del = (float)errorForward['d'][i];
			float amb = (float)errorForward['a'][i];
			float total = match + subst + insert + del + amb;
			if(total == 0){	break;	}
			errorForwardFile << i+1 << '\t' << total << '\t' << match/total  << '\t' << subst/total  << '\t' << insert/total  << '\t' << del/total  << '\t' << amb/total << endl;
		}
		errorForwardFile.close();

		string errorReverseFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.seq.reverse";
		ofstream errorReverseFile;
		m->openOutputFile(errorReverseFileName, errorReverseFile);
		outputNames.push_back(errorReverseFileName);  outputTypes["error.reverse"].push_back(errorReverseFileName);

		errorReverseFile << "position\ttotalseqs\tmatch\tsubstitution\tinsertion\tdeletion\tambiguous" << endl;
		for(int i=0;i<maxLength;i++){
			float match = (float)errorReverse['m'][i];
			float subst = (float)errorReverse['s'][i];
			float insert = (float)errorReverse['i'][i];
			float del = (float)errorReverse['d'][i];
			float amb = (float)errorReverse['a'][i];
			float total = match + subst + insert + del + amb;
			if(total == 0){	break;	}
			errorReverseFile << i+1 << '\t' << total << '\t' << match/total  << '\t' << subst/total  << '\t' << insert/total  << '\t' << del/total  << '\t' << amb/total << endl;
		}
		errorReverseFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorFRFile");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::printErrorQuality(map<char, vector<int> > qScoreErrorMap){
	try{

		string errorQualityFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.quality";
		ofstream errorQualityFile;
		m->openOutputFile(errorQualityFileName, errorQualityFile);
		outputNames.push_back(errorQualityFileName);  outputTypes["error.quality"].push_back(errorQualityFileName);

		errorQualityFile << "qscore\tmatches\tsubstitutions\tinsertions\tambiguous" << endl;
		for(int i=0;i<41;i++){
			errorQualityFile << i << '\t' << qScoreErrorMap['m'][i] << '\t' << qScoreErrorMap['s'][i] << '\t' << qScoreErrorMap['i'][i] << '\t'<< qScoreErrorMap['a'][i] << endl;
		}
		errorQualityFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorFRFile");
		exit(1);
	}
}


//***************************************************************************************************************

void SeqErrorCommand::printQualityFR(vector<vector<int> > qualForwardMap, vector<vector<int> > qualReverseMap){

	try{
		int numRows = 0;
		int numColumns = qualForwardMap[0].size();

		for(int i=0;i<qualForwardMap.size();i++){
			for(int j=0;j<numColumns;j++){
				if(qualForwardMap[i][j] != 0){
					if(numRows < i)		{	numRows = i+20;		}
				}
			}
		}

		string qualityForwardFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.qual.forward";
		ofstream qualityForwardFile;
		m->openOutputFile(qualityForwardFileName, qualityForwardFile);
		outputNames.push_back(qualityForwardFileName);  outputTypes["error.qual.forward"].push_back(qualityForwardFileName);

		for(int i=0;i<numColumns;i++){	qualityForwardFile << '\t' << i;	}	qualityForwardFile << endl;

		for(int i=0;i<numRows;i++){
			qualityForwardFile << i+1;
			for(int j=0;j<numColumns;j++){
				qualityForwardFile << '\t' << qualForwardMap[i][j];
			}

			qualityForwardFile << endl;
		}
		qualityForwardFile.close();

		
		string qualityReverseFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.qual.reverse";
		ofstream qualityReverseFile;
		m->openOutputFile(qualityReverseFileName, qualityReverseFile);
		outputNames.push_back(qualityReverseFileName);  outputTypes["error.qual.reverse"].push_back(qualityReverseFileName);
		
		for(int i=0;i<numColumns;i++){	qualityReverseFile << '\t' << i;	}	qualityReverseFile << endl;
		for(int i=0;i<numRows;i++){
			
			qualityReverseFile << i+1;
			for(int j=0;j<numColumns;j++){
				qualityReverseFile << '\t' << qualReverseMap[i][j];
			}
			qualityReverseFile << endl;
		}
		qualityReverseFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorFRFile");
		exit(1);
	}
	
}

//***************************************************************************************************************
