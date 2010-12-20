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

//**********************************************************************************************************************
vector<string> SeqErrorCommand::getValidParameters(){	
	try {
		string Array[] =  {"query", "reference", "name", "qfile", "report", "threshold", "inputdir", "outputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SeqErrorCommand::SeqErrorCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["error"] = tempOutNames;
		outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "SeqErrorCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SeqErrorCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"query","reference"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SeqErrorCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getRequiredFiles");
		exit(1);
	}
}
//***************************************************************************************************************

SeqErrorCommand::SeqErrorCommand(string option)  {
	try {
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			string temp;
			
			//valid paramters for this command
			string AlignArray[] =  {"query", "reference", "name", "qfile", "report", "threshold", "inputdir", "outputdir"};
			
//need to implement name file option
			
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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
			outputTypes["error"] = tempOutNames;
			outputTypes["count"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("query");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["query"] = inputDir + it->second;		}
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
			queryFileName = validParameter.validFile(parameters, "query", true);
			if (queryFileName == "not found") { m->mothurOut("query is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (queryFileName == "not open") { abort = true; }	
			
			referenceFileName = validParameter.validFile(parameters, "reference", true);
			if (referenceFileName == "not found") { m->mothurOut("reference is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (referenceFileName == "not open") { abort = true; }	
			

			//check for optional parameters
			namesFileName = validParameter.validFile(parameters, "name", true);
			if(namesFileName == "not found"){	namesFileName = "";	}
			
			qualFileName = validParameter.validFile(parameters, "qfile", true);
			if(qualFileName == "not found"){	qualFileName = "";	}

			reportFileName = validParameter.validFile(parameters, "report", true);
			if(reportFileName == "not found"){	reportFileName = "";	}
			
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
						
			errorSummaryFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.summary";
			m->openOutputFile(errorSummaryFileName, errorSummaryFile);
			outputNames.push_back(errorSummaryFileName); outputTypes["error.summary"].push_back(errorSummaryFileName);
			printErrorHeader();

			errorSeqFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.seq";
			m->openOutputFile(errorSeqFileName, errorSeqFile);
			outputNames.push_back(errorSeqFileName); outputTypes["error.seq"].push_back(errorSeqFileName);
			printErrorHeader();
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "SeqErrorCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void SeqErrorCommand::help(){
	try {
		m->mothurOut("The seq.error command reads a query alignment file and a reference alignment file and creates .....\n");
		
		
		
		m->mothurOut("Example seq.error(...).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n");
		m->mothurOut("For more details please check out the wiki http://www.mothur.org/wiki/seq.error .\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

SeqErrorCommand::~SeqErrorCommand(){
	errorSummaryFile.close();	
	errorSeqFile.close();	
}

//***************************************************************************************************************

int SeqErrorCommand::execute(){
	try{
		if (abort == true) { return 0; }

		getReferences();	//read in reference sequences - make sure there's no ambiguous bases

		map<string, int> weights;
		if(namesFileName != ""){	weights = getWeights();	}
		
		ifstream queryFile;
		m->openInputFile(queryFileName, queryFile);
		
		ifstream reportFile;
		ifstream qualFile;

		ReportFile report;
		QualityScores quality;

		if(qualFileName != "" && reportFileName != ""){
			m->openInputFile(qualFileName, qualFile);
			report = ReportFile(reportFile, reportFileName);
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
		
		while(queryFile){
			Compare minCompare;
			Sequence query(queryFile);
			
			for(int i=0;i<numRefs;i++){
				Compare currCompare = getErrors(query, referenceSeqs[i]);
				
				if(currCompare.errorRate < minCompare.errorRate){
					minCompare = currCompare;
				}
			}

			if(namesFileName != ""){
				it = weights.find(query.getName());
				minCompare.weight = it->second;
			}
			else	{	minCompare.weight = 1;	}

			printErrorData(minCompare);

			if(qualFileName != "" && reportFileName != ""){
				report = ReportFile(reportFile);
				
				int origLength = report.getQueryLength();
				int startBase = report.getQueryStart();
				int endBase = report.getQueryEnd();

				quality = QualityScores(qualFile, origLength);
				quality.updateQScoreErrorMap(qScoreErrorMap, minCompare.sequence, startBase, endBase, minCompare.weight);
			}			
			
			if(minCompare.errorRate < threshold){
				totalBases += (minCompare.total * minCompare.weight);
				totalMatches += minCompare.matches * minCompare.weight;
				if(minCompare.mismatches > maxMismatch){
					maxMismatch = minCompare.mismatches;
					misMatchCounts.resize(maxMismatch + 1, 0);
				}				
				misMatchCounts[minCompare.mismatches] += minCompare.weight;
				numSeqs++;
			}
			
			
		}
		queryFile.close();
		
		int total = 0;
		
		if(qualFileName != "" && reportFileName != ""){
			string errorQualityFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".error.quality";
			ofstream errorQualityFile;
			m->openOutputFile(errorQualityFileName, errorQualityFile);
			outputNames.push_back(errorQualityFileName);  outputTypes["error.quality"].push_back(errorQualityFileName);
			
			errorQualityFile << "qscore\tmatches\tsubstitutions\tinsertions\tambiguous" << endl;
			for(int i=0;i<41;i++){
				errorQualityFile << i << '\t' << qScoreErrorMap['m'][i] << '\t' << qScoreErrorMap['s'][i] << '\t' << qScoreErrorMap['i'][i] << '\t'<< qScoreErrorMap['a'][i] << endl;
			}
		}
		
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
		
		while(referenceFile){
			Sequence currentSeq(referenceFile);
			int numAmbigs = currentSeq.getAmbigBases();
			
			if(numAmbigs != 0){
				m->mothurOut("Warning: " + toString(currentSeq.getName()) + " has " + toString(numAmbigs) + " ambiguous bases, these bases will be removed\n");
				currentSeq.removeAmbigBases();
			}
			referenceSeqs.push_back(currentSeq);
			m->gobble(referenceFile);
		}
		numRefs = referenceSeqs.size();
		
		referenceFile.close();
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
			if(q[i] != '.' && r[i] != '.' && (q[i] != '-' || r[i] != '-')){			//	no missing data and no double gaps
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
				m->mothurOut("Warning: " + toString(query.getName()) + " extend beyond " + toString(reference.getName()) + ".  Ignoring the extra bases in the query\n");
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
		errorSummaryFile << "insertions\tdeletions\tsubstitutions\tambig\tmatches\tmismatches\ttotal\terror\n";
		
		errorSummaryFile << setprecision(6);
		errorSummaryFile.setf(ios::fixed);
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorHeader");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::printErrorData(Compare error){
	try {
		errorSummaryFile << error.queryName << '\t' << error.refName << '\t' << error.weight << '\t';
		errorSummaryFile << error.AA << '\t' << error.AT << '\t' << error.AG << '\t' << error.AC << '\t';
		errorSummaryFile << error.TA << '\t' << error.TT << '\t' << error.TG << '\t' << error.TC << '\t';
		errorSummaryFile << error.GA << '\t' << error.GT << '\t' << error.GG << '\t' << error.GC << '\t';
		errorSummaryFile << error.CA << '\t' << error.CT << '\t' << error.CG << '\t' << error.CC << '\t';
		errorSummaryFile << error.NA << '\t' << error.NT << '\t' << error.NG << '\t' << error.NC << '\t';
		errorSummaryFile << error.Ai << '\t' << error.Ti << '\t' << error.Gi << '\t' << error.Ci << '\t' << error.Ni << '\t' ;
		errorSummaryFile << error.dA << '\t' << error.dT << '\t' << error.dG << '\t' << error.dC << '\t';
		
		errorSummaryFile << error.Ai + error.Ti + error.Gi + error.Ci << '\t';			//insertions
		errorSummaryFile << error.dA + error.dT + error.dG + error.dC << '\t';			//deletions
		errorSummaryFile << error.mismatches - (error.Ai + error.Ti + error.Gi + error.Ci) - (error.dA + error.dT + error.dG + error.dC) - (error.NA + error.NT + error.NG + error.NC + error.Ni) << '\t';	//substitutions
		errorSummaryFile << error.NA + error.NT + error.NG + error.NC + error.Ni << '\t';	//ambiguities
		errorSummaryFile << error.matches << '\t' << error.mismatches << '\t' << error.total << '\t' << error.errorRate << endl;
		
		errorSeqFile << '>' << error.queryName << "\tref:" << error.refName << '\n' << error.sequence << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorData");
		exit(1);
	}
}

//***************************************************************************************************************







