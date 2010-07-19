/*
 *  seqerrorcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "seqerrorcommand.h"

//***************************************************************************************************************

SeqErrorCommand::SeqErrorCommand(string option)  {
	try {
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			string temp;
			
			//valid paramters for this command
			string AlignArray[] =  {"query", "reference", "name", "threshold"};
			
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
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("query");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["query"] = inputDir + it->second;		}
				}
				
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
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
			queryFileName = validParameter.validFile(parameters, "query", true);
			if (queryFileName == "not found") { m->mothurOut("query is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (queryFileName == "not open") { abort = true; }	
			
			referenceFileName = validParameter.validFile(parameters, "reference", true);
			if (referenceFileName == "not found") { m->mothurOut("reference is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (referenceFileName == "not open") { abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			namesFileName = validParameter.validFile(parameters, "name", true);
			if(namesFileName == "not found"){	namesFileName = "";	}
			cout << namesFileName << endl;
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);
			if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(queryFileName); //if user entered a file with a path then preserve it	
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found") { temp = "1.00"; }
			convert(temp, threshold);  
						
			errorFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".errors";
			openOutputFile(errorFileName, errorFile);
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

SeqErrorCommand::~SeqErrorCommand(){	errorFile.close();	}

//***************************************************************************************************************

int SeqErrorCommand::execute(){
	try{
		if (abort == true) { return 0; }

		getReferences();	//read in reference sequences - make sure there's no ambiguous bases

		map<string, int> weights;
		if(namesFileName != ""){	weights = getWeights();	}
		
		ifstream queryFile;
		openInputFile(queryFileName, queryFile);
				
		int totalBases = 0;
		int totalMatches = 0;
		
		vector<int> misMatchCounts(11, 0);
		int maxMismatch = 0;
		int numSeqs = 0;
		
		map<string, int>::iterator it;
		
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
			else {
				minCompare.weight = 1;
			}

			printErrorData(minCompare);

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
		
		
		string errorCountFileName = queryFileName.substr(0,queryFileName.find_last_of('.')) + ".count";
		ofstream errorCountFile;
		openOutputFile(errorCountFileName, errorCountFile);
		
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
		openInputFile(referenceFileName, referenceFile);
		
		while(referenceFile){
			Sequence currentSeq(referenceFile);
			int numAmbigs = currentSeq.getAmbigBases();
			
			if(numAmbigs != 0){
				m->mothurOut("Warning: " + toString(currentSeq.getName()) + " has " + toString(numAmbigs) + " ambiguous bases, these bases will be removed\n");
				currentSeq.removeAmbigBases();
			}
			referenceSeqs.push_back(currentSeq);
			gobble(referenceFile);
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
					if(r[i] == 'A'){	errors.AA++;	errors.matches++;	}
					if(r[i] == 'T'){	errors.AT++;	}
					if(r[i] == 'G'){	errors.AG++;	}
					if(r[i] == 'C'){	errors.AC++;	}
					if(r[i] == '-'){	errors.Ai++;	}
				}
				else if(q[i] == 'T'){
					if(r[i] == 'A'){	errors.TA++;	}
					if(r[i] == 'T'){	errors.TT++;	errors.matches++;	}
					if(r[i] == 'G'){	errors.TG++;	}
					if(r[i] == 'C'){	errors.TC++;	}
					if(r[i] == '-'){	errors.Ti++;	}
				}
				else if(q[i] == 'G'){
					if(r[i] == 'A'){	errors.GA++;	}
					if(r[i] == 'T'){	errors.GT++;	}
					if(r[i] == 'G'){	errors.GG++;	errors.matches++;	}
					if(r[i] == 'C'){	errors.GC++;	}
					if(r[i] == '-'){	errors.Gi++;	}
				}
				else if(q[i] == 'C'){
					if(r[i] == 'A'){	errors.CA++;	}
					if(r[i] == 'T'){	errors.CT++;	}
					if(r[i] == 'G'){	errors.CG++;	}
					if(r[i] == 'C'){	errors.CC++;	errors.matches++;	}
					if(r[i] == '-'){	errors.Ci++;	}
				}
				else if(q[i] == 'N'){
					if(r[i] == 'A'){	errors.NA++;	}
					if(r[i] == 'T'){	errors.NT++;	}
					if(r[i] == 'G'){	errors.NG++;	}
					if(r[i] == 'C'){	errors.NC++;	}
					if(r[i] == '-'){	errors.Ni++;	}
				}
				else if(q[i] == '-' && r[i] != '-'){
					if(r[i] == 'A'){	errors.dA++;	}
					if(r[i] == 'T'){	errors.dT++;	}
					if(r[i] == 'G'){	errors.dG++;	}
					if(r[i] == 'C'){	errors.dC++;	}
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
	openInputFile(namesFileName, nameFile);
	
	string seqName;
	string redundantSeqs;
	map<string, int> nameCountMap;
	
	while(nameFile){
		nameFile >> seqName >> redundantSeqs;
		nameCountMap[seqName] = getNumNames(redundantSeqs); 
		gobble(nameFile);
	}
	return nameCountMap;
}


//***************************************************************************************************************

void SeqErrorCommand::printErrorHeader(){
	try {
		errorFile << "query\treference\tweight\t";
		errorFile << "AA\tAT\tAG\tAC\tTA\tTT\tTG\tTC\tGA\tGT\tGG\tGC\tCA\tCT\tCG\tCC\tNA\tNT\tNG\tNC\tAi\tTi\tGi\tCi\tNi\tdA\tdT\tdG\tdC\t";
		errorFile << "insertions\tdeletions\tsubstitutions\tambig\tmatches\tmismatches\ttotal\terror\n";
		
		errorFile << setprecision(6);
		errorFile.setf(ios::fixed);
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorHeader");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::printErrorData(Compare error){
	try {
		errorFile << error.queryName << '\t' << error.refName << '\t' << error.weight << '\t';
		errorFile << error.AA << '\t' << error.AT << '\t' << error.AG << '\t' << error.AC << '\t';
		errorFile << error.TA << '\t' << error.TT << '\t' << error.TG << '\t' << error.TC << '\t';
		errorFile << error.GA << '\t' << error.GT << '\t' << error.GG << '\t' << error.GC << '\t';
		errorFile << error.CA << '\t' << error.CT << '\t' << error.CG << '\t' << error.CC << '\t';
		errorFile << error.NA << '\t' << error.NT << '\t' << error.NG << '\t' << error.NC << '\t';
		errorFile << error.Ai << '\t' << error.Ti << '\t' << error.Gi << '\t' << error.Ci << '\t' << error.Ni << '\t' ;
		errorFile << error.dA << '\t' << error.dT << '\t' << error.dG << '\t' << error.dC << '\t';
		
		errorFile << error.Ai + error.Ti + error.Gi + error.Ci << '\t';			//insertions
		errorFile << error.dA + error.dT + error.dG + error.dC << '\t';			//deletions
		errorFile << error.mismatches - (error.Ai + error.Ti + error.Gi + error.Ci) - (error.dA + error.dT + error.dG + error.dC) - (error.NA + error.NT + error.NG + error.NC + error.Ni) << '\t';	//substitutions
		errorFile << error.NA + error.NT + error.NG + error.NC + error.Ni << '\t';	//ambiguities
		errorFile << error.matches << '\t' << error.mismatches << '\t' << error.total << '\t' << error.errorRate << endl;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorData");
		exit(1);
	}
}

//***************************************************************************************************************







