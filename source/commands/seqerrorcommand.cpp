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
		CommandParameter pquery("fasta", "InputTypes", "", "", "none", "none", "none","errorType",false,true,true); parameters.push_back(pquery);
		CommandParameter preference("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(preference);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "QualReport","",false,false); parameters.push_back(pqfile);
		CommandParameter preport("report", "InputTypes", "", "", "none", "none", "QualReport","",false,false); parameters.push_back(preport);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pignorechimeras("ignorechimeras", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pignorechimeras);
		CommandParameter pthreshold("threshold", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pthreshold);
		CommandParameter paligned("aligned", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(paligned);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The fasta parameter...\n";
		helpString += "The reference parameter...\n";
		helpString += "The qfile parameter...\n";
		helpString += "The report parameter...\n";
		helpString += "The name parameter allows you to provide a name file associated with the fasta file.\n";
        helpString += "The count parameter allows you to provide a count file associated with the fasta file.\n";
		helpString += "The ignorechimeras parameter...\n";
		helpString += "The threshold parameter...\n";
		helpString += "The processors parameter...\n";
		helpString += "Example seq.error(...).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/seq.error .\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************

string SeqErrorCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "errorsummary")            {   pattern = "[filename],error.summary";   }
        else if (type == "errorseq")            {   pattern = "[filename],error.seq";   }
        else if (type == "errorquality")            {   pattern = "[filename],error.quality";   }
        else if (type == "errorqualforward")            {   pattern = "[filename],error.qual.forward";   }
        else if (type == "errorqualreverse")            {   pattern = "[filename],error.qual.reverse";   }
        else if (type == "errorforward")            {   pattern = "[filename],error.seq.forward";   }
        else if (type == "errorreverse")            {   pattern = "[filename],error.seq.reverse";   }
        else if (type == "errorcount")            {   pattern = "[filename],error.count";   }
        else if (type == "errormatrix")            {   pattern = "[filename],error.matrix";   }
        else if (type == "errorchimera")            {   pattern = "[filename],error.chimera";   }
        else if (type == "errorref-query")            {   pattern = "[filename],error.ref-query";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SeqErrorCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************

SeqErrorCommand::SeqErrorCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["errorsummary"] = tempOutNames;
		outputTypes["errorseq"] = tempOutNames;
		outputTypes["errorquality"] = tempOutNames;
		outputTypes["errorqualforward"] = tempOutNames;
		outputTypes["errorqualreverse"] = tempOutNames;
		outputTypes["errorforward"] = tempOutNames;
		outputTypes["errorreverse"] = tempOutNames;
		outputTypes["errorcount"] = tempOutNames;
		outputTypes["errormatrix"] = tempOutNames;
        outputTypes["errorchimera"] = tempOutNames;
        outputTypes["errorref-query"] = tempOutNames;
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
			outputTypes["errorsummary"] = tempOutNames;
			outputTypes["errorseq"] = tempOutNames;
			outputTypes["errorquality"] = tempOutNames;
			outputTypes["errorqualforward"] = tempOutNames;
			outputTypes["errorqualreverse"] = tempOutNames;
			outputTypes["errorforward"] = tempOutNames;
			outputTypes["errorreverse"] = tempOutNames;
			outputTypes["errorcount"] = tempOutNames;
			outputTypes["errormatrix"] = tempOutNames;
            outputTypes["errorchimera"] = tempOutNames;
            outputTypes["errorref-query"] = tempOutNames;

			
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
                
                it = parameters.find("count");
				//user has given a names file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
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
			else if (queryFileName == "not open") { queryFileName = ""; abort = true; }	
			else { m->setFastaFile(queryFileName); }
			
			referenceFileName = validParameter.validFile(parameters, "reference", true);
			if (referenceFileName == "not found") { m->mothurOut("reference is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (referenceFileName == "not open") { abort = true; }	
			
			//check for optional parameters
			namesFileName = validParameter.validFile(parameters, "name", true);
			if(namesFileName == "not found"){	namesFileName = "";	}
			else if (namesFileName == "not open") { namesFileName = ""; abort = true; }	
			else { m->setNameFile(namesFileName); }
            
            //check for optional parameters
			countfile = validParameter.validFile(parameters, "count", true);
			if(countfile == "not found"){	countfile = "";	}
			else if (countfile == "not open") { countfile = ""; abort = true; }
			else { m->setCountTableFile(countfile); }
			
			qualFileName = validParameter.validFile(parameters, "qfile", true);
			if(qualFileName == "not found"){	qualFileName = "";	}
			else if (qualFileName == "not open") { qualFileName = ""; abort = true; }	
			else { m->setQualFile(qualFileName); }
			
			reportFileName = validParameter.validFile(parameters, "report", true);
			if(reportFileName == "not found"){	reportFileName = "";	}
			else if (reportFileName == "not open") { reportFileName = ""; abort = true; }	
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);
			if (outputDir == "not found"){ //if user entered a file with a path then preserve it
				outputDir = m->hasPath(queryFileName); }
			
            if ((countfile != "") && (namesFileName != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found") { temp = "1.00"; }
			m->mothurConvert(temp, threshold);
			
			temp = validParameter.validFile(parameters, "ignorechimeras", false);	if (temp == "not found") { temp = "T"; }
			ignoreChimeras = m->isTrue(temp);
			
            temp = validParameter.validFile(parameters, "aligned", false);			if (temp == "not found"){	temp = "t";				}
			aligned = m->isTrue(temp); 

			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 

			if ((namesFileName == "") && (queryFileName != "")){
				vector<string> files; files.push_back(queryFileName); 
				parser.getNameFile(files);
			}

            if(aligned ){
                if((reportFileName != "" && qualFileName == "") || (reportFileName == "" && qualFileName != "")){
                    m->mothurOut("if you use either a qual file or a report file, you have to have both.");
                    m->mothurOutEndLine();
                    abort = true; 
                }
			}
            else{
                if(reportFileName != ""){
                    m->mothurOut("we are ignoring the report file if your sequences are not aligned.  we will check that the sequences in your fasta and and qual file are the same length.");
                    m->mothurOutEndLine();
                }
            }
            

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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		int start = time(NULL);
		maxLength = 5000;
		totalBases = 0;
		totalMatches = 0;
        substitutionMatrix.resize(6);
        for(int i=0;i<6;i++){	substitutionMatrix[i].resize(6,0);	}
		
        string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string errorSummaryFileName = getOutputFileName("errorsummary",variables);
		outputNames.push_back(errorSummaryFileName); outputTypes["errorsummary"].push_back(errorSummaryFileName);
			
		string errorSeqFileName = getOutputFileName("errorseq",variables);
		outputNames.push_back(errorSeqFileName); outputTypes["errorseq"].push_back(errorSeqFileName);
		
		string errorChimeraFileName = getOutputFileName("errorchimera",variables);
		outputNames.push_back(errorChimeraFileName); outputTypes["errorchimera"].push_back(errorChimeraFileName);
		
		setLines(queryFileName, qualFileName, reportFileName);
		
		if (m->getControl_pressed()) { return 0; }
        
        getReferences();	//read in reference sequences - make sure there's no ambiguous bases
        
        if(namesFileName != "")     {	weights = getWeights();         }
        else if (countfile != "")   {
            CountTable ct;
            ct.readTable(countfile, false, false);
            weights = ct.getNameMap();
        }
        
		if (m->getControl_pressed()) { return 0; }
        
        int numSeqs = 0;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		if(processors == 1){
			numSeqs = driver(queryFileName, qualFileName, reportFileName, errorSummaryFileName, errorSeqFileName, errorChimeraFileName, lines[0], qLines[0], rLines[0]);
            
		}else{
			numSeqs = createProcesses(queryFileName, qualFileName, reportFileName, errorSummaryFileName, errorSeqFileName, errorChimeraFileName);
		}	
#else
		numSeqs = driver(queryFileName, qualFileName, reportFileName, errorSummaryFileName, errorSeqFileName, errorChimeraFileName, lines[0], qLines[0], rLines[0]);
#endif

		if(qualFileName != ""){		
			printErrorQuality(qScoreErrorMap);
			printQualityFR(qualForwardMap, qualReverseMap);
		}
		
		printErrorFRFile(errorForward, errorReverse);
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }

		string errorCountFileName = getOutputFileName("errorcount",variables);
		ofstream errorCountFile;
		m->openOutputFile(errorCountFileName, errorCountFile);
		outputNames.push_back(errorCountFileName);  outputTypes["errorcount"].push_back(errorCountFileName);
        m->mothurOut("\nMultiply error rate by 100 to obtain the percent sequencing errors.\n");
		m->mothurOut("Overall error rate:\t" + toString((double)(totalBases - totalMatches) / (double)totalBases) + "\n");
		m->mothurOut("Errors\tSequences\n");
		errorCountFile << "Errors\tSequences\n";		
		for(int i=0;i<misMatchCounts.size();i++){
			m->mothurOut(toString(i) + '\t' + toString(misMatchCounts[i]) + '\n');
			errorCountFile << i << '\t' << misMatchCounts[i] << endl;
		}
		errorCountFile.close();
		
//		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }

		printSubMatrix();
        
		string megAlignmentFileName = getOutputFileName("errorref-query",variables);
		ofstream megAlignmentFile;
		m->openOutputFile(megAlignmentFileName, megAlignmentFile);
		outputNames.push_back(megAlignmentFileName);  outputTypes["errorref-query"].push_back(megAlignmentFileName);
        
		for(int i=0;i<numRefs;i++){
			megAlignmentFile << referenceSeqs[i].getInlineSeq() << endl;
			megAlignmentFile << megaAlignVector[i] << endl;
		}
        megAlignmentFile.close();
        
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");
		m->mothurOutEndLine();
		
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

//**********************************************************************************************************************

int SeqErrorCommand::createProcesses(string filename, string qFileName, string rFileName, string summaryFileName, string errorOutputFileName, string chimeraOutputFileName) {	
	try {
		int process = 1;
		processIDS.clear();
		map<char, vector<int> >::iterator it;
		int num = 0;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				
				num = driver(filename, qFileName, rFileName, summaryFileName + m->mothurGetpid(process) + ".temp", errorOutputFileName+ m->mothurGetpid(process) + ".temp", chimeraOutputFileName + m->mothurGetpid(process) + ".temp", lines[process], qLines[process], rLines[process]);
				
				//pass groupCounts to parent
				ofstream out;
				string tempFile = filename + m->mothurGetpid(process) + ".info.temp";
				m->openOutputFile(tempFile, out);
				
				//output totalBases and totalMatches
				out << num << '\t' << totalBases << '\t' << totalMatches << endl << endl;
				
				//output substitutionMatrix
				for(int i = 0; i < substitutionMatrix.size(); i++) {
					for (int j = 0; j < substitutionMatrix[i].size(); j++) {
						out << substitutionMatrix[i][j] << '\t';
					}
					out << endl;
				}
				out << endl;
				
				//output qScoreErrorMap
				for (it = qScoreErrorMap.begin(); it != qScoreErrorMap.end(); it++) {
					vector<int> thisScoreErrorMap = it->second;
					out << it->first << '\t';
					for (int i = 0; i < thisScoreErrorMap.size(); i++) {
						out << thisScoreErrorMap[i] << '\t';
					}
					out << endl;
				}
				out << endl;
				
				//output qualForwardMap
				for(int i = 0; i < qualForwardMap.size(); i++) {
					for (int j = 0; j < qualForwardMap[i].size(); j++) {
						out << qualForwardMap[i][j] << '\t';
					}
					out << endl;
				}
				out << endl;
				
				//output qualReverseMap
				for(int i = 0; i < qualReverseMap.size(); i++) {
					for (int j = 0; j < qualReverseMap[i].size(); j++) {
						out << qualReverseMap[i][j] << '\t';
					}
					out << endl;
				}
				out << endl;
				
				
				//output errorForward
				for (it = errorForward.begin(); it != errorForward.end(); it++) {
					vector<int> thisErrorForward = it->second;
					out << it->first << '\t';
					for (int i = 0; i < thisErrorForward.size(); i++) {
						out << thisErrorForward[i] << '\t';
					}
					out << endl;
				}
				out << endl;
				
				//output errorReverse
				for (it = errorReverse.begin(); it != errorReverse.end(); it++) {
					vector<int> thisErrorReverse = it->second;
					out << it->first << '\t';
					for (int i = 0; i < thisErrorReverse.size(); i++) {
						out << thisErrorReverse[i] << '\t';
					}
					out << endl;
				}
				out << endl;
				
				//output misMatchCounts
				out << misMatchCounts.size() << endl;
				for (int j = 0; j < misMatchCounts.size(); j++) {
					out << misMatchCounts[j] << '\t';
				}
				out << endl;
				
				
				//output megaAlignVector
				for (int j = 0; j < megaAlignVector.size(); j++) {
					out << megaAlignVector[j] << endl;
				}
				out << endl;
				
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		num = driver(filename, qFileName, rFileName, summaryFileName, errorOutputFileName, chimeraOutputFileName, lines[0], qLines[0], rLines[0]);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//append files
		for(int i=0;i<processIDS.size();i++){
			
			m->mothurOut("Appending files from process " + toString(processIDS[i])); m->mothurOutEndLine();
			
			m->appendFiles((summaryFileName + toString(processIDS[i]) + ".temp"), summaryFileName);
			m->mothurRemove((summaryFileName + toString(processIDS[i]) + ".temp"));
			m->appendFiles((errorOutputFileName + toString(processIDS[i]) + ".temp"), errorOutputFileName);
			m->mothurRemove((errorOutputFileName + toString(processIDS[i]) + ".temp"));
			m->appendFiles((chimeraOutputFileName + toString(processIDS[i]) + ".temp"), chimeraOutputFileName);
			m->mothurRemove((chimeraOutputFileName + toString(processIDS[i]) + ".temp"));
			
			ifstream in;
			string tempFile =  filename + toString(processIDS[i]) + ".info.temp";
			m->openInputFile(tempFile, in);
			
			//input totalBases and totalMatches
			int tempBases, tempMatches, tempNumSeqs;
			in >> tempNumSeqs >> tempBases >> tempMatches; m->gobble(in);
			totalBases += tempBases; totalMatches += tempMatches; num += tempNumSeqs;
			
			//input substitutionMatrix
			int tempNum;
			for(int i = 0; i < substitutionMatrix.size(); i++) {
				for (int j = 0; j < substitutionMatrix[i].size(); j++) {
					in >> tempNum; substitutionMatrix[i][j] += tempNum;
				}
				m->gobble(in);
			}
			m->gobble(in);
			
			//input qScoreErrorMap
			char first;
			for (int i = 0; i < qScoreErrorMap.size(); i++) {
				in >> first;
				vector<int> thisScoreErrorMap = qScoreErrorMap[first];
				
				for (int i = 0; i < thisScoreErrorMap.size(); i++) {
					in >> tempNum; thisScoreErrorMap[i] += tempNum;
				}
				qScoreErrorMap[first] = thisScoreErrorMap;
				m->gobble(in);
			}
			m->gobble(in);
			
			//input qualForwardMap
			for(int i = 0; i < qualForwardMap.size(); i++) {
				for (int j = 0; j < qualForwardMap[i].size(); j++) {
					in >> tempNum; qualForwardMap[i][j] += tempNum;
				}
				m->gobble(in);
			}
			m->gobble(in);
			
			//input qualReverseMap
			for(int i = 0; i < qualReverseMap.size(); i++) {
				for (int j = 0; j < qualReverseMap[i].size(); j++) {
					in >> tempNum; qualReverseMap[i][j] += tempNum;
				}
				m->gobble(in);
			}
			m->gobble(in);
			
			//input errorForward
			for (int i = 0; i < errorForward.size(); i++) {
				in >> first;
				vector<int> thisErrorForward = errorForward[first];
				
				for (int i = 0; i < thisErrorForward.size(); i++) {
					in >> tempNum; thisErrorForward[i] += tempNum;
				}
				errorForward[first] = thisErrorForward;
				m->gobble(in);
			}
			m->gobble(in);
			
			//input errorReverse
			for (int i = 0; i < errorReverse.size(); i++) {
				in >> first;
				vector<int> thisErrorReverse = errorReverse[first];
				
				for (int i = 0; i < thisErrorReverse.size(); i++) {
					in >> tempNum; thisErrorReverse[i] += tempNum;
				}
				errorReverse[first] = thisErrorReverse;
				m->gobble(in);
			}
			m->gobble(in);
			
			//input misMatchCounts
			int misMatchSize;
			in >> misMatchSize; m->gobble(in);
			if (misMatchSize > misMatchCounts.size()) {	misMatchCounts.resize(misMatchSize, 0);	}
			for (int j = 0; j < misMatchSize; j++) {
				in >> tempNum; misMatchCounts[j] += tempNum;
			}
			m->gobble(in);
			
			//input megaAlignVector
			string thisLine;
			for (int j = 0; j < megaAlignVector.size(); j++) {
				thisLine = m->getline(in); m->gobble(in); megaAlignVector[j] += thisLine + '\n';
			}
			m->gobble(in);
			
			in.close(); m->mothurRemove(tempFile);
			
		}
#endif		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "createProcesses");
		exit(1);
	}
}

//**********************************************************************************************************************

int SeqErrorCommand::driver(string filename, string qFileName, string rFileName, string summaryFileName, string errorOutputFileName, string chimeraOutputFileName, linePair line, linePair qline, linePair rline) {	
	
	try {
        ReportFile report;
		QualityScores quality;
		
		misMatchCounts.resize(11, 0);
		int maxMismatch = 0;
		int numSeqs = 0;
		
		map<string, int>::iterator it;
		qScoreErrorMap['m'].assign(101, 0);
		qScoreErrorMap['s'].assign(101, 0);
		qScoreErrorMap['i'].assign(101, 0);
		qScoreErrorMap['a'].assign(101, 0);
		
		errorForward['m'].assign(maxLength,0);
		errorForward['s'].assign(maxLength,0);
		errorForward['i'].assign(maxLength,0);
		errorForward['d'].assign(maxLength,0);
		errorForward['a'].assign(maxLength,0);
		
		errorReverse['m'].assign(maxLength,0);
		errorReverse['s'].assign(maxLength,0);
		errorReverse['i'].assign(maxLength,0);
		errorReverse['d'].assign(maxLength,0);
		errorReverse['a'].assign(maxLength,0);	
		
		//open inputfiles and go to beginning place for this processor
		ifstream queryFile;
		m->openInputFile(filename, queryFile);
        
		queryFile.seekg(line.start);
		
		ifstream reportFile;
		ifstream qualFile;
		if((qFileName != "" && rFileName != "" && aligned)){
			m->openInputFile(qFileName, qualFile);
			qualFile.seekg(qline.start);  
			
			//gobble headers
			if (rline.start == 0) {  report.readHeaders(reportFile, rFileName); }
			else{
				m->openInputFile(rFileName, reportFile);
				reportFile.seekg(rline.start); 
			}
			
			qualForwardMap.resize(maxLength);
			qualReverseMap.resize(maxLength);
			for(int i=0;i<maxLength;i++){
				qualForwardMap[i].assign(101,0);
				qualReverseMap[i].assign(101,0);
			}	
		}
		else if(qFileName != "" && !aligned){

            m->openInputFile(qFileName, qualFile);
			qualFile.seekg(qline.start);  
			
			qualForwardMap.resize(maxLength);
			qualReverseMap.resize(maxLength);
			for(int i=0;i<maxLength;i++){
				qualForwardMap[i].assign(101,0);
				qualReverseMap[i].assign(101,0);
			}	
        }
        
		ofstream outChimeraReport;
		m->openOutputFile(chimeraOutputFileName, outChimeraReport);
		
        RefChimeraTest chimeraTest = RefChimeraTest(referenceSeqs, aligned);
        if (line.start == 0) { chimeraTest.printHeader(outChimeraReport); }        
        
		ofstream errorSummaryFile;
		m->openOutputFile(summaryFileName, errorSummaryFile);
		if (line.start == 0) { printErrorHeader(errorSummaryFile); }
		
		ofstream errorSeqFile;
		m->openOutputFile(errorOutputFileName, errorSeqFile);
		
		megaAlignVector.assign(numRefs, "");
		
		int index = 0;
		bool ignoreSeq = 0;
		
		bool moreSeqs = 1;
		while (moreSeqs) {
						
			Sequence query(queryFile);
            int numParentSeqs = -1;
            int closestRefIndex = -1;
                        
            string querySeq = query.getAligned();
            if (!aligned) {  querySeq = query.getUnaligned();  }
            
            numParentSeqs = chimeraTest.analyzeQuery(query.getName(), querySeq, outChimeraReport);
            
            closestRefIndex = chimeraTest.getClosestRefIndex();
            
            Sequence reference = referenceSeqs[closestRefIndex];
            
            reference.setAligned(chimeraTest.getClosestRefAlignment());
            query.setAligned(chimeraTest.getQueryAlignment());
            
			if(numParentSeqs > 1 && ignoreChimeras == 1)	{	ignoreSeq = 1;	}
            else											{	ignoreSeq = 0;	}

			Compare minCompare = getErrors(query, reference);
			
			if((namesFileName != "") || (countfile != "")){
				it = weights.find(query.getName());
				minCompare.weight = it->second;
			}
			else{	minCompare.weight = 1;	}
			
            
			printErrorData(minCompare, numParentSeqs, errorSummaryFile, errorSeqFile);
			
			if(!ignoreSeq){
				for(int i=0;i<minCompare.sequence.length();i++){
					char letter = minCompare.sequence[i];
					if(letter != 'r'){
						errorForward[letter][i] += minCompare.weight;
						errorReverse[letter][minCompare.total-i-1] += minCompare.weight;	
					}
				}                
			}
			
			if(aligned && qualFileName != "" && reportFileName != ""){
				report.read(reportFile);
				
				//				int origLength = report.getQueryLength();
				int startBase = report.getQueryStart();
				int endBase = report.getQueryEnd();
				
				quality.read(qualFile);
				
				if(!ignoreSeq){
					quality.updateQScoreErrorMap(qScoreErrorMap, minCompare.sequence, startBase, endBase, minCompare.weight);
					quality.updateForwardMap(qualForwardMap, startBase, endBase, minCompare.weight);
					quality.updateReverseMap(qualReverseMap, startBase, endBase, minCompare.weight);
				}
			}
            else if(aligned == false && qualFileName != ""){

                quality.read(qualFile);
                int qualityLength = quality.getLength();
                
                if(qualityLength != query.getNumBases()){   cout << "warning - quality and fasta sequence files do not match at " << query.getName() << '\t' << qualityLength <<'\t' << query.getNumBases() << endl;   }
                
                int startBase = 1;
                int endBase = qualityLength;

                if(!ignoreSeq){
					quality.updateQScoreErrorMap(qScoreErrorMap, minCompare.sequence, startBase, endBase, minCompare.weight);
					quality.updateForwardMap(qualForwardMap, startBase, endBase, minCompare.weight);
					quality.updateReverseMap(qualReverseMap, startBase, endBase, minCompare.weight);
				}
            }
            
			if(minCompare.errorRate <= threshold && !ignoreSeq){                
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
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = queryFile.tellg();
				if ((pos == -1) || (pos >= line.end)) { break; }
			#else
				if (queryFile.eof()) { break; }
			#endif
			
			if(index % 100 == 0){	m->mothurOutJustToScreen(toString(index)+"\n");	 }
		}
		queryFile.close();
		outChimeraReport.close();
        errorSummaryFile.close();	
		errorSeqFile.close();
        
        if(qFileName != "" && rFileName != "")      {   reportFile.close(); qualFile.close();   }
		else if(qFileName != "" && aligned == false){   qualFile.close();                       }
        
		//report progress
		m->mothurOutJustToScreen(toString(index)+"\n");	
		
		return index;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "driver");
		exit(1);
	}
}

//***************************************************************************************************************

void SeqErrorCommand::getReferences(){
	try {
		int numAmbigSeqs = 0;
		
		int maxStartPos = 0;
		int minEndPos = 100000;
		
		
        int start = time(NULL);
        
        ifstream referenceFile;
        m->openInputFile(referenceFileName, referenceFile);
        
        while(referenceFile){
            Sequence currentSeq(referenceFile);
            int numAmbigs = currentSeq.getAmbigBases();
            if(numAmbigs > 0){	numAmbigSeqs++;	}
            
            //			int startPos = currentSeq.getStartPos();
            //			if(startPos > maxStartPos)	{	maxStartPos = startPos;	}
            //
            //			int endPos = currentSeq.getEndPos();
            //			if(endPos < minEndPos)		{	minEndPos = endPos;		}
            if (currentSeq.getNumBases() == 0) {
                m->mothurOut("[WARNING]: " + currentSeq.getName() + " is blank, ignoring.");m->mothurOutEndLine();
            }else {
                referenceSeqs.push_back(currentSeq);
            }
            
            m->gobble(referenceFile);
        }
        referenceFile.close();
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " to read " + toString(referenceSeqs.size()) + " sequences.");m->mothurOutEndLine();
        
		
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
        Compare errors;
        
		if(query.getAlignLength() != reference.getAlignLength()){
			m->mothurOut("Warning: " + toString(query.getName()) + " and " + toString(reference.getName()) + " are different lengths\n");
		}
		int alignLength = query.getAlignLength();
	
		string q = query.getAligned();
		string r = reference.getAligned();

		int started = 0;
        
        errors.sequence = "";
		for(int i=0;i<alignLength;i++){

			if(q[i] != '.' && r[i] != '.' && (q[i] != '-' || r[i] != '-')){			//	no missing data and no double gaps
				if(r[i] != 'N'){
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
				else{
					
					if(q[i] == '-'){
						errors.sequence += 'd';	errors.total++;
					}						
					else{
						errors.sequence += 'r';
					}
				}
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
        if(errors.total != 0){  errors.errorRate = (double)(errors.total-errors.matches) / (double)errors.total;    }
        else{   errors.errorRate = 0;   }

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
    
    nameFile.close();
    
	return nameCountMap;
}

//***************************************************************************************************************

void SeqErrorCommand::printErrorHeader(ofstream& errorSummaryFile){
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

void SeqErrorCommand::printErrorData(Compare error, int numParentSeqs, ofstream& errorSummaryFile, ofstream& errorSeqFile){
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
        string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string subMatrixFileName = getOutputFileName("errormatrix",variables);
		ofstream subMatrixFile;
		m->openOutputFile(subMatrixFileName, subMatrixFile);
		outputNames.push_back(subMatrixFileName);  outputTypes["errormatrix"].push_back(subMatrixFileName);
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
        string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string errorForwardFileName = getOutputFileName("errorforward",variables);
		ofstream errorForwardFile;
		m->openOutputFile(errorForwardFileName, errorForwardFile);
		outputNames.push_back(errorForwardFileName);  outputTypes["errorforward"].push_back(errorForwardFileName);

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

		string errorReverseFileName = getOutputFileName("errorreverse",variables);
		ofstream errorReverseFile;
		m->openOutputFile(errorReverseFileName, errorReverseFile);
		outputNames.push_back(errorReverseFileName);  outputTypes["errorreverse"].push_back(errorReverseFileName);

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
        string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string errorQualityFileName = getOutputFileName("errorquality",variables);
		ofstream errorQualityFile;
		m->openOutputFile(errorQualityFileName, errorQualityFile);
		outputNames.push_back(errorQualityFileName);  outputTypes["errorquality"].push_back(errorQualityFileName);

		errorQualityFile << "qscore\tmatches\tsubstitutions\tinsertions\tambiguous" << endl;
		for(int i=0;i<101;i++){
			errorQualityFile << i << '\t' << qScoreErrorMap['m'][i] << '\t' << qScoreErrorMap['s'][i] << '\t' << qScoreErrorMap['i'][i] << '\t'<< qScoreErrorMap['a'][i] << endl;
		}
		errorQualityFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "printErrorQuality");
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
        string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string qualityForwardFileName = getOutputFileName("errorqualforward",variables);
		ofstream qualityForwardFile;
		m->openOutputFile(qualityForwardFileName, qualityForwardFile);
		outputNames.push_back(qualityForwardFileName);  outputTypes["errorqualforward"].push_back(qualityForwardFileName);

		for(int i=0;i<numColumns;i++){	qualityForwardFile << '\t' << i;	}	qualityForwardFile << endl;

		for(int i=0;i<numRows;i++){
			qualityForwardFile << i+1;
			for(int j=0;j<numColumns;j++){
				qualityForwardFile << '\t' << qualForwardMap[i][j];
			}

			qualityForwardFile << endl;
		}
		qualityForwardFile.close();

		
		string qualityReverseFileName = getOutputFileName("errorqualreverse",variables);
		ofstream qualityReverseFile;
		m->openOutputFile(qualityReverseFileName, qualityReverseFile);
		outputNames.push_back(qualityReverseFileName);  outputTypes["errorqualreverse"].push_back(qualityReverseFileName);
		
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
		m->errorOut(e, "SeqErrorCommand", "printQualityFR");
		exit(1);
	}
	
}

/**************************************************************************************************/

int SeqErrorCommand::setLines(string filename, string qfilename, string rfilename) {
	try {
        
        vector<unsigned long long> fastaFilePos;
        vector<unsigned long long> qfileFilePos;
        vector<unsigned long long> rfileFilePos;

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		//set file positions for fasta file
		fastaFilePos = m->divideFile(filename, processors);
		
		if (qfilename == "") { }
        else {
            //get name of first sequence in each chunk
            map<string, int> firstSeqNames;
            for (int i = 0; i < (fastaFilePos.size()-1); i++) {
                ifstream in;
                m->openInputFile(filename, in);
                in.seekg(fastaFilePos[i]);
                
                //adjust start if null strings
                if (i == 0) {  m->zapGremlins(in); m->gobble(in);  }
                
                Sequence temp(in);
                firstSeqNames[temp.getName()] = i;
                
                in.close();
            }
            
            //make copy to use below
            map<string, int> firstSeqNamesReport = firstSeqNames;
            
            
            if (qfilename != "") {
                //seach for filePos of each first name in the qfile and save in qfileFilePos
                ifstream inQual;
                m->openInputFile(qfilename, inQual);
                
                string input;
                while(!inQual.eof()){
                    input = m->getline(inQual);
                    
                    if (input.length() != 0) {
                        if(input[0] == '>'){ //this is a sequence name line
                            istringstream nameStream(input);
                            
                            string sname = "";  nameStream >> sname;
                            sname = sname.substr(1);
                            
                            m->checkName(sname);
                            
                            map<string, int>::iterator it = firstSeqNames.find(sname);
                            
                            if(it != firstSeqNames.end()) { //this is the start of a new chunk
                                unsigned long long pos = inQual.tellg();
                                qfileFilePos.push_back(pos - input.length() - 1);
                                firstSeqNames.erase(it);
                            }
                        }
                    }
                    
                    if (firstSeqNames.size() == 0) { break; }
                }
                inQual.close();
                
                if (firstSeqNames.size() != 0) {
                    for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                        m->mothurOut(it->first + " is in your fasta file and not in your quality file, aborting."); m->mothurOutEndLine();
                    }
                    m->setControl_pressed(true);
                    return processors;
                }
                
                //get last file position of qfile
                FILE * pFile;
                unsigned long long size;
                
                //get num bytes in file
                qfilename = m->getFullPathName(qfilename);
                pFile = fopen (qfilename.c_str(),"rb");
                if (pFile==NULL) perror ("Error opening file");
                else{
                    fseek (pFile, 0, SEEK_END);
                    size=ftell (pFile);
                    fclose (pFile);
                }
                
                qfileFilePos.push_back(size);
            }
            
            if(aligned){
                //seach for filePos of each first name in the rfile and save in rfileFilePos
                string junk, input;
                ifstream inR;
                
                m->openInputFile(rfilename, inR);
                
                //read column headers
                for (int i = 0; i < 16; i++) {
                    if (!inR.eof())	{	inR >> junk;	}
                    else			{	break;			}
                }
                
                while(!inR.eof()){
                    
                    input = m->getline(inR);
                    
                    if (input.length() != 0) {
                        
                        istringstream nameStream(input);
                        string sname = "";  nameStream >> sname;
                        
                        m->checkName(sname);
                        
                        map<string, int>::iterator it = firstSeqNamesReport.find(sname);
                        
                        if(it != firstSeqNamesReport.end()) { //this is the start of a new chunk
                            unsigned long long pos = inR.tellg();
                            rfileFilePos.push_back(pos - input.length() - 1);
                            firstSeqNamesReport.erase(it);
                        }
                    }
                    
                    if (firstSeqNamesReport.size() == 0) { break; }
                    m->gobble(inR);
                }
                inR.close();
                
                if (firstSeqNamesReport.size() != 0) {
                    for (map<string, int>::iterator it = firstSeqNamesReport.begin(); it != firstSeqNamesReport.end(); it++) {
                        m->mothurOut(it->first + " is in your fasta file and not in your report file, aborting."); m->mothurOutEndLine();
                    }
                    m->setControl_pressed(true);
                    return processors;
                }
                
                //get last file position of qfile
                FILE * rFile;
                unsigned long long sizeR;
                
                //get num bytes in file
                rfilename = m->getFullPathName(rfilename);
                rFile = fopen (rfilename.c_str(),"rb");
                if (rFile==NULL) perror ("Error opening file");
                else{
                    fseek (rFile, 0, SEEK_END);
                    sizeR=ftell (rFile);
                    fclose (rFile);
                }
                
                rfileFilePos.push_back(sizeR);
            }
        }
#else
		
		fastaFilePos.push_back(0); qfileFilePos.push_back(0);
		//get last file position of fastafile
		FILE * pFile;
		unsigned long long size;
		
		//get num bytes in file
        filename = m->getFullPathName(filename);
		pFile = fopen (filename.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
		fastaFilePos.push_back(size);
		
        if (qfilename != "") {
            //get last file position of qualfile
            FILE * qFile;
            
            //get num bytes in file
            qfilename = m->getFullPathName(qfilename);
            qFile = fopen (qfilename.c_str(),"rb");
            if (qFile==NULL) perror ("Error opening file");
            else{
                fseek (qFile, 0, SEEK_END);
                size=ftell (qFile);
                fclose (qFile);
            }
            qfileFilePos.push_back(size);
        }
        
         if (reportFileName != "" && aligned ) {
             rfileFilePos.push_back(0);
             
             //get last file position of qualfile
             FILE * rFile;
             
             //get num bytes in file
             rfilename = m->getFullPathName(rfilename);
             rFile = fopen (rfilename.c_str(),"rb");
             if (rFile==NULL) perror ("Error opening file");
             else{
                 fseek (rFile, 0, SEEK_END);
                 size=ftell (rFile);
                 fclose (rFile);
             }
             rfileFilePos.push_back(size);
         }
        
        processors = 1;
#endif
        
        if (m->getControl_pressed()) { return 0; }
        
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            lines.push_back(linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
            if (qualFileName != "") {  qLines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));  }
            if (reportFileName != "" && aligned ) {  rLines.push_back(linePair(rfileFilePos[i], rfileFilePos[(i+1)]));  }
        }
        if(qualFileName == "")	{	qLines = lines;	rLines = lines; } //fills with duds
        if(aligned == false){   rLines = lines; }
        
        
        return processors;

	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "setLines");
		exit(1);
	}
}

//***************************************************************************************************************
