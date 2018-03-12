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
		helpString += "The fasta parameter allows you to provide your dataset's fasta file.\n";
		helpString += "The reference parameter contains the Mock sequences.\n";
		helpString += "The qfile parameter ...\n";
		helpString += "The report parameter...\n";
		helpString += "The name parameter allows you to provide a name file associated with the fasta file.\n";
        helpString += "The count parameter allows you to provide a count file associated with the fasta file.\n";
		helpString += "The ignorechimeras parameter...\n";
		helpString += "The threshold parameter...\n";
		//helpString += "The processors parameter...\n";
		helpString += "Example seq.error(...).\n";
		;
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
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a names file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a names file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}

				it = parameters.find("qfile");
				//user has given a quality score file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("report");
				//user has given a alignment report file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["report"] = inputDir + it->second;		}
				}
				
			}
			//check for required parameters
			queryFileName = validParameter.validFile(parameters, "fasta");
			if (queryFileName == "not found") { 
				queryFileName = current->getFastaFile(); 
				if (queryFileName != "") { m->mothurOut("Using " + queryFileName + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fasta file and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (queryFileName == "not open") { queryFileName = ""; abort = true; }	
			else { current->setFastaFile(queryFileName); }
			
			referenceFileName = validParameter.validFile(parameters, "reference");
			if (referenceFileName == "not found") { m->mothurOut("reference is a required parameter for the seq.error command."); m->mothurOutEndLine(); abort = true; }
			else if (referenceFileName == "not open") { abort = true; }	
			
			//check for optional parameters
			namesFileName = validParameter.validFile(parameters, "name");
			if(namesFileName == "not found"){	namesFileName = "";	}
			else if (namesFileName == "not open") { namesFileName = ""; abort = true; }	
			else { current->setNameFile(namesFileName); }
            
            //check for optional parameters
			countfile = validParameter.validFile(parameters, "count");
			if(countfile == "not found"){	countfile = "";	}
			else if (countfile == "not open") { countfile = ""; abort = true; }
			else { current->setCountFile(countfile); }
			
			qualFileName = validParameter.validFile(parameters, "qfile");
			if(qualFileName == "not found"){	qualFileName = "";	}
			else if (qualFileName == "not open") { qualFileName = ""; abort = true; }	
			else { current->setQualFile(qualFileName); }
			
			reportFileName = validParameter.validFile(parameters, "report");
			if(reportFileName == "not found"){	reportFileName = "";	}
			else if (reportFileName == "not open") { reportFileName = ""; abort = true; }	
			
			outputDir = validParameter.valid(parameters, "outputdir");
			if (outputDir == "not found"){ //if user entered a file with a path then preserve it
				outputDir = util.hasPath(queryFileName); }
			
            if ((countfile != "") && (namesFileName != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			temp = validParameter.valid(parameters, "threshold");	if (temp == "not found") { temp = "1.00"; }
			util.mothurConvert(temp, threshold);
			
			temp = validParameter.valid(parameters, "ignorechimeras");	if (temp == "not found") { temp = "T"; }
			ignoreChimeras = util.isTrue(temp);
			
            temp = validParameter.valid(parameters, "aligned");			if (temp == "not found"){	temp = "t";				}
			aligned = util.isTrue(temp); 

			//temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			//processors = current->setProcessors(temp);
            processors=1;

			if ((namesFileName == "") && (queryFileName != "")){
				vector<string> files; files.push_back(queryFileName); 
				if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
		
		long start = time(NULL);
		maxLength = 5000;
		totalBases = 0;
		totalMatches = 0;
        vector<vector<int> > substitutionMatrix;
        vector<vector<int> > qualForwardMap;
        vector<vector<int> > qualReverseMap;
        vector<int> misMatchCounts;
        map<char, vector<int> > qScoreErrorMap;
        map<char, vector<int> > errorForward;
        map<char, vector<int> > errorReverse;
        vector<string> megaAlignVector;
        
        substitutionMatrix.resize(6);
        for(int i=0;i<6;i++){	substitutionMatrix[i].resize(6,0);	}
		
        string fileNameRoot = outputDir + util.getRootName(util.getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string errorSummaryFileName = getOutputFileName("errorsummary",variables);
		outputNames.push_back(errorSummaryFileName); outputTypes["errorsummary"].push_back(errorSummaryFileName);
			
		string errorSeqFileName = getOutputFileName("errorseq",variables);
		outputNames.push_back(errorSeqFileName); outputTypes["errorseq"].push_back(errorSeqFileName);
		
		string errorChimeraFileName = getOutputFileName("errorchimera",variables);
		outputNames.push_back(errorChimeraFileName); outputTypes["errorchimera"].push_back(errorChimeraFileName);
		
        vector<Sequence> referenceSeqs = getReferences(referenceFileName);	//read in reference sequences - make sure there's no ambiguous bases
        
        if (m->getControl_pressed()) { return 0; }
        
        long long numSeqs = createProcesses(queryFileName, qualFileName, reportFileName, errorSummaryFileName, errorSeqFileName, errorChimeraFileName, substitutionMatrix, qualForwardMap, qualReverseMap, misMatchCounts, qScoreErrorMap, errorForward, errorReverse, megaAlignVector, referenceSeqs);

		if(qualFileName != ""){		
			printErrorQuality(qScoreErrorMap);
			printQualityFR(qualForwardMap, qualReverseMap);
		}
		
		printErrorFRFile(errorForward, errorReverse);
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }

		string errorCountFileName = getOutputFileName("errorcount",variables);
		ofstream errorCountFile;
		util.openOutputFile(errorCountFileName, errorCountFile);
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
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }

		printSubMatrix(substitutionMatrix);
        
		string megAlignmentFileName = getOutputFileName("errorref-query",variables);
		ofstream megAlignmentFile;
		util.openOutputFile(megAlignmentFileName, megAlignmentFile);
		outputNames.push_back(megAlignmentFileName);  outputTypes["errorref-query"].push_back(megAlignmentFileName);
        
		for(int i=0;i<numRefs;i++){
			megAlignmentFile << referenceSeqs[i].getInlineSeq() << endl;
			megAlignmentFile << megaAlignVector[i] << endl;
		}
        megAlignmentFile.close();
      
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.\n");
		 
		m->mothurOut("\nOutput File Names: \n");
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]+"\n");  }
		m->mothurOutEndLine();
		
		return 0;	
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
struct seqErrorData {
    OutputWriter* summaryFile;
    OutputWriter* errorFile;
    OutputWriter* chimeraFile;
    vector<vector<int> > substitutionMatrix;
    vector<vector<int> > qualForwardMap;
    vector<vector<int> > qualReverseMap;
    vector<int> misMatchCounts;
    map<char, vector<int> > qScoreErrorMap;
    map<char, vector<int> > errorForward;
    map<char, vector<int> > errorReverse;
    vector<string> megaAlignVector;
    map<string, int> weights;
    vector<Sequence> referenceSeqs;
    
    string filename, qFileName, rFileName, chimeraHeader, errorHeader;
    bool aligned, ignoreChimeras, hasNameMap;
    double threshold;
    int maxLength, totalBases, totalMatches, maxMismatch;
    long long count;
    MothurOut* m;
    Utils util;
    
    seqErrorData(){  }
    seqErrorData(string fname, string qname, string rname, OutputWriter* sum,  OutputWriter* err, OutputWriter* chim, vector<Sequence> refseqs, bool al, int mxl, map<string, int> nam, double thre) {
        m = MothurOut::getInstance();
        weights = nam;
        hasNameMap = false;
        if (weights.size() > 0) { hasNameMap = true;  }
        maxLength = mxl;
        filename = fname;
        qFileName = qname;
        rFileName = rname;
        summaryFile = sum;
        errorFile = err;
        chimeraFile = chim;
        referenceSeqs = refseqs;
        aligned = al;
        threshold = thre;
        chimeraHeader = "";
        errorHeader = "";
        count = 0; totalBases = 0; totalMatches = 0; maxMismatch = 0;
        substitutionMatrix.resize(6);
        for(int i=0;i<6;i++){	substitutionMatrix[i].resize(6,0);	}
    }
};
//***************************************************************************************************************

Compare getErrors(Sequence query, Sequence reference, MothurOut* m){
    try {
        Compare errors;
        
        if(query.getAlignLength() != reference.getAlignLength()){ m->mothurOut("[WARNING]: " + toString(query.getName()) + " and " + toString(reference.getName()) + " are different lengths\n"); }
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

void printErrorData(Compare error, int numParentSeqs, OutputWriter* summaryFile, OutputWriter* errorFile, vector<vector<int> >& substitutionMatrix, bool ignoreChimeras){
    string summaryOutput = error.queryName + '\t' + error.refName + '\t' + toString(error.weight) + '\t';
    summaryOutput += toString(error.AA)  + '\t' +  toString(error.AT) + '\t' + toString(error.AG)  + '\t' +  toString(error.AC) + '\t';
    summaryOutput += toString(error.TA)  + '\t' +  toString(error.TT) + '\t' + toString(error.TG)  + '\t' +  toString(error.TC) + '\t';
    summaryOutput += toString(error.GA)  + '\t' +  toString(error.GT) + '\t' + toString(error.GG)  + '\t' +  toString(error.GC) + '\t';
    summaryOutput += toString(error.CA)  + '\t' +  toString(error.CT) + '\t' + toString(error.CG)  + '\t' +  toString(error.CC) + '\t';
    summaryOutput += toString(error.NA)  + '\t' +  toString(error.NT) + '\t' + toString(error.NG)  + '\t' +  toString(error.NC) + '\t';
    summaryOutput += toString(error.Ai)  + '\t' +  toString(error.Ti) + '\t' + toString(error.Gi)  + '\t' +  toString(error.Ci)  + '\t' +  toString(error.Ni) + '\t';
    summaryOutput += toString(error.dA)  + '\t' +  toString(error.dT) + '\t' + toString(error.dG)  + '\t' +  toString(error.dC) + '\t';
    
    summaryOutput += toString(error.Ai + error.Ti + error.Gi + error.Ci) + '\t';			//insertions
    summaryOutput += toString(error.dA + error.dT + error.dG + error.dC) + '\t';			//deletions
    summaryOutput += toString(error.mismatches - (error.Ai + error.Ti + error.Gi + error.Ci) - (error.dA + error.dT + error.dG + error.dC) - (error.NA + error.NT + error.NG + error.NC + error.Ni)) + '\t';	//substitutions
    summaryOutput += toString(error.NA + error.NT + error.NG + error.NC + error.Ni) + '\t';	//ambiguities
    summaryOutput += toString(error.matches) + '\t' + toString(error.mismatches)  + '\t' + toString(error.total) + '\t' + toString(error.errorRate)  + '\t' + toString(numParentSeqs) + "\n";
    summaryFile->write(summaryOutput);
    
    summaryOutput = '>' + error.queryName + "\tref:" + error.refName + '\n' + error.sequence + '\n';
    errorFile->write(summaryOutput);
    
    int a=0;		int t=1;		int g=2;		int c=3;
    int gap=4;		int n=5;
    
    if(numParentSeqs == 1 || !ignoreChimeras){
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
/**********************************************************************************************************************/
void driverSeqError(seqErrorData* params) {
    try {
        ReportFile report;
        QualityScores quality;
        
        params->misMatchCounts.resize(11, 0);
        map<string, int>::iterator it;
        params->qScoreErrorMap['m'].assign(101, 0);
        params->qScoreErrorMap['s'].assign(101, 0);
        params->qScoreErrorMap['i'].assign(101, 0);
        params->qScoreErrorMap['a'].assign(101, 0);
        
        params->errorForward['m'].assign(params->maxLength,0);
        params->errorForward['s'].assign(params->maxLength,0);
        params->errorForward['i'].assign(params->maxLength,0);
        params->errorForward['d'].assign(params->maxLength,0);
        params->errorForward['a'].assign(params->maxLength,0);
        
        params->errorReverse['m'].assign(params->maxLength,0);
        params->errorReverse['s'].assign(params->maxLength,0);
        params->errorReverse['i'].assign(params->maxLength,0);
        params->errorReverse['d'].assign(params->maxLength,0);
        params->errorReverse['a'].assign(params->maxLength,0);
        
        //open inputfiles and go to beginning place for this processor
        ifstream queryFile;
        params->util.openInputFile(params->filename, queryFile);
        
        
        ifstream reportFile; ifstream qualFile;
        if((params->qFileName != "" && params->rFileName != "" && params->aligned)){
            params->util.openInputFile(params->qFileName, qualFile);
            
            //gobble headers
            report.readHeaders(reportFile, params->rFileName);
            
            params->qualForwardMap.resize(params->maxLength);
            params->qualReverseMap.resize(params->maxLength);
            for(int i=0;i<params->maxLength;i++){
                params->qualForwardMap[i].assign(101,0);
                params->qualReverseMap[i].assign(101,0);
            }
        }
        else if(params->qFileName != "" && !params->aligned){
            
            params->util.openInputFile(params->qFileName, qualFile);
            
            params->qualForwardMap.resize(params->maxLength);
            params->qualReverseMap.resize(params->maxLength);
            for(int i=0;i<params->maxLength;i++){
                params->qualForwardMap[i].assign(101,0);
                params->qualReverseMap[i].assign(101,0);
            }
        }
        
        RefChimeraTest chimeraTest = RefChimeraTest(params->referenceSeqs, params->aligned);
        params->chimeraHeader = chimeraTest.getHeader();
        
        int numRefs = params->referenceSeqs.size();
        params->megaAlignVector.assign(numRefs, "");
        
        int index = 0;
        bool ignoreSeq = 0;
        
        bool moreSeqs = 1;
        while (moreSeqs) {
            
            Sequence query(queryFile);
            int numParentSeqs = -1;
            int closestRefIndex = -1;
            
            string querySeq = query.getAligned();
            if (!params->aligned) {  querySeq = query.getUnaligned();  }
            
            string output = "";
            numParentSeqs = chimeraTest.analyzeQuery(query.getName(), querySeq, output);
            params->chimeraFile->write(output);
            
            closestRefIndex = chimeraTest.getClosestRefIndex();
            
            Sequence reference = params->referenceSeqs[closestRefIndex];
            
            reference.setAligned(chimeraTest.getClosestRefAlignment());
            query.setAligned(chimeraTest.getQueryAlignment());
            
            if(numParentSeqs > 1 && params->ignoreChimeras == 1)	{	ignoreSeq = 1;	}
            else                                                    {	ignoreSeq = 0;	}
            
            Compare minCompare = getErrors(query, reference, params->m);
            
            if(params->hasNameMap){
                it = params->weights.find(query.getName());
                minCompare.weight = it->second;
            }
            else{	minCompare.weight = 1;	}
            
            printErrorData(minCompare, numParentSeqs, params->summaryFile, params->errorFile, params->substitutionMatrix, params->ignoreChimeras);
            
            if(!ignoreSeq){
                for(int i=0;i<minCompare.sequence.length();i++){
                    char letter = minCompare.sequence[i];
                    if(letter != 'r'){
                        params->errorForward[letter][i] += minCompare.weight;
                        params->errorReverse[letter][minCompare.total-i-1] += minCompare.weight;
                    }
                }
            }
            
            if(params->aligned && params->qFileName != "" && params->rFileName != ""){
                report.read(reportFile);
                
                int startBase = report.getQueryStart();
                int endBase = report.getQueryEnd();
                
                quality.read(qualFile);
                
                if(!ignoreSeq){
                    quality.updateQScoreErrorMap(params->qScoreErrorMap, minCompare.sequence, startBase, endBase, minCompare.weight);
                    quality.updateForwardMap(params->qualForwardMap, startBase, endBase, minCompare.weight);
                    quality.updateReverseMap(params->qualReverseMap, startBase, endBase, minCompare.weight);
                }
            }
            else if(!params->aligned && params->qFileName != ""){
                
                quality.read(qualFile);
                int qualityLength = quality.getLength();
                
                if(qualityLength != query.getNumBases()){   params->m->mothurOut("[WARNING]: - quality and fasta sequence files do not match at " + query.getName() + '\t' + toString(qualityLength) + '\t' + toString(query.getNumBases()) +'\n');  }
                
                int startBase = 1;
                int endBase = qualityLength;
                
                if(!ignoreSeq){
                    quality.updateQScoreErrorMap(params->qScoreErrorMap, minCompare.sequence, startBase, endBase, minCompare.weight);
                    quality.updateForwardMap(params->qualForwardMap, startBase, endBase, minCompare.weight);
                    quality.updateReverseMap(params->qualReverseMap, startBase, endBase, minCompare.weight);
                }
            }
            
            if(minCompare.errorRate <= params->threshold && !ignoreSeq){
                params->totalBases += (minCompare.total * minCompare.weight);
                params->totalMatches += minCompare.matches * minCompare.weight;
                if(minCompare.mismatches > params->maxMismatch){
                    params->maxMismatch = minCompare.mismatches;
                    params->misMatchCounts.resize(params->maxMismatch + 1, 0);
                }				
                params->misMatchCounts[minCompare.mismatches] += minCompare.weight;
                params->megaAlignVector[closestRefIndex] += query.getInlineSeq() + '\n';
            }
            
            index++;

            if (queryFile.eof()) { break; }
            
            if(index % 100 == 0){	params->m->mothurOutJustToScreen(toString(index)+"\n");	 }
        }
        queryFile.close();
        if(params->qFileName != "" && params->rFileName != "")      {   reportFile.close(); qualFile.close();   }
        else if(params->qFileName != "" && params->aligned == false){   qualFile.close();                       }
        
        //report progress
        params->m->mothurOutJustToScreen(toString(index)+"\n");
        
        params->count = index;
    }
    catch(exception& e) {
        params->m->errorOut(e, "SeqErrorCommand", "driver");
        exit(1);
    }
}
//**********************************************************************************************************************
long long SeqErrorCommand::createProcesses(string filename, string qFileName, string rFileName, string summaryFileName, string errorOutputFileName, string chimeraOutputFileName, vector<vector<int> >& substitutionMatrix, vector<vector<int> >& qualForwardMap, vector<vector<int> >& qualReverseMap, vector<int>& misMatchCounts, map<char, vector<int> >& qScoreErrorMap, map<char, vector<int> >& errorForward, map<char, vector<int> >& errorReverse, vector<string>& megaAlignVector, vector<Sequence>& referenceSeqs) {
	try {
        map<string, int> weights;
        if(namesFileName != "")     {	weights = util.readNames(namesFileName);                                            }
        else if (countfile != "")   {   CountTable ct;  ct.readTable(countfile, false, false);  weights = ct.getNameMap();  }
        
        ofstream out;
        util.openOutputFile(summaryFileName, out);
        out << "query\treference\tweight\tAA\tAT\tAG\tAC\tTA\tTT\tTG\tTC\tGA\tGT\tGG\tGC\tCA\tCT\tCG\tCC\tNA\tNT\tNG\tNC\tAi\tTi\tGi\tCi\tNi\tdA\tdT\tdG\tdC\tinsertions\tdeletions\tsubstitutions\tambig\tmatches\tmismatches\ttotal\terror\tnumparents\n";
        out.close();

        auto synchronizedOutputFastaTrimFile = std::make_shared<SynchronizedOutputFile>(summaryFileName, true); //append to add headers
        auto synchronizedOutputErrorFile = std::make_shared<SynchronizedOutputFile>(errorOutputFileName); synchronizedOutputErrorFile->setFixedShowPoint(); synchronizedOutputErrorFile->setPrecision(6);
        auto synchronizedOutputChimeraFile = std::make_shared<SynchronizedOutputFile>(chimeraOutputFileName);
        
        OutputWriter* threadSummaryWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
        OutputWriter* threadErrorWriter = new OutputWriter(synchronizedOutputErrorFile);
        OutputWriter* threadChimeraWriter = new OutputWriter(synchronizedOutputChimeraFile);
        
        seqErrorData* dataBundle = new seqErrorData(filename, qFileName, rFileName, threadSummaryWriter, threadErrorWriter, threadChimeraWriter, referenceSeqs, aligned, maxLength, weights, threshold);
        driverSeqError(dataBundle);


        totalMatches = dataBundle->totalMatches;
        
        totalBases = dataBundle->totalBases;
       
        int maxMisMatch = dataBundle->maxMismatch;
       
        long long numSeqs = dataBundle->count;
     
        substitutionMatrix = dataBundle->substitutionMatrix;
      
        qualForwardMap = dataBundle->qualForwardMap;
  
        qualReverseMap = dataBundle->qualReverseMap;
   
        misMatchCounts = dataBundle->misMatchCounts;
  
        qScoreErrorMap = dataBundle->qScoreErrorMap;
    
        errorForward = dataBundle->errorForward;

        errorReverse = dataBundle->errorReverse;

        megaAlignVector = dataBundle->megaAlignVector;

        delete threadErrorWriter; delete threadChimeraWriter; delete threadSummaryWriter;
        delete dataBundle;

		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "createProcesses");
		exit(1);
	}
}

//***************************************************************************************************************
vector<Sequence> SeqErrorCommand::getReferences(string refFileName){
	try {
		int numAmbigSeqs = 0;
		int maxStartPos = 0;
		int minEndPos = 100000;
        
        ifstream referenceFile;
        util.openInputFile(refFileName, referenceFile);
        
        vector<Sequence> referenceSeqs;
        while(!referenceFile.eof()){
            if (m->getControl_pressed()) { break; }
            
            Sequence currentSeq(referenceFile);
            int numAmbigs = currentSeq.getAmbigBases();
            if(numAmbigs > 0){	numAmbigSeqs++;	}
            
            			int startPos = currentSeq.getStartPos();
            			if(startPos > maxStartPos)	{	maxStartPos = startPos;	}
            
            			int endPos = currentSeq.getEndPos();
            			if(endPos < minEndPos)		{	minEndPos = endPos;		}
            if (currentSeq.getNumBases() == 0) {
                m->mothurOut("[WARNING]: " + currentSeq.getName() + " is blank, ignoring.");m->mothurOutEndLine();
            }else {
                referenceSeqs.push_back(currentSeq);
            }
            
            util.gobble(referenceFile);
        }
        referenceFile.close();
        
		numRefs = referenceSeqs.size();
		
		for(int i=0;i<numRefs;i++){
			referenceSeqs[i].padToPos(maxStartPos);
			referenceSeqs[i].padFromPos(minEndPos);
        }
		
		if(numAmbigSeqs != 0){ m->mothurOut("Warning: " + toString(numAmbigSeqs) + " reference sequences have ambiguous bases, these bases will be ignored\n"); }
		
        return referenceSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorCommand", "getReferences");
		exit(1);
	}
}
//***************************************************************************************************************

void SeqErrorCommand::printSubMatrix(vector<vector<int> >& substitutionMatrix){
	try {
        string fileNameRoot = outputDir + util.getRootName(util.getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string subMatrixFileName = getOutputFileName("errormatrix",variables);
		ofstream subMatrixFile;
		util.openOutputFile(subMatrixFileName, subMatrixFile);
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

void SeqErrorCommand::printErrorFRFile(map<char, vector<int> >& errorForward, map<char, vector<int> >& errorReverse){
	try{
        string fileNameRoot = outputDir + util.getRootName(util.getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string errorForwardFileName = getOutputFileName("errorforward",variables);
		ofstream errorForwardFile;
		util.openOutputFile(errorForwardFileName, errorForwardFile);
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
		util.openOutputFile(errorReverseFileName, errorReverseFile);
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

void SeqErrorCommand::printErrorQuality(map<char, vector<int> >& qScoreErrorMap){
	try{
        string fileNameRoot = outputDir + util.getRootName(util.getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string errorQualityFileName = getOutputFileName("errorquality",variables);
		ofstream errorQualityFile;
		util.openOutputFile(errorQualityFileName, errorQualityFile);
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

void SeqErrorCommand::printQualityFR(vector<vector<int> >& qualForwardMap, vector<vector<int> >& qualReverseMap){

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
        string fileNameRoot = outputDir + util.getRootName(util.getSimpleName(queryFileName));
        map<string, string> variables; 
		variables["[filename]"] = fileNameRoot;
		string qualityForwardFileName = getOutputFileName("errorqualforward",variables);
		ofstream qualityForwardFile;
		util.openOutputFile(qualityForwardFileName, qualityForwardFile);
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
		util.openOutputFile(qualityReverseFileName, qualityReverseFile);
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

int SeqErrorCommand::setLines(string filename, string qfilename, string rfilename, vector<linePair>& lines, vector<linePair>& qLines, vector<linePair>& rLines) {
	try {
        
        vector<unsigned long long> fastaFilePos;
        vector<unsigned long long> qfileFilePos;
        vector<unsigned long long> rfileFilePos;

#if defined NON_WINDOWS
		//set file positions for fasta file
		fastaFilePos = util.divideFile(filename, processors);
		
		if (qfilename == "") { }
        else {
            //get name of first sequence in each chunk
            map<string, int> firstSeqNames;
            for (int i = 0; i < (fastaFilePos.size()-1); i++) {
                ifstream in;
                util.openInputFile(filename, in);
                in.seekg(fastaFilePos[i]);
                
                //adjust start if null strings
                if (i == 0) {  util.zapGremlins(in); util.gobble(in);  }
                
                Sequence temp(in);
                firstSeqNames[temp.getName()] = i;
                
                in.close();
            }
            
            //make copy to use below
            map<string, int> firstSeqNamesReport = firstSeqNames;
            
            
            if (qfilename != "") {
                //seach for filePos of each first name in the qfile and save in qfileFilePos
                ifstream inQual;
                util.openInputFile(qfilename, inQual);
                
                string input;
                while(!inQual.eof()){
                    input = util.getline(inQual);
                    
                    if (input.length() != 0) {
                        if(input[0] == '>'){ //this is a sequence name line
                            istringstream nameStream(input);
                            
                            string sname = "";  nameStream >> sname;
                            sname = sname.substr(1);
                            
                            util.checkName(sname);
                            
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
                qfilename = util.getFullPathName(qfilename);
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
                
                util.openInputFile(rfilename, inR);
                
                //read column headers
                for (int i = 0; i < 16; i++) {
                    if (!inR.eof())	{	inR >> junk;	}
                    else			{	break;			}
                }
                
                while(!inR.eof()){
                    
                    input = util.getline(inR);
                    
                    if (input.length() != 0) {
                        
                        istringstream nameStream(input);
                        string sname = "";  nameStream >> sname;
                        
                        util.checkName(sname);
                        
                        map<string, int>::iterator it = firstSeqNamesReport.find(sname);
                        
                        if(it != firstSeqNamesReport.end()) { //this is the start of a new chunk
                            unsigned long long pos = inR.tellg();
                            rfileFilePos.push_back(pos - input.length() - 1);
                            firstSeqNamesReport.erase(it);
                        }
                    }
                    
                    if (firstSeqNamesReport.size() == 0) { break; }
                    util.gobble(inR);
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
                rfilename = util.getFullPathName(rfilename);
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
        filename = util.getFullPathName(filename);
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
            qfilename = util.getFullPathName(qfilename);
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
             rfilename = util.getFullPathName(rfilename);
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
