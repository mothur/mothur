//
//  primerdesigncommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/18/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "primerdesigncommand.h"

//**********************************************************************************************************************
vector<string> PrimerDesignCommand::setParameters(){	
	try {
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","summary-list",false,true,true); parameters.push_back(plist);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","",false,true, true); parameters.push_back(pfasta);
 		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pcount);
        CommandParameter plength("length", "Number", "", "18", "", "", "","",false,false); parameters.push_back(plength);
        CommandParameter pmintm("mintm", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmintm);
        CommandParameter pmaxtm("maxtm", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxtm);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pprocessors);
        CommandParameter potunumber("otunumber", "Number", "", "-1", "", "", "","",false,true,true); parameters.push_back(potunumber);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
        CommandParameter pcutoff("cutoff", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PrimerDesignCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The primer.design allows you to identify sequence fragments that are specific to particular OTUs.\n";
		helpString += "The primer.design command parameters are: list, fasta, name, count, otunumber, cutoff, length, pdiffs, mintm, maxtm, processors and label.\n";
		helpString += "The list parameter allows you to provide a list file and is required.\n";
        helpString += "The fasta parameter allows you to provide a fasta file and is required.\n";
        helpString += "The name parameter allows you to provide a name file associated with your fasta file.\n";
        helpString += "The count parameter allows you to provide a count file associated with your fasta file.\n";
        helpString += "The label parameter is used to indicate the label you want to use from your list file.\n";
        helpString += "The otunumber parameter is used to indicate the otu you want to use from your list file. It is required.\n";
        helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The length parameter is used to indicate the length of the primer. The default is 18.\n";
        helpString += "The mintm parameter is used to indicate minimum melting temperature.\n";
        helpString += "The maxtm parameter is used to indicate maximum melting temperature.\n";
        helpString += "The processors parameter allows you to indicate the number of processors you want to use. Default=1.\n";
        helpString += "The cutoff parameter allows you set a percentage of sequences that support the base. For example: cutoff=97 would only return a sequence that only showed ambiguities for bases that were not supported by at least 97% of sequences.\n";
		helpString += "The primer.desing command should be in the following format: primer.design(list=yourListFile, fasta=yourFastaFile, name=yourNameFile)\n";
		helpString += "primer.design(list=final.an.list, fasta=final.fasta, name=final.names, label=0.03)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string PrimerDesignCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],[distance],otu.cons.fasta"; } 
        else if (type == "summary") {  pattern = "[filename],[distance],primer.summary"; }
        else if (type == "list") {  pattern = "[filename],pick,[extension]"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PrimerDesignCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
PrimerDesignCommand::PrimerDesignCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames; 
        outputTypes["fasta"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
        
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "PrimerDesignCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
PrimerDesignCommand::PrimerDesignCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames; 
            outputTypes["fasta"] = tempOutNames;
            outputTypes["list"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
 				string path;
				it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
                
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
            }
                        
			//check for parameters
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { m->setCountTableFile(countfile); }
            
            //get fastafile - it is required
            fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort=true;  }
			else if (fastafile == "not found") {  
                fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
            }else  { m->setFastaFile(fastafile); }
            
            //get listfile - it is required
            listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort=true;  }
			else if (listfile == "not found") {  
                listfile = m->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current listfile and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
            }else  { m->setListFile(listfile); }

            
			if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = m->hasPath(listfile); //if user entered a file with a path then preserve it	
			}
            
            string temp = validParameter.validFile(parameters, "cutoff", false);  if (temp == "not found") { temp = "100"; }
			m->mothurConvert(temp, cutoff); 
            
            temp = validParameter.validFile(parameters, "pdiffs", false);  if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, pdiffs); 
            
            temp = validParameter.validFile(parameters, "length", false);  if (temp == "not found") { temp = "18"; }
			m->mothurConvert(temp, length); 
            
            temp = validParameter.validFile(parameters, "mintm", false);  if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, minTM); 
            
            temp = validParameter.validFile(parameters, "maxtm", false);  if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxTM); 
            
            temp = validParameter.validFile(parameters, "otunumber", false);  if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, otunumber); 
            if (otunumber < 1) {  m->mothurOut("[ERROR]: You must provide an OTU number, aborting.\n"); abort = true; }
            
            temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
            label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }
        
            if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    parser.getNameFile(files);
                }
            }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "PrimerDesignCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int PrimerDesignCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        int start = time(NULL);
        //////////////////////////////////////////////////////////////////////////////
        //              get file inputs                                             //
        //////////////////////////////////////////////////////////////////////////////
        
        //reads list file and selects the label the users specified or the first label
        getListVector();
        if (otunumber > list->getNumBins()) { m->mothurOut("[ERROR]: You selected an OTU number larger than the number of OTUs you have in your list file, quitting.\n"); return 0; }
        
        map<string, int> nameMap;
        unsigned long int numSeqs;  //used to sanity check the files. numSeqs = total seqs for namefile and uniques for count.
                                    //list file should have all seqs if namefile was used to create it and only uniques in count file was used.
        
        if (namefile != "")         {  nameMap = m->readNames(namefile, numSeqs);       }
        else if (countfile != "")   {  nameMap = readCount(numSeqs);                    }
        else { numSeqs = list->getNumSeqs();  }
        
        //sanity check
        if (numSeqs != list->getNumSeqs()) {
            if (namefile != "")         {  m->mothurOut("[ERROR]: Your list file contains " + toString(list->getNumSeqs()) + " sequences, and your name file contains " + toString(numSeqs) + " sequences, aborting. Do you have the correct files? Perhaps you forgot to include the name file when you clustered? \n");   }
            else if (countfile != "") {
                m->mothurOut("[ERROR]: Your list file contains " + toString(list->getNumSeqs()) + " sequences, and your count file contains " + toString(numSeqs) + " unique sequences, aborting. Do you have the correct files? Perhaps you forgot to include the count file when you clustered? \n");  
            }
            m->control_pressed = true;
        }
        
        if (m->control_pressed) { delete list; return 0; }
        
        //////////////////////////////////////////////////////////////////////////////
        //              process data                                                //
        //////////////////////////////////////////////////////////////////////////////
        m->mothurOut("\nFinding consensus sequences for each otu..."); cout.flush();
        
        vector<Sequence> conSeqs = createProcessesConSeqs(nameMap, numSeqs);
        
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listfile));
        variables["[distance]"] = list->getLabel();
        string consFastaFile = getOutputFileName("fasta", variables);
        outputNames.push_back(consFastaFile); outputTypes["fasta"].push_back(consFastaFile);
        ofstream out;
        m->openOutputFile(consFastaFile, out);
        for (int i = 0; i < conSeqs.size(); i++) {  conSeqs[i].printSequence(out);  }
        out.close();
        
        m->mothurOut("Done.\n\n");
        
        set<string> primers = getPrimer(conSeqs[otunumber-1]);  
        
        if (m->control_pressed) { delete list; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
        string consSummaryFile = getOutputFileName("summary", variables);
        outputNames.push_back(consSummaryFile); outputTypes["summary"].push_back(consSummaryFile);
        ofstream outSum;
        m->openOutputFile(consSummaryFile, outSum);
        
        outSum << "PrimerOtu: " << otunumber << " Members: " << list->get(otunumber-1) << endl << "Primers\tminTm\tmaxTm" << endl;
        
        //find min and max melting points
        vector<double> minTms;
        vector<double> maxTms;
        string primerString = "";
        for (set<string>::iterator it = primers.begin(); it != primers.end();) {
            
            double minTm, maxTm;
            findMeltingPoint(*it, minTm, maxTm);
            if ((minTM == -1) && (maxTM == -1)) { //user did not set min or max Tm so save this primer
                minTms.push_back(minTm);
                maxTms.push_back(maxTm);
                outSum << *it << '\t' << minTm << '\t' << maxTm << endl;
                it++;
            }else if ((minTM == -1) && (maxTm <= maxTM)){ //user set max and no min, keep if below max
                minTms.push_back(minTm);
                maxTms.push_back(maxTm);
                outSum << *it << '\t' << minTm << '\t' << maxTm << endl;
                it++;
            }else if ((maxTM == -1) && (minTm >= minTM)){ //user set min and no max, keep if above min
                minTms.push_back(minTm);
                maxTms.push_back(maxTm);
                outSum << *it << '\t' << minTm << '\t' << maxTm << endl;
                it++;
            }else if ((maxTm <= maxTM) && (minTm >= minTM)) { //keep if above min and below max
                minTms.push_back(minTm);
                maxTms.push_back(maxTm);
                outSum << *it << '\t' << minTm << '\t' << maxTm << endl;
                it++;
            }else { primers.erase(it++);  } //erase because it didn't qualify
        }
        
        outSum << "\nOTUNumber\tPrimer\tStart\tEnd\tLength\tMismatches\tminTm\tmaxTm\n";
        outSum.close();
        
        //check each otu's conseq for each primer in otunumber
        set<int> otuToRemove = createProcesses(consSummaryFile, minTms, maxTms, primers, conSeqs);
        
        if (m->control_pressed) { delete list; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
        //print new list file
        map<string, string> mvariables; 
        mvariables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listfile));
        mvariables["[extension]"] = m->getExtension(listfile);
        string newListFile = getOutputFileName("list", mvariables);
        outputNames.push_back(newListFile); outputTypes["list"].push_back(newListFile);
        ofstream outList;
        m->openOutputFile(newListFile, outList);
        
        outList << list->getLabel() << '\t' << (list->getNumBins()-otuToRemove.size()) << '\t';
        for (int j = 0; j < list->getNumBins(); j++) {
            if (m->control_pressed) { break; }
            //good otus
            if (otuToRemove.count(j) == 0) {  
                string bin = list->get(j);
                if (bin != "") {  outList << bin << '\t';  } 
            }
        }
        outList << endl;
        outList.close();
        
        if (m->control_pressed) { delete list; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
        delete list;
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process " + toString(list->getNumBins()) + " OTUs.\n");
        
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "execute");
		exit(1);
	}
}
//********************************************************************/
//used http://www.biophp.org/minitools/melting_temperature/ as a reference to substitute degenerate bases 
// in order to find the min and max Tm values.
//Tm =  64.9°C + 41°C x (number of G’s and C’s in the primer – 16.4)/N

/* A = adenine
 * C = cytosine
 * G = guanine
 * T = thymine
 * R = G A (purine)
 * Y = T C (pyrimidine)
 * K = G T (keto)
 * M = A C (amino)
 * S = G C (strong bonds)
 * W = A T (weak bonds)
 * B = G T C (all but A)
 * D = G A T (all but C)
 * H = A C T (all but G)
 * V = G C A (all but T)
 * N = A G C T (any) */

int PrimerDesignCommand::findMeltingPoint(string primer, double& minTm, double& maxTm){
    try {
        string minTmprimer = primer;
        string maxTmprimer = primer;
        
        //find minimum Tm string substituting for degenerate bases
        for (int i = 0; i < minTmprimer.length(); i++) {
            minTmprimer[i] = toupper(minTmprimer[i]);
            
            if (minTmprimer[i] == 'Y') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'R') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'W') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'K') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'M') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'D') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'V') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'H') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'B') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'N') { minTmprimer[i] = 'A'; }
            else if (minTmprimer[i] == 'S') { minTmprimer[i] = 'G'; }
        }
        
        //find maximum Tm string substituting for degenerate bases
        for (int i = 0; i < maxTmprimer.length(); i++) {
            maxTmprimer[i] = toupper(maxTmprimer[i]);
            
            if (maxTmprimer[i] == 'Y') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'R') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'W') { maxTmprimer[i] = 'A'; }
            else if (maxTmprimer[i] == 'K') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'M') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'D') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'V') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'H') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'B') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'N') { maxTmprimer[i] = 'G'; }
            else if (maxTmprimer[i] == 'S') { maxTmprimer[i] = 'G'; }
        }
        
        int numGC = 0;
        for (int i = 0; i < minTmprimer.length(); i++) {
            if (minTmprimer[i] == 'G')       { numGC++; }
            else if (minTmprimer[i] == 'C')  { numGC++; }
        }
        
        minTm = 64.9 + 41 * (numGC - 16.4) / (double) minTmprimer.length();
        
        numGC = 0;
        for (int i = 0; i < maxTmprimer.length(); i++) {
            if (maxTmprimer[i] == 'G')       { numGC++; }
            else if (maxTmprimer[i] == 'C')  { numGC++; }
        }
        
        maxTm = 64.9 + 41 * (numGC - 16.4) / (double) maxTmprimer.length();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "findMeltingPoint");
		exit(1);
	}
}
//********************************************************************/
//search for a primer over the sequence string
bool PrimerDesignCommand::findPrimer(string rawSequence, string primer, vector<int>& primerStart, vector<int>& primerEnd, vector<int>& mismatches){
	try {
        bool foundAtLeastOne = false;  //innocent til proven guilty
        
        //look for exact match
        if(rawSequence.length() < primer.length()) {  return false;  }
			
        //search for primer
        for (int j = 0; j < rawSequence.length()-length; j++){
            
            if (m->control_pressed) {  return foundAtLeastOne; }
            
            string rawChunk = rawSequence.substr(j, length);
            
            int numDiff = countDiffs(primer, rawChunk);
           
            if(numDiff <= pdiffs){
                primerStart.push_back(j);
                primerEnd.push_back(j+length);
                mismatches.push_back(numDiff);
                foundAtLeastOne = true;
            }
        }
		
		return foundAtLeastOne;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "findPrimer");
		exit(1);
	}
}
//********************************************************************/
//find all primers for the given sequence
set<string> PrimerDesignCommand::getPrimer(Sequence primerSeq){
	try {
        set<string> primers;
        
        string rawSequence = primerSeq.getUnaligned();
        
        for (int j = 0; j < rawSequence.length()-length; j++){
            if (m->control_pressed) { break; }
            
            string primer = rawSequence.substr(j, length);
            primers.insert(primer);
        }
        
        return primers;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "getPrimer");
		exit(1);
	}
}
/**************************************************************************************************/
set<int> PrimerDesignCommand::createProcesses(string newSummaryFile, vector<double>& minTms, vector<double>& maxTms, set<string>& primers, vector<Sequence>& conSeqs) {
	try {
		
		vector<int> processIDS;
		int process = 1;
        set<int> otusToRemove;
        int numBinsProcessed = 0;
		
		//sanity check
        int numBins = conSeqs.size();
		if (numBins < processors) { processors = numBins; }
		
		//divide the otus between the processors
		vector<linePair> lines;
		int numOtusPerProcessor = numBins / processors;
		for (int i = 0; i < processors; i++) {
			int startIndex =  i * numOtusPerProcessor;
			int endIndex = (i+1) * numOtusPerProcessor;
			if(i == (processors - 1)){	endIndex = numBins; 	}
			lines.push_back(linePair(startIndex, endIndex));
		}
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
                //clear old file because we append in driver
                m->mothurRemove(newSummaryFile + toString(getpid()) + ".temp");
                
				otusToRemove = driver(newSummaryFile + toString(getpid()) + ".temp", minTms, maxTms, primers, conSeqs, lines[process].start, lines[process].end, numBinsProcessed);
                
                string tempFile = toString(getpid()) + ".otus2Remove.temp";
                ofstream outTemp;
                m->openOutputFile(tempFile, outTemp);
                
                outTemp << numBinsProcessed << endl;
                outTemp << otusToRemove.size() << endl;
                for (set<int>::iterator it = otusToRemove.begin(); it != otusToRemove.end(); it++) { outTemp << *it << endl; }
                outTemp.close();
                
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		otusToRemove = driver(newSummaryFile, minTms, maxTms, primers, conSeqs, lines[0].start, lines[0].end, numBinsProcessed);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            string tempFile = toString(processIDS[i]) +  ".otus2Remove.temp";
            ifstream intemp;
            m->openInputFile(tempFile, intemp);
            
            int num;
            intemp >> num; m->gobble(intemp);
            if (num != (lines[i+1].end - lines[i+1].start)) { m->mothurOut("[ERROR]: process " + toString(processIDS[i]) + " did not complete processing all OTUs assigned to it, quitting.\n"); m->control_pressed = true; }
            intemp >> num; m->gobble(intemp);
            for (int k = 0; k < num; k++) {
                int otu;
                intemp >> otu; m->gobble(intemp);
                otusToRemove.insert(otu); 
            }
            intemp.close();
            m->mothurRemove(tempFile);
        }

        
    #else
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the primerDesignData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<primerDesignData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
			// Allocate memory for thread data.
			string extension = toString(i) + ".temp";
			m->mothurRemove(newSummaryFile+extension);
            
			primerDesignData* tempPrimer = new primerDesignData((newSummaryFile+extension), m, lines[i].start, lines[i].end, minTms, maxTms, primers, conSeqs, pdiffs, otunumber, length, i);
			pDataArray.push_back(tempPrimer);
			processIDS.push_back(i);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyPrimerThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
        
		//using the main process as a worker saves time and memory
		otusToRemove = driver(newSummaryFile, minTms, maxTms, primers, conSeqs, lines[0].start, lines[0].end, numBinsProcessed);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            for (set<int>::iterator it = pDataArray[i]->otusToRemove.begin(); it != pDataArray[i]->otusToRemove.end(); it++) { 
				otusToRemove.insert(*it);  
			}
            int num = pDataArray[i]->numBinsProcessed;
            if (num != (lines[processIDS[i]].end - lines[processIDS[i]].start)) { m->mothurOut("[ERROR]: process " + toString(processIDS[i]) + " did not complete processing all OTUs assigned to it, quitting.\n"); m->control_pressed = true; }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
#endif		
		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((newSummaryFile + toString(processIDS[i]) + ".temp"), newSummaryFile);
			m->mothurRemove((newSummaryFile + toString(processIDS[i]) + ".temp"));
		}
		
		return otusToRemove;	
		
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************
set<int> PrimerDesignCommand::driver(string summaryFileName, vector<double>& minTms, vector<double>& maxTms, set<string>& primers, vector<Sequence>& conSeqs, int start, int end, int& numBinsProcessed){
	try {
        set<int> otuToRemove;
        
        ofstream outSum;
        m->openOutputFileAppend(summaryFileName, outSum);
        
        for (int i = start; i < end; i++) {
        
            if (m->control_pressed) { break; }
            
            if (i != (otunumber-1)) {
                int primerIndex = 0;
                for (set<string>::iterator it = primers.begin(); it != primers.end(); it++) {
                    vector<int> primerStarts;
                    vector<int> primerEnds;
                    vector<int> mismatches;
                    
                    bool found = findPrimer(conSeqs[i].getUnaligned(), (*it), primerStarts, primerEnds, mismatches);
                    
                    //if we found it report to the table
                    if (found) {
                        for (int j = 0; j < primerStarts.size(); j++) {
                            outSum << (i+1) << '\t' << *it << '\t' << primerStarts[j] << '\t' << primerEnds[j] << '\t' << length << '\t' << mismatches[j] << '\t' << minTms[primerIndex] << '\t' << maxTms[primerIndex] << endl;
                        }
                        otuToRemove.insert(i);
                    }
                    primerIndex++;
                }
            }
            numBinsProcessed++;
        }
        outSum.close();
        
        
        return otuToRemove;
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/ 
vector< vector< vector<unsigned int> > > PrimerDesignCommand::driverGetCounts(map<string, int>& nameMap, unsigned long int& fastaCount, vector<unsigned int>& otuCounts, unsigned long long& start, unsigned long long& end){
    try {
        vector< vector< vector<unsigned int> > > counts;
        map<string, int> seq2Bin;
        alignedLength = 0;
        
        ifstream in;
		m->openInputFile(fastafile, in);
        
		in.seekg(start);
        
		bool done = false;
		fastaCount = 0;
        
		while (!done) {
            if (m->control_pressed) { in.close(); return counts; }
            
			Sequence seq(in); m->gobble(in);
            
			if (seq.getName() != "") {
                if (fastaCount == 0) { alignedLength = seq.getAligned().length(); initializeCounts(counts, alignedLength, seq2Bin, nameMap, otuCounts); }
                else if (alignedLength != seq.getAligned().length()) {
                    m->mothurOut("[ERROR]: your sequences are not all the same length. primer.design requires sequences to be aligned."); m->mothurOutEndLine(); m->control_pressed = true; break;
                }
                
                int num = 1;
                map<string, int>::iterator itCount;
                if (namefile != "") { 
                    itCount = nameMap.find(seq.getName());
                    if (itCount == nameMap.end()) {  m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your name file, aborting."); m->mothurOutEndLine(); m->control_pressed = true; break; }
                    else { num = itCount->second; }
                    fastaCount+=num;
                }else if (countfile != "") {
                    itCount = nameMap.find(seq.getName());
                    if (itCount == nameMap.end()) {  m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your count file, aborting."); m->mothurOutEndLine(); m->control_pressed = true; break; }
                    else { num = itCount->second; }
                    fastaCount++;
                }else {
                    fastaCount++;
                }
                
                //increment counts
                itCount = seq2Bin.find(seq.getName());
                if (itCount == seq2Bin.end()) {
                    if ((namefile != "") || (countfile != "")) {
                        m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your list file, aborting. Perhaps you forgot to include your name or count file while clustering.\n"); m->mothurOutEndLine(); m->control_pressed = true; break;
                    }else{
                        m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your list file, aborting."); m->mothurOutEndLine(); m->control_pressed = true; break;
                    }
                }else {
                    otuCounts[itCount->second] += num;
                    string aligned = seq.getAligned();
                    for (int i = 0; i < alignedLength; i++) {
                        char base = toupper(aligned[i]);
                        if (base == 'A') { counts[itCount->second][i][0]+=num; }
                        else if (base == 'T') { counts[itCount->second][i][1]+=num; }
                        else if (base == 'G') { counts[itCount->second][i][2]+=num; }
                        else if (base == 'C') { counts[itCount->second][i][3]+=num; }
                        else { counts[itCount->second][i][4]+=num; }
                    }
                }

            }

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= end)) { break; }
#else
            if (in.eof()) { break; }
#endif
		}
				
		in.close();
        
        return counts;
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "driverGetCounts");
		exit(1);
	}
}
/**************************************************************************************************/
vector<Sequence> PrimerDesignCommand::createProcessesConSeqs(map<string, int>& nameMap, unsigned long int& numSeqs) {
	try {
        vector< vector< vector<unsigned int> > > counts;
        vector<unsigned int> otuCounts;
		vector<int> processIDS;
		int process = 1;
        unsigned long int fastaCount = 0;
 		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		
        vector<unsigned long long> positions; 
        vector<fastaLinePair> lines;
        positions = m->divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(fastaLinePair(positions[i], positions[(i+1)]));	}

		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
 				counts = driverGetCounts(nameMap, fastaCount, otuCounts, lines[process].start, lines[process].end);
                
                string tempFile = toString(getpid()) + ".cons_counts.temp";
                ofstream outTemp;
                m->openOutputFile(tempFile, outTemp);
                
                outTemp << fastaCount << endl;
                //pass counts
                outTemp << counts.size() << endl;
                for (int i = 0; i < counts.size(); i++) { 
                    outTemp << counts[i].size() << endl;
                    for (int j = 0; j < counts[i].size(); j++) { 
                        for (int k = 0; k < 5; k++) {  outTemp << counts[i][j][k] << '\t'; }
                        outTemp << endl;
                    }
                }
                //pass otuCounts
                outTemp << otuCounts.size() << endl;
                for (int i = 0; i < otuCounts.size(); i++) { outTemp << otuCounts[i] << '\t'; }
                outTemp << endl;
                outTemp.close();
                
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		counts = driverGetCounts(nameMap, fastaCount, otuCounts, lines[0].start, lines[0].end);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            string tempFile = toString(processIDS[i]) +  ".cons_counts.temp";
            ifstream intemp;
            m->openInputFile(tempFile, intemp);
            
            unsigned long int num;
            intemp >> num; m->gobble(intemp); fastaCount += num;
            intemp >> num; m->gobble(intemp);
            if (num != counts.size()) { m->mothurOut("[ERROR]: " + tempFile + " was not built correctly by the child process, quitting.\n"); m->control_pressed = true; }
            else {
                //read counts
                for (int k = 0; k < num; k++) {
                    int alength;
                    intemp >> alength; m->gobble(intemp);
                    if (alength != alignedLength) {  m->mothurOut("[ERROR]: your sequences are not all the same length. primer.design requires sequences to be aligned."); m->mothurOutEndLine(); m->control_pressed = true; }
                    else {
                        for (int j = 0; j < alength; j++) {
                            for (int l = 0; l < 5; l++) {  unsigned int numTemp; intemp >> numTemp; m->gobble(intemp); counts[k][j][l] += numTemp;  }
                        }
                    }
                }
                //read otuCounts
                intemp >> num; m->gobble(intemp);
                for (int k = 0; k < num; k++) {  
                    unsigned int numTemp; intemp >> numTemp; m->gobble(intemp); 
                    otuCounts[k] += numTemp; 
                }
            }
            intemp.close();
            m->mothurRemove(tempFile);
        }
        
        
#else
        unsigned long long start = 0;
        unsigned long long end = 1000;
		counts = driverGetCounts(nameMap, fastaCount, otuCounts, start, end);	
#endif		
        
        //you will have a nameMap error if there is a namefile or countfile, but if those aren't given we want to make sure the fasta and list file match.
        if (fastaCount != numSeqs) {
            if ((namefile == "") && (countfile == ""))   {  m->mothurOut("[ERROR]: Your list file contains " + toString(list->getNumSeqs()) + " sequences, and your fasta file contains " + toString(fastaCount) + " sequences, aborting. Do you have the correct files? Perhaps you forgot to include the name or count file? \n");   }
            m->control_pressed = true;
        }
        
		vector<Sequence> conSeqs;
        
        if (m->control_pressed) { return conSeqs; }
        
		//build consensus seqs
        string snumBins = toString(counts.size());
        for (int i = 0; i < counts.size(); i++) {
            if (m->control_pressed) { break; }
            
            string otuLabel = "Otu";
            string sbinNumber = toString(i+1);
            if (sbinNumber.length() < snumBins.length()) { 
                int diff = snumBins.length() - sbinNumber.length();
                for (int h = 0; h < diff; h++) { otuLabel += "0"; }
            }
            otuLabel += sbinNumber; 
            
            string cons = "";
            for (int j = 0; j < counts[i].size(); j++) {
                cons += getBase(counts[i][j], otuCounts[i]);
            }
            Sequence consSeq(otuLabel, cons);
            conSeqs.push_back(consSeq);
        }
        
        if (m->control_pressed) { conSeqs.clear(); return conSeqs; }
        
        return conSeqs;
	
		
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "createProcessesConSeqs");
		exit(1);
	}
}
//***************************************************************************************************************

char PrimerDesignCommand::getBase(vector<unsigned int> counts, int size){  //A,T,G,C,Gap
	try{
		/* A = adenine
         * C = cytosine
         * G = guanine
         * T = thymine
         * R = G A (purine)
         * Y = T C (pyrimidine)
         * K = G T (keto)
         * M = A C (amino)
         * S = G C (strong bonds)
         * W = A T (weak bonds)
         * B = G T C (all but A)
         * D = G A T (all but C)
         * H = A C T (all but G)
         * V = G C A (all but T)
         * N = A G C T (any) */
		
		char conBase = 'N';
		
		//zero out counts that don't make the cutoff
		float percentage = (100.0 - cutoff) / 100.0;
        
		for (int i = 0; i < counts.size(); i++) {
            float countPercentage = counts[i] / (float) size;
			if (countPercentage < percentage) { counts[i] = 0; }
		}
		
		//any
		if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'n'; }
		//any no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'N'; }
		//all but T
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'v'; }	
		//all but T no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'V'; }	
		//all but G
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'h'; }	
		//all but G no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'H'; }	
		//all but C
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'd'; }	
		//all but C no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'D'; }	
		//all but A
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'b'; }	
		//all but A no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'B'; }	
		//W = A T (weak bonds)
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'w'; }	
		//W = A T (weak bonds) no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'W'; }	
		//S = G C (strong bonds)
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 's'; }	
		//S = G C (strong bonds) no gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'S'; }	
		//M = A C (amino)
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'm'; }	
		//M = A C (amino) no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'M'; }	
		//K = G T (keto)
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'k'; }	
		//K = G T (keto) no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'K'; }	
		//Y = T C (pyrimidine)
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'y'; }	
		//Y = T C (pyrimidine) no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'Y'; }	
		//R = G A (purine)
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'r'; }	
		//R = G A (purine) no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'R'; }	
		//only A
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'a'; }	
		//only A no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'A'; }	
		//only T
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 't'; }	
		//only T no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'T'; }	
		//only G
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'g'; }	
		//only G no gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'G'; }	
		//only C
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'c'; }	
		//only C no gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'C'; }	
		//only gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = '-'; }
		//cutoff removed all counts
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'N'; }
		else{ m->mothurOut("[ERROR]: cannot find consensus base."); m->mothurOutEndLine(); }
		
		return conBase;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "getBase");
		exit(1);
	}
}

//**********************************************************************************************************************
int PrimerDesignCommand::initializeCounts(vector< vector< vector<unsigned int> > >& counts, int length, map<string, int>& seq2Bin, map<string, int>& nameMap, vector<unsigned int>& otuCounts){
	try {
        counts.clear();
        otuCounts.clear();
        seq2Bin.clear();
        
        //vector< vector< vector<unsigned int> > > counts - otu < spot_in_alignment < counts_for_A,T,G,C,Gap > > >
        for (int i = 0; i < list->getNumBins(); i++) {
            string binNames = list->get(i);
            vector<string> names;
            m->splitAtComma(binNames, names);
            otuCounts.push_back(0);
            
            //lets be smart and only map the unique names if a name or count file was given to save search time and memory
            if ((namefile != "") || (countfile != "")) {
                for (int j = 0; j < names.size(); j++) {
                    map<string, int>::iterator itNames = nameMap.find(names[j]);
                    if (itNames != nameMap.end()) { //add name because its a unique one
                        seq2Bin[names[j]] = i;
                    }
                }
            }else { //map everyone
                for (int j = 0; j < names.size(); j++) { seq2Bin[names[j]] = i;  }
            }
            
            vector<unsigned int> temp; temp.resize(5, 0); //A,T,G,C,Gap
            vector< vector<unsigned int> > temp2;
            for (int j = 0; j < length; j++) {
                temp2.push_back(temp);
            }
            counts.push_back(temp2);
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "initializeCounts");
		exit(1);
	}
}
//**********************************************************************************************************************
map<string, int> PrimerDesignCommand::readCount(unsigned long int& numSeqs){
	try {
        map<string, int> nameMap;
        
        CountTable ct;
        ct.readTable(countfile);
        vector<string> namesOfSeqs = ct.getNamesOfSeqs();
        numSeqs = ct.getNumUniqueSeqs();
        
        for (int i = 0; i < namesOfSeqs.size(); i++) {
            if (m->control_pressed) { break; }
            
            nameMap[namesOfSeqs[i]] = ct.getNumSeqs(namesOfSeqs[i]);
        }
        
        return nameMap;
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "readCount");
		exit(1);
	}
}
//**********************************************************************************************************************
int PrimerDesignCommand::getListVector(){
	try {
		InputData input(listfile, "list");
		list = input.getListVector();
		string lastLabel = list->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			if (m->control_pressed) {  return 0;  }
			
			if(labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input.getListVector(lastLabel);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
				break;
			}
			
			lastLabel = list->getLabel();			
			
			//get next line to process
			//prevent memory leak
			delete list; 
			list = input.getListVector();
		}
		
		
		if (m->control_pressed) {  return 0;  }
		
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
			delete list; 
			list = input.getListVector(lastLabel);
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "getListVector");	
		exit(1);
	}
}
//********************************************************************/
/* A = adenine
 * C = cytosine
 * G = guanine
 * T = thymine
 * R = G A (purine)
 * Y = T C (pyrimidine)
 * K = G T (keto)
 * M = A C (amino)
 * S = G C (strong bonds)
 * W = A T (weak bonds)
 * B = G T C (all but A)
 * D = G A T (all but C)
 * H = A C T (all but G)
 * V = G C A (all but T)
 * N = A G C T (any) */
int PrimerDesignCommand::countDiffs(string oligo, string seq){
	try {
		
		int length = oligo.length();
		int countDiffs = 0;
		
		for(int i=0;i<length;i++){
            
			oligo[i] = toupper(oligo[i]);
            seq[i] = toupper(seq[i]);
            
			if(oligo[i] != seq[i]){
                if(oligo[i] == 'A' && (seq[i] != 'A' && seq[i] != 'M' && seq[i] != 'R' && seq[i] != 'W' && seq[i] != 'D' && seq[i] != 'H' && seq[i] != 'V'))       {	countDiffs++;	}
                else if(oligo[i] == 'C' && (seq[i] != 'C' && seq[i] != 'Y' && seq[i] != 'M' && seq[i] != 'S' && seq[i] != 'B' && seq[i] != 'H' && seq[i] != 'V'))       {	countDiffs++;	}
                else if(oligo[i] == 'G' && (seq[i] != 'G' && seq[i] != 'R' && seq[i] != 'K' && seq[i] != 'S' && seq[i] != 'B' && seq[i] != 'D' && seq[i] != 'V'))       {	countDiffs++;	}
                else if(oligo[i] == 'T' && (seq[i] != 'T' && seq[i] != 'Y' && seq[i] != 'K' && seq[i] != 'W' && seq[i] != 'B' && seq[i] != 'D' && seq[i] != 'H'))       {	countDiffs++;	}
                else if((oligo[i] == '.' || oligo[i] == '-'))           {	countDiffs++;	}
				else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))                         {	countDiffs++;	}
				else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))                        {	countDiffs++;	}
				else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))                        {	countDiffs++;	}
				else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))                        {	countDiffs++;	}
				else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))                        {	countDiffs++;	}
				else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))                        {	countDiffs++;	}
				else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))                        {	countDiffs++;	}
				else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))       {	countDiffs++;	}
				else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))       {	countDiffs++;	}
				else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))       {	countDiffs++;	}
				else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))       {	countDiffs++;	}	
            }
			
		}
		
		return countDiffs;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "countDiffs");
		exit(1);
	}
}
//**********************************************************************************************************************



