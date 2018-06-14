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
        CommandParameter potunumber("otulabel", "String", "", "", "", "", "","",false,true,true); parameters.push_back(potunumber);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
        CommandParameter pcutoff("cutoff", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
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
		helpString += "The primer.design command parameters are: list, fasta, name, count, otulabel, cutoff, length, pdiffs, mintm, maxtm, processors and label.\n";
		helpString += "The list parameter allows you to provide a list file and is required.\n";
        helpString += "The fasta parameter allows you to provide a fasta file and is required.\n";
        helpString += "The name parameter allows you to provide a name file associated with your fasta file.\n";
        helpString += "The count parameter allows you to provide a count file associated with your fasta file.\n";
        helpString += "The label parameter is used to indicate the label you want to use from your list file.\n";
        helpString += "The otulabel parameter is used to indicate the otu you want to use from your list file. It is required.\n";
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
 				string path;
				it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
                
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
            }
                        
			//check for parameters
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            //get fastafile - it is required
            fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort=true;  }
			else if (fastafile == "not found") {  
                fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
            }else  { current->setFastaFile(fastafile); }
            
            //get listfile - it is required
            listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { listfile = ""; abort=true;  }
			else if (listfile == "not found") {  
                listfile = current->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current listfile and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
            }else  { current->setListFile(listfile); }

            
			if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = util.hasPath(listfile); //if user entered a file with a path then preserve it	
			}
            
            string temp = validParameter.valid(parameters, "cutoff");  if (temp == "not found") { temp = "100"; }
			util.mothurConvert(temp, cutoff); 
            
            temp = validParameter.valid(parameters, "pdiffs");  if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, pdiffs); 
            
            temp = validParameter.valid(parameters, "length");  if (temp == "not found") { temp = "18"; }
			util.mothurConvert(temp, length); 
            
            temp = validParameter.valid(parameters, "mintm");  if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, minTM); 
            
            temp = validParameter.valid(parameters, "maxtm");  if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxTM); 
            
            otulabel = validParameter.valid(parameters, "otulabel");  if (otulabel == "not found") { otulabel = ""; }
            if (otulabel == "") {  m->mothurOut("[ERROR]: You must provide an OTU label, aborting.\n"); abort = true; }
            
            temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
            
            label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }
        
            if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        long start = time(NULL);
        //////////////////////////////////////////////////////////////////////////////
        //              get file inputs                                             //
        //////////////////////////////////////////////////////////////////////////////
        
        //reads list file and selects the label the users specified or the first label
        ListVector* list = getListVector();
        vector<string> binLabels = list->getLabels();
        int binIndex = findIndex(otulabel, binLabels);
        if (binIndex == -1) { m->mothurOut("[ERROR]: You selected an OTU label that is not in your in your list file, quitting.\n"); return 0; }
        
        map<string, int> nameMap;
        long long numSeqs;  //used to sanity check the files. numSeqs = total seqs for namefile and uniques for count.
                                    //list file should have all seqs if namefile was used to create it and only uniques in count file was used.
        
        if (namefile != "")         {  unsigned long int temp; nameMap = util.readNames(namefile, temp);   numSeqs = temp;      }
        else if (countfile != "")   {
            CountTable ct; ct.readTable(countfile, false, false);
            numSeqs = ct.getNumUniqueSeqs();
            nameMap = ct.getNameMap();
        }else { numSeqs = list->getNumSeqs();  }
        
        //sanity check
        if (numSeqs != list->getNumSeqs()) {
            if (namefile != "")         {  m->mothurOut("[ERROR]: Your list file contains " + toString(list->getNumSeqs()) + " sequences, and your name file contains " + toString(numSeqs) + " sequences, aborting. Do you have the correct files? Perhaps you forgot to include the name file when you clustered? \n");   }
            else if (countfile != "") {
                m->mothurOut("[ERROR]: Your list file contains " + toString(list->getNumSeqs()) + " sequences, and your count file contains " + toString(numSeqs) + " unique sequences, aborting. Do you have the correct files? Perhaps you forgot to include the count file when you clustered? \n");  
            }
            m->setControl_pressed(true);
        }
        
        if (m->getControl_pressed()) { delete list; return 0; }
        
        //////////////////////////////////////////////////////////////////////////////
        //              process data                                                //
        //////////////////////////////////////////////////////////////////////////////
        m->mothurOut("\nFinding consensus sequences for each otu...\n");
        
        map<string, int> seq2Bin = getSequenceBinAssignments(list, nameMap);
        
        vector<Sequence> conSeqs = createProcessesConSeqs(nameMap, seq2Bin, binLabels);
        
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[distance]"] = list->getLabel();
        string consFastaFile = getOutputFileName("fasta", variables);
        outputNames.push_back(consFastaFile); outputTypes["fasta"].push_back(consFastaFile);
        ofstream out;
        util.openOutputFile(consFastaFile, out);
        for (int i = 0; i < conSeqs.size(); i++) {  conSeqs[i].printSequence(out);  }
        out.close();
        
        set<string> primers = getPrimer(conSeqs[binIndex]);
        
        if (m->getControl_pressed()) { delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        
        string consSummaryFile = getOutputFileName("summary", variables);
        outputNames.push_back(consSummaryFile); outputTypes["summary"].push_back(consSummaryFile);
        ofstream outSum;
        util.openOutputFile(consSummaryFile, outSum);
        
        outSum << "PrimerOtu: " << otulabel << " Members: " << list->get(binIndex) << endl << "Primers\tminTm\tmaxTm" << endl;
        
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
        
        outSum << "\nOTU\tPrimer\tStart\tEnd\tLength\tMismatches\tminTm\tmaxTm\n";
        outSum.close();
        
        m->mothurOut("\nProcessing OTUs...\n");
        
        //check each otu's conseq for each primer in otunumber
        set<int> otuToRemove = createProcesses(consSummaryFile, minTms, maxTms, primers, conSeqs, binIndex, binLabels);
        
        if (m->getControl_pressed()) { delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        
        //print new list file
        map<string, string> mvariables; 
        mvariables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile));
        mvariables["[extension]"] = util.getExtension(listfile);
        string newListFile = getOutputFileName("list", mvariables);
        ofstream outListTemp;
        util.openOutputFile(newListFile+".temp", outListTemp);
        
        outListTemp << list->getLabel() << '\t' << (list->getNumBins()-otuToRemove.size());
        string headers = "label\tnumOtus";
        for (int j = 0; j < list->getNumBins(); j++) {
            if (m->getControl_pressed()) { break; }
            //good otus
            if (otuToRemove.count(j) == 0) {  
                string bin = list->get(j);
                if (bin != "") {  outListTemp << '\t' << bin;  headers += '\t' + binLabels[j]; }
            }
        }
        outListTemp << endl;
        outListTemp.close();
        
        ofstream outList;
        util.openOutputFile(newListFile, outList);
        outList << headers << endl;
        outList.close();
        util.appendFiles(newListFile+".temp", newListFile);
        util.mothurRemove(newListFile+".temp");
        outputNames.push_back(newListFile); outputTypes["list"].push_back(newListFile);
        delete list;
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process " + toString(list->getNumBins()) + " OTUs.\n");
        
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        
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
//find all primers for the given sequence
set<string> PrimerDesignCommand::getPrimer(Sequence primerSeq){
	try {
        set<string> primers;
        
        string rawSequence = primerSeq.getUnaligned();
        
        for (int j = 0; j < rawSequence.length()-length; j++){
            if (m->getControl_pressed()) { break; }
            
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
int countDiffs(string oligo, string seq, MothurOut* m){
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
//********************************************************************/
//search for a primer over the sequence string
bool findPrimer(string rawSequence, string primer, vector<int>& primerStart, vector<int>& primerEnd, vector<int>& mismatches, int length, int pdiffs, MothurOut* m){
    try {
        bool foundAtLeastOne = false;  //innocent til proven guilty
        
        //look for exact match
        if(rawSequence.length() < primer.length()) {  return false;  }
        
        //search for primer
        for (int j = 0; j < rawSequence.length()-length; j++){
            
            if (m->getControl_pressed()) {  return foundAtLeastOne; }
            
            string rawChunk = rawSequence.substr(j, length);
            
            int numDiff = countDiffs(primer, rawChunk, m);
            
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

/**************************************************************************************************/
struct primerDesignData {
    OutputWriter* summaryFile;
    MothurOut* m;
    int start;
    int end;
    int pdiffs, length, binIndex;
    set<string> primers;
    vector<double> minTms, maxTms;
    set<int> otusToRemove;
    vector<Sequence> consSeqs;
    vector<string> binLabels;
    int numBinsProcessed;
    Utils util;
    
    primerDesignData(){ delete summaryFile; }
    primerDesignData(OutputWriter* sf, int st, int en, vector<double> min, vector<double> max, set<string> pri, vector<Sequence> seqs, int d, int otun, int l, vector<string> bl) {
        summaryFile = sf;
        m = MothurOut::getInstance();
        start = st;
        end = en;
        pdiffs = d;
        minTms = min;
        maxTms = max;
        primers = pri;
        consSeqs = seqs;
        binIndex = otun;
        length = l;
        binLabels = bl;
        numBinsProcessed = 0;
    }
};
//**********************************************************************************************************************
void driverPDesign(primerDesignData* params){
    try {
        
        for (int i = params->start; i < params->end; i++) {
            
            if (params->m->getControl_pressed()) { break; }
            
            if (i != (params->binIndex)) {
                int primerIndex = 0;
                string output = "";
                for (set<string>::iterator it = params->primers.begin(); it != params->primers.end(); it++) {
                    vector<int> primerStarts;
                    vector<int> primerEnds;
                    vector<int> mismatches;
                    
                    bool found = findPrimer(params->consSeqs[i].getUnaligned(), (*it), primerStarts, primerEnds, mismatches, params->length, params->pdiffs, params->m);
                    
                    //if we found it report to the table
                    if (found) {
                        for (int j = 0; j < primerStarts.size(); j++) {
                            output += params->binLabels[i] + '\t' + *it + '\t' + toString(primerStarts[j]) + '\t' + toString(primerEnds[j]) + '\t' + toString(params->length) + '\t' + toString(mismatches[j]) + '\t' + toString(params->minTms[primerIndex]) + '\t' + toString(params->maxTms[primerIndex]) + '\n';
                        }
                        params->otusToRemove.insert(i);
                    }
                    primerIndex++;
                }
                params->summaryFile->write(output);
            }
            params->numBinsProcessed++;
            
            if((params->numBinsProcessed) % 100 == 0){	params->m->mothurOutJustToScreen(toString(params->numBinsProcessed)+"\n"); 		}
        }
        
        if((params->numBinsProcessed) % 100 != 0){	params->m->mothurOutJustToScreen(toString(params->numBinsProcessed)+"\n"); 		}
     }
    catch(exception& e) {
        params->m->errorOut(e, "PrimerDesignCommand", "driver");
        exit(1);
    }
}
/**************************************************************************************************/
set<int> PrimerDesignCommand::createProcesses(string newSummaryFile, vector<double>& minTms, vector<double>& maxTms, set<string>& primers, vector<Sequence>& conSeqs, int binIndex, vector<string>& binLabels) {
	try {
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<primerDesignData*> data;
		
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
		
        auto synchronizedOutputSummaryFile = std::make_shared<SynchronizedOutputFile>(newSummaryFile, true); //open append
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadSummaryWriter = new OutputWriter(synchronizedOutputSummaryFile);
            
            primerDesignData* dataBundle = new primerDesignData(threadSummaryWriter, lines[i+1].start, lines[i+1].end, minTms, maxTms, primers, conSeqs, pdiffs, binIndex, length, binLabels);
            
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverPDesign, dataBundle));
        }
        
        OutputWriter* threadSummaryWriter = new OutputWriter(synchronizedOutputSummaryFile);
        primerDesignData* dataBundle = new primerDesignData(threadSummaryWriter, lines[0].start, lines[0].end, minTms, maxTms, primers, conSeqs, pdiffs, binIndex, length, binLabels);
        
        driverPDesign(dataBundle);
        set<int> otusToRemove = dataBundle->otusToRemove;
    
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
    
            otusToRemove.insert(data[i]->otusToRemove.begin(), data[i]->otusToRemove.end());
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
		return otusToRemove;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************
int initializeCounts(vector< vector< vector<unsigned int> > >& counts, int length, int numBins, MothurOut* m){
    try {
        counts.clear();
        
        //vector< vector< vector<unsigned int> > > counts - otu < spot_in_alignment < counts_for_A,T,G,C,Gap > > >
        for (int i = 0; i < numBins; i++) {
            
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
/**************************************************************************************************/
struct primerCountsData {
    map<string, int> seq2Bin;
    vector< vector< vector<unsigned int> > > counts;  // - otu < spot_in_alignment < counts_for_A,T,G,C,Gap > > >
    MothurOut* m;
    long long start, end, count, total, alignedLength, numBins;
    string fastafile;
    Utils util;
    map<string, int> nameMap;
    vector<long long> otuCounts;
    bool hasNameMap;
    
    primerCountsData(){  }
    primerCountsData(string ff, map<string, int> nmp, long long st, long long en, map<string, int> seq2B, int nb) {
        fastafile = ff;
        m = MothurOut::getInstance();
        start = st;
        end = en;
        hasNameMap = false;
        nameMap = nmp;
        if (nameMap.size() != 0) { hasNameMap = true; }
        seq2Bin = seq2B;
        numBins = nb;
        count = 0; total = 0;
    }
};
//**********************************************************************************************************************
map<string, int> PrimerDesignCommand::getSequenceBinAssignments(ListVector* list, map<string, int>& nameMap){
    try {
        map<string, int> seq2Bin;
        bool hasNameMap = false;
        if (nameMap.size() != 0) { hasNameMap = true; }
        
        for (int i = 0; i < list->getNumBins(); i++) {
            string binNames = list->get(i);
            vector<string> names;
            util.splitAtComma(binNames, names);
            
            //lets be smart and only map the unique names if a name or count file was given to save search time and memory
            if (hasNameMap) {
                for (int j = 0; j < names.size(); j++) {
                    map<string, int>::iterator itNames = nameMap.find(names[j]);
                    if (itNames != nameMap.end()) { //add name because its a unique one
                        seq2Bin[names[j]] = i;
                    }
                }
            }else { for (int j = 0; j < names.size(); j++) { seq2Bin[names[j]] = i;  } } //map everyone
        }

        return seq2Bin;
    }
    catch(exception& e) {
        m->errorOut(e, "PrimerDesignCommand", "getSequenceBinAssignments");
        exit(1);
    }
}
/**************************************************************************************************/
void driverGetCounts(primerCountsData* params){
    try {
        params->otuCounts.resize(params->numBins, 0);
        params->alignedLength = 0;
        
        ifstream in; params->util.openInputFile(params->fastafile, in);
		in.seekg(params->start);
        
        //adjust start if null strings
        if (params->start == 0) {  params->util.zapGremlins(in); params->util.gobble(in);  }
        
		bool done = false;
    
		while (!done) {
            if (params->m->getControl_pressed()) { break; }
            
			Sequence seq(in); params->util.gobble(in);
            
			if (seq.getName() != "") {
                if (params->count == 0) { params->alignedLength = seq.getAligned().length(); initializeCounts(params->counts, params->alignedLength, params->numBins, params->m); }
                else if (params->alignedLength != seq.getAligned().length()) {
                    params->m->mothurOut("[ERROR]: your sequences are not all the same length. primer.design requires sequences to be aligned.\n"); params->m->setControl_pressed(true); break;
                }
                
                int num = 1;
                if (params->hasNameMap) {
                    map<string, int>::iterator itCount = params->nameMap.find(seq.getName());
                    if (itCount == params->nameMap.end()) {  params->m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your name or count file, aborting.\n");  params->m->setControl_pressed(true); break; }
                    else { params->total += itCount->second; num = itCount->second; }
                }else { params->total++; }
                
                //increment counts
                map<string, int>::iterator itCount = params->seq2Bin.find(seq.getName());
                if (itCount == params->seq2Bin.end()) { params->m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your list file, aborting.\n");  params->m->setControl_pressed(true); break;
                }else {
                    params->otuCounts[itCount->second] += num;
                    string aligned = seq.getAligned();
                    for (int i = 0; i < params->alignedLength; i++) {
                        char base = toupper(aligned[i]);
                        if (base == 'A') { params->counts[itCount->second][i][0]+=num; }
                        else if (base == 'T') { params->counts[itCount->second][i][1]+=num; }
                        else if (base == 'G') { params->counts[itCount->second][i][2]+=num; }
                        else if (base == 'C') { params->counts[itCount->second][i][3]+=num; }
                        else { params->counts[itCount->second][i][4]+=num; }
                    }
                }

            }
            params->count++;
            
            if((params->count) % 1000 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 		}
            
#if defined NON_WINDOWS
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= params->end)) { break; }
#else
            if ((params->count == params->end) || (in.eof())) { break; }
#endif
		}
		
        if((params->count) % 1000 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n");		}
        
		in.close();
    }
	catch(exception& e) {
		params->m->errorOut(e, "PrimerDesignCommand", "driverGetCounts");
		exit(1);
	}
}
/**************************************************************************************************/
vector<Sequence> PrimerDesignCommand::createProcessesConSeqs(map<string, int>& nameMap, map<string, int>& seq2Bin, vector<string>& binLabels) {
	try {
        int numBins = binLabels.size();
        vector<linePair> lines;
#if defined NON_WINDOWS
        vector<unsigned long long> positions;  positions = util.divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        long long numSeqs;
        vector<unsigned long long> positions = util.setFilePosFasta(fastafile, numSeqs);
        if (numSeqs < processors) { processors = numSeqs; m->mothurOut("Reducing processors to " + toString(numSeqs) + ".\n");  }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<primerCountsData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            if (m->getControl_pressed()) {  break; }
            
            primerCountsData* dataBundle = new primerCountsData(fastafile, nameMap, lines[i+1].start, lines[i+1].end, seq2Bin, numBins);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverGetCounts, dataBundle));
        }
        
        primerCountsData* dataBundle = new primerCountsData(fastafile, nameMap, lines[0].start, lines[0].end, seq2Bin, numBins);
        driverGetCounts(dataBundle);
        vector< vector< vector<unsigned int> > > counts = dataBundle->counts;
        vector<long long> otuCounts = dataBundle->otuCounts;
        long long total = dataBundle->total;
        int alignedLength = dataBundle->alignedLength;
        delete dataBundle;
     
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            total += data[i]->total;
            
            if (m->getControl_pressed()) {  break; }
            
            if (data[i]->alignedLength != alignedLength) {  m->mothurOut("[ERROR]: your sequences are not all the same length. primer.design requires sequences to be aligned.\n"); m->setControl_pressed(true); }
            
            for (int k = 0; k < numBins; k++) { //for each bin
                for (int j = 0; j < data[i]->alignedLength; j++) { //for each position
                    for (int l = 0; l < 5; l++) {  counts[k][j][l] += data[i]->counts[k][j][l];  }  //for each base
                }
                otuCounts[k] += data[i]->otuCounts[k];
            }
            delete data[i];
            delete workerThreads[i];
        }
        vector<Sequence> conSeqs;
        if (m->getControl_pressed()) {  return conSeqs; }
        
		//build consensus seqs
        for (int i = 0; i < counts.size(); i++) {
            if (m->getControl_pressed()) {  break; }
            
            string otuLabel = binLabels[i];
            
            string cons = "";
            for (int j = 0; j < counts[i].size(); j++) { cons += getBase(counts[i][j], otuCounts[i]); }
            
            Sequence consSeq(otuLabel, cons);
            
            conSeqs.push_back(consSeq);
        }
        
        if (m->getControl_pressed()) { conSeqs.clear(); return conSeqs; }
        
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
ListVector* PrimerDesignCommand::getListVector(){
	try {
		InputData input(listfile, "list", nullVector);
		ListVector* list = input.getListVector();
		string lastLabel = list->getLabel();
		
		if (label == "") { label = lastLabel;  return list; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {  return list;  }
			
			if(labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
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
		
		
		if (m->getControl_pressed()) {  return list;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			delete list; 
			list = input.getListVector(lastLabel);
		}	
		
		return list;
	}
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "getListVector");	
		exit(1);
	}
}
//**********************************************************************************************************************
int PrimerDesignCommand::findIndex(string binLabel, vector<string> binLabels){
	try {
        int index = -1;
        for (int i = 0; i < binLabels.size(); i++){
            if (m->getControl_pressed()) { return index; }
            if (util.isLabelEquivalent(binLabel, binLabels[i])) { index = i; break; }
        }
        return index;
    }
	catch(exception& e) {
		m->errorOut(e, "PrimerDesignCommand", "findIndex");
		exit(1);
	}
}
//**********************************************************************************************************************



