//
//  prcseqscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "pcrseqscommand.h"

//**********************************************************************************************************************
vector<string> PcrSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
		CommandParameter poligos("oligos", "InputTypes", "", "", "ecolioligos", "none", "none","",false,false,true); parameters.push_back(poligos);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false,true); parameters.push_back(pgroup);
        CommandParameter ptax("taxonomy", "InputTypes", "", "", "none", "none", "none","taxonomy",false,false,true); parameters.push_back(ptax);
        CommandParameter pecoli("ecoli", "InputTypes", "", "", "ecolioligos", "none", "none","",false,false); parameters.push_back(pecoli);
		CommandParameter pstart("start", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pstart);
		CommandParameter pend("end", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pend);
 		CommandParameter pnomatch("nomatch", "Multiple", "reject-keep", "reject", "", "", "","",false,false); parameters.push_back(pnomatch);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
        CommandParameter prdiffs("rdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(prdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pkeepprimer("keepprimer", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkeepprimer);
        CommandParameter pkeepdots("keepdots", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pkeepdots);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PcrSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The pcr.seqs command reads a fasta file.\n";
        helpString += "The pcr.seqs command parameters are fasta, oligos, name, group, count, taxonomy, ecoli, start, end, nomatch, pdiffs, rdiffs, processors, keepprimer and keepdots.\n";
		helpString += "The ecoli parameter is used to provide a fasta file containing a single reference sequence (e.g. for e. coli) this must be aligned. Mothur will trim to the start and end positions of the reference sequence.\n";
        helpString += "The start parameter allows you to provide a starting position to trim to.\n";
        helpString += "The end parameter allows you to provide a ending position to trim from.\n";
        helpString += "The nomatch parameter allows you to decide what to do with sequences where the primer is not found. Default=reject, meaning remove from fasta file.  if nomatch=true, then do nothing to sequence.\n";
        helpString += "The processors parameter allows you to use multiple processors.\n";
        helpString += "The keepprimer parameter allows you to keep the primer, default=false.\n";
        helpString += "The keepdots parameter allows you to keep the leading and trailing .'s, default=true.\n";
        helpString += "The pdiffs parameter is used to specify the number of differences allowed in the forward primer. The default is 0.\n";
        helpString += "The rdiffs parameter is used to specify the number of differences allowed in the reverse primer. The default is 0.\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/Pcr.seqs .\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string PcrSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],pcr,[extension]-[filename],[tag],pcr,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "accnos")      {   pattern = "[filename],bad.accnos";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PcrSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

PcrSeqsCommand::PcrSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "PcrSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

PcrSeqsCommand::PcrSeqsCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
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
			outputTypes["fasta"] = tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
            outputTypes["accnos"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
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
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("ecoli");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ecoli"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
				
			}
            
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else if (fastafile == "not open") { fastafile = ""; abort = true; }	
			else { current->setFastaFile(fastafile); }
			
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(fastafile);	}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.valid(parameters, "keepprimer");  if (temp == "not found")    {	temp = "f";	}
			keepprimer = util.isTrue(temp);	
            
            temp = validParameter.valid(parameters, "keepdots");  if (temp == "not found")    {	temp = "t";	}
			keepdots = util.isTrue(temp);	
            
			temp = validParameter.validFile(parameters, "oligos");
			if (temp == "not found"){	oligosfile = "";		}
			else if(temp == "not open"){	oligosfile = ""; abort = true;	} 
			else					{	oligosfile = temp; current->setOligosFile(oligosfile);		}
			
            ecolifile = validParameter.validFile(parameters, "ecoli");
			if (ecolifile == "not found"){	ecolifile = "";		}
			else if(ecolifile == "not open"){	ecolifile = ""; abort = true;	} 
			
            namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found"){	namefile = "";		}
			else if(namefile == "not open"){	namefile = ""; abort = true;	} 
            else { current->setNameFile(namefile); }
            
            groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found"){	groupfile = "";		}
			else if(groupfile == "not open"){	groupfile = ""; abort = true;	} 
            else { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
            
            taxfile = validParameter.validFile(parameters, "taxonomy");
			if (taxfile == "not found"){	taxfile = "";		}
			else if(taxfile == "not open"){	taxfile = ""; abort = true;	} 
            else { current->setTaxonomyFile(taxfile); }
			 			
			temp = validParameter.valid(parameters, "start");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, start);
            
            temp = validParameter.valid(parameters, "end");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, end);
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
                    
            temp = validParameter.valid(parameters, "pdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, pdiffs);
            
            temp = validParameter.valid(parameters, "rdiffs");		if (temp == "not found") { temp = "0"; }
            util.mothurConvert(temp, rdiffs);
			
            nomatch = validParameter.valid(parameters, "nomatch");	if (nomatch == "not found") { nomatch = "reject"; }
			
            if ((nomatch != "reject") && (nomatch != "keep")) { m->mothurOut("[ERROR]: " + nomatch + " is not a valid entry for nomatch. Choices are reject and keep.\n");  abort = true; }
            
            //didnt set anything
			if ((oligosfile == "") && (ecolifile == "") && (start == -1) && (end == -1)) {
                m->mothurOut("[ERROR]: You did not set any options. Please provide an oligos or ecoli file, or set start or end.\n"); abort = true;
            }
            
            if ((oligosfile == "") && (ecolifile == "") && (start < 0) && (end == -1)) { m->mothurOut("[ERROR]: Invalid start value.\n"); abort = true; }
            
            if ((ecolifile != "") && (start != -1) && (end != -1)) {
                m->mothurOut("[ERROR]: You provided an ecoli file , but set the start or end parameters. Unsure what you intend.  When you provide the ecoli file, mothur thinks you want to use the start and end of the sequence in the ecoli file.\n"); abort = true;
            }

            
            if ((oligosfile != "") && (ecolifile != "")) {
                 m->mothurOut("[ERROR]: You can not use an ecoli file at the same time as an oligos file.\n"); abort = true;
            }
			
			//check to make sure you didn't forget the name file by mistake			
			if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "PcrSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int PcrSeqsCommand::execute(){
	try{
        
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        long start = time(NULL);
        fileAligned = true; pairedOligos = false;
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string trimSeqFile = getOutputFileName("fasta",variables);
		outputNames.push_back(trimSeqFile); outputTypes["fasta"].push_back(trimSeqFile);
        variables["[tag]"] = "scrap";
        string badSeqFile = getOutputFileName("fasta",variables);
		
		
        length = 0;
		if(oligosfile != ""){    readOligos();     if (m->getDebug()) { m->mothurOut("[DEBUG]: read oligos file. numprimers = " + toString(numFPrimers) + ", revprimers = " + toString(numRPrimers) + ".\n"); } }  if (m->getControl_pressed()) {  return 0; }
        if(ecolifile != "") {    readEcoli();      }  if (m->getControl_pressed()) {  return 0; }
        
        vector<unsigned long long> positions; 
        long long numFastaSeqs = 0;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = util.divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        if (processors == 1) {
            lines.push_back(linePair(0, 1000));
        }else {
            positions = m->setFilePosFasta(fastafile, numFastaSeqs); 
            if (numFastaSeqs < processors) { processors = numFastaSeqs; }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif
        if (m->getControl_pressed()) {  return 0; }

        set<string> badNames;
        numFastaSeqs = createProcesses(fastafile, trimSeqFile, badSeqFile, badNames);  	
		
		if (m->getControl_pressed()) {  return 0; }		
        
        //don't write or keep if blank
        if (badNames.size() != 0)   { writeAccnos(badNames);        }   
        if (util.isBlank(badSeqFile)) { util.mothurRemove(badSeqFile);  }
        else { outputNames.push_back(badSeqFile); outputTypes["fasta"].push_back(badSeqFile); }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        if (namefile != "")			{		readName(badNames);		}   
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        if (groupfile != "")		{		readGroup(badNames);    }
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		if (taxfile != "")			{		readTax(badNames);		}
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
  		if (countfile != "")			{		readCount(badNames);		}
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
      
        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}
		
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
        
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to screen " + toString(numFastaSeqs) + " sequences.");
		m->mothurOutEndLine();

		
		return 0;	
        
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int PcrSeqsCommand::createProcesses(string filename, string goodFileName, string badFileName, set<string>& badSeqNames) {
	try {
        
        vector<int> processIDS;   
        int process = 1;
		int num = 0;
        int pstart = -1; int pend = -1;
        bool adjustNeeded = false;
        bool recalc = false;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
                string locationsFile = toString(process) + ".temp";
				num = driverPcr(filename, goodFileName + toString(process) + ".temp", badFileName + toString(process) + ".temp", locationsFile, badSeqNames, lines[process], pstart, adjustNeeded);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename + toString(process) + ".num.temp";
				util.openOutputFile(tempFile, out);
                out << pstart << '\t' << adjustNeeded << endl;
				out << num << '\t' << badSeqNames.size() << endl;
                for (set<string>::iterator it = badSeqNames.begin(); it != badSeqNames.end(); it++) {
                    out << (*it) << endl;
                }
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->setControl_pressed(false);
                for (int i=0;i<processIDS.size();i++) {
                    util.mothurRemove(filename + (toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;

			}
		}
		
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->setControl_pressed(false);  for (int i=0;i<processIDS.size();i++) {util.mothurRemove(filename + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            lines.clear();
            vector<unsigned long long> positions = util.divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    string locationsFile = toString(process) + ".temp";
                    num = driverPcr(filename, goodFileName + toString(process) + ".temp", badFileName + toString(process) + ".temp", locationsFile, badSeqNames, lines[process], pstart, adjustNeeded);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = filename + toString(process) + ".num.temp";
                    util.openOutputFile(tempFile, out);
                    out << pstart << '\t' << adjustNeeded << endl;
                    out << num << '\t' << badSeqNames.size() << endl;
                    for (set<string>::iterator it = badSeqNames.begin(); it != badSeqNames.end(); it++) {
                        out << (*it) << endl;
                    }
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
            
        }

        string locationsFile = toString(process) + ".temp";
        num = driverPcr(filename, goodFileName, badFileName, locationsFile, badSeqNames, lines[0], pstart, adjustNeeded);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
			util.openInputFile(tempFile, in);
            int numBadNames = 0; string name = "";
            int tpstart = -1; bool tempAdjust = false;
            
			if (!in.eof()) {
                in >> tpstart >> tempAdjust; util.gobble(in);
                
                if (tempAdjust) { adjustNeeded = true; }
                if (tpstart != -1)   {
                    if (tpstart != pstart) { adjustNeeded = true; }
                    if (tpstart < pstart) { pstart = tpstart; } //smallest start
                } 
                int tempNum = 0; in >> tempNum >> numBadNames; num += tempNum; util.gobble(in);
            }
            for (int j = 0; j < numBadNames; j++) {
                in >> name; util.gobble(in);
                badSeqNames.insert(name);
            }
			in.close(); util.mothurRemove(tempFile);
            
            util.appendFiles((goodFileName + toString(processIDS[i]) + ".temp"), goodFileName);
            util.mothurRemove((goodFileName + toString(processIDS[i]) + ".temp"));
            
            util.appendFiles((badFileName + toString(processIDS[i]) + ".temp"), badFileName);
            util.mothurRemove((badFileName + toString(processIDS[i]) + ".temp"));
            
            util.appendFiles((toString(processIDS[i]) + ".temp"), locationsFile);
            util.mothurRemove((toString(processIDS[i]) + ".temp"));
		}
    #else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the sumScreenData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to badSeqNames.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<pcrData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
        string locationsFile = "locationsFile.txt";
        util.mothurRemove(locationsFile);
        util.mothurRemove(goodFileName);
        util.mothurRemove(badFileName);
        
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
            string extension = "";
            if (i!=0) {extension += toString(i) + ".temp"; processIDS.push_back(i); }
            
			// Allocate memory for thread data.
			pcrData* tempPcr = new pcrData(filename, goodFileName+extension, badFileName+extension, locationsFile+extension, m, oligosfile, ecolifile, nomatch, keepprimer, keepdots, start, end, length, pdiffs, rdiffs, lines[i].start, lines[i].end);
			pDataArray.push_back(tempPcr);
			
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyPcrThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        //do your part
        num = driverPcr(filename, (goodFileName+toString(processors-1)+".temp"), (badFileName+toString(processors-1)+".temp"), (locationsFile+toString(processors-1)+".temp"), badSeqNames, lines[processors-1], pstart, adjustNeeded);
        processIDS.push_back(processors-1);
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->fend) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->fend) + " sequences assigned to it, quitting. \n"); m->setControl_pressed(true); 
            }
            if (pDataArray[i]->adjustNeeded) { adjustNeeded = true; }
            if (pDataArray[i]->pstart != -1)   {
                if (pDataArray[i]->pstart != pstart) { adjustNeeded = true; }
                if (pDataArray[i]->pstart < pstart) { pstart = pDataArray[i]->pstart; }
            } //smallest start
            
            for (set<string>::iterator it = pDataArray[i]->badSeqNames.begin(); it != pDataArray[i]->badSeqNames.end(); it++) {	badSeqNames.insert(*it);       }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            util.appendFiles((goodFileName + toString(processIDS[i]) + ".temp"), goodFileName);
            util.mothurRemove((goodFileName + toString(processIDS[i]) + ".temp"));
            
            util.appendFiles((badFileName + toString(processIDS[i]) + ".temp"), badFileName);
            util.mothurRemove((badFileName + toString(processIDS[i]) + ".temp"));
            
            util.appendFiles((locationsFile+toString(processIDS[i]) + ".temp"), locationsFile);
            util.mothurRemove((locationsFile+toString(processIDS[i]) + ".temp"));
		}
        
#endif	
        
        

        if (fileAligned && adjustNeeded) {
            //find pend - pend is the biggest ending value, but we must account for when we adjust the start.  That adjustment may make the "new" end larger then the largest end. So lets find out what that "new" end will be.
            ifstream inLocations;
            util.openInputFile(locationsFile, inLocations);
            
            while(!inLocations.eof()) {
                
                if (m->getControl_pressed()) { break; }
                
                string name = "";
                int thisStart = -1; int thisEnd = -1;
                if (numFPrimers != 0)    { inLocations >> name >> thisStart; util.gobble(inLocations); }
                if (numRPrimers != 0)  { inLocations >> name >> thisEnd;   util.gobble(inLocations); }
                else { pend = -1; break; }
                
                int myDiff = 0;
                if (pstart != -1) {
                    if (thisStart != -1) {
                        if (thisStart != pstart) { myDiff += (thisStart - pstart); }
                    }
                }
                
                int myEnd = thisEnd + myDiff;
                //cout << name << '\t' << thisStart << '\t' << thisEnd << " diff = " << myDiff << '\t' << myEnd << endl;
                
                if (thisEnd != -1) { if (myEnd > pend) { pend = myEnd; } }
                
            }
            inLocations.close();
            
            adjustDots(goodFileName, locationsFile, pstart, pend);
        }else { util.mothurRemove(locationsFile); }
        
        return num;
        
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "createProcesses");
		exit(1);
	}
}

//**********************************************************************************************************************
int PcrSeqsCommand::driverPcr(string filename, string goodFasta, string badFasta, string locationsName, set<string>& badSeqNames, linePair filePos, int& pstart, bool& adjustNeeded){
	try {
		ofstream goodFile;
		util.openOutputFile(goodFasta, goodFile);
        
        ofstream badFile;
		util.openOutputFile(badFasta, badFile);
        
        ofstream locationsFile;
		util.openOutputFile(locationsName, locationsFile);
		
		ifstream inFASTA;
		util.openInputFile(filename, inFASTA);
        
		inFASTA.seekg(filePos.start);
        
		bool done = false;
		int count = 0;
        set<int> lengths;
        set<int> locations; //locations[0] = beginning locations, 
        
        //pdiffs, bdiffs, primers, barcodes, revPrimers
        map<string, int> primers;
        map<string, int> barcodes; //not used
        vector<string> revPrimer;
        if (pairedOligos) {
            map<int, oligosPair> primerPairs = oligos.getPairedPrimers();
            for (map<int, oligosPair>::iterator it = primerPairs.begin(); it != primerPairs.end(); it++) {
                primers[(it->second).forward] = it->first;
                revPrimer.push_back((it->second).reverse);
                //cout << (it->second).forward << '\t' << (it->second).reverse << endl;
            }
        }else{
            primers = oligos.getPrimers();
            revPrimer = oligos.getReversePrimers();
            //cout << (primers.begin())->first << '\t' << revPrimer[0] << endl;
        }
        
        TrimOligos trim(pdiffs, rdiffs, 0, primers, barcodes, revPrimer);
        
		while (!done) {
            
			if (m->getControl_pressed()) {  break; }
			
			Sequence currSeq(inFASTA); util.gobble(inFASTA);
            
            if (fileAligned) { //assume aligned until proven otherwise
                lengths.insert(currSeq.getAligned().length());
                if (lengths.size() > 1) { fileAligned = false; }
            }
            
            string trashCode = "";
            string locationsString = "";
            int thisPStart = -1;
            int thisPEnd = -1;
            int totalDiffs = 0;
            string commentString = "";
            
            if (m->getControl_pressed()) {  break; }
            
			if (currSeq.getName() != "") {
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: seq name = " + currSeq.getName() + ".\n"); } 
                
                bool goodSeq = true;
                if (oligosfile != "") {
                    map<int, int> mapAligned;
                    bool aligned = isAligned(currSeq.getAligned(), mapAligned);
                    
                    
                    //process primers
                    if (primers.size() != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        vector<int> results = trim.findForward(currSeq, primerStart, primerEnd);
                        bool good = true;
                        if (results[0] > pdiffs) { good = false; }
                        totalDiffs += results[0];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trim.getCodeValue(results[1], pdiffs) + ") ";
                        
                        if(!good){	if (nomatch == "reject") { goodSeq = false; } trashCode += "f";	}
                        else{
                            //are you aligned
                            if (aligned) { 
                                if (!keepprimer)    {
                                    if (keepdots)   { currSeq.filterToPos(mapAligned[primerEnd-1]+1);   } //mapAligned[primerEnd-1] is the location of the last base in the primer. we want to trim to the space just after that.  The -1 & +1 ensures if the primer is followed by gaps they are not trimmed causing an aligned sequence dataset to become unaligned.
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerEnd-1]+1));
                                        if (fileAligned) {
                                            thisPStart = mapAligned[primerEnd-1]+1; //locations[0].insert(mapAligned[primerEnd-1]+1);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerEnd-1]+1) + "\n";
                                        }
                                    }
                                } 
                                else                {  
                                    if (keepdots)   { currSeq.filterToPos(mapAligned[primerStart]);  }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerStart]));
                                        if (fileAligned) {
                                            thisPStart = mapAligned[primerStart]; //locations[0].insert(mapAligned[primerStart]);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerStart]) + "\n";
                                        }
                                    }
                                }
                                isAligned(currSeq.getAligned(), mapAligned);
                            }else { 
                                if (!keepprimer)    { currSeq.setAligned(currSeq.getUnaligned().substr(primerEnd)); } 
                                else                { currSeq.setAligned(currSeq.getUnaligned().substr(primerStart)); } 
                            }
                        }
                    }
                    
                    //process reverse primers
                    if (revPrimer.size() != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        vector<int> results = trim.findReverse(currSeq, primerStart, primerEnd);
                        bool good = true;
                        if (results[0] > rdiffs) { good = false; }
                        totalDiffs += results[0];
                        commentString += "rpdiffs=" + toString(results[0]) + "(" + trim.getCodeValue(results[1], rdiffs) + ") ";
                        
                        if(!good){	if (nomatch == "reject") { goodSeq = false; } trashCode += "r";	}
                        else{
                            //are you aligned
                            if (aligned) { 
                                if (!keepprimer)    {  
                                    if (keepdots)   { currSeq.filterFromPos(mapAligned[primerStart]); }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerStart]));
                                        if (fileAligned) {
                                            thisPEnd = mapAligned[primerStart]; //locations[1].insert(mapAligned[primerStart]);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerStart]) + "\n";
                                        }
                                    }
                                } 
                                else                {  
                                    if (keepdots)   { currSeq.filterFromPos(mapAligned[primerEnd-1]+1); }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerEnd-1]+1));
                                        if (fileAligned) {
                                            thisPEnd = mapAligned[primerEnd-1]+1; //locations[1].insert(mapAligned[primerEnd-1]+1);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerEnd-1]+1) + "\n";
                                        }
                                    }
                                } 
                            }
                            else { 
                                if (!keepprimer)    { currSeq.setAligned(currSeq.getUnaligned().substr(0, primerStart));   } 
                                else                { currSeq.setAligned(currSeq.getUnaligned().substr(0, primerEnd));     }
                            }
                        }
                    }
                }else if (ecolifile != "") {
                    //make sure the seqs are aligned
                    if (!fileAligned) { m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); m->setControl_pressed(true); break; }
                    else if (currSeq.getAligned().length() != length) {
                        m->mothurOut("[ERROR]: seqs are not the same length as ecoli seq. When using ecoli option your sequences must be aligned and the same length as the ecoli sequence.\n"); m->setControl_pressed(true); break; 
                    }else {
                        if (keepdots)   {
                            currSeq.filterFromPos(end);
                            currSeq.filterToPos(start-1);
                        }else {
                            string seqString = currSeq.getAligned().substr(0, end);
                            seqString = seqString.substr(start);
                            currSeq.setAligned(seqString); 
                        }
                    }
                }else{ //using start and end to trim
                    //make sure the seqs are aligned
                    if (!fileAligned) { m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); m->setControl_pressed(true); break; }
                    else {
                        if (end != -1) {
                            if (end > currSeq.getAligned().length()) {  m->mothurOut("[ERROR]: end is longer than your sequence length, aborting.\n"); m->setControl_pressed(true); break; }
                            else {
                                if (keepdots)   { currSeq.filterFromPos(end); }
                                else {
                                    string seqString = currSeq.getAligned().substr(0, end);
                                    currSeq.setAligned(seqString);
                                }
                            }
                        }
                        
                        if (start != -1) { 
                            if (keepdots)   {  currSeq.filterToPos(start-1);  }
                            else {
                                string seqString = currSeq.getAligned().substr(start);
                                currSeq.setAligned(seqString);
                            }
                        }
                    }
                }
                
                if (commentString != "") {
                    string seqComment = currSeq.getComment();
                    currSeq.setComment("\t" + commentString + "\t" + seqComment);
                }
                
                if (totalDiffs > (pdiffs + rdiffs)) { trashCode += "t"; goodSeq = false; }
                
                //trimming removed all bases
                if (currSeq.getUnaligned() == "") { goodSeq = false; }
                
				if(goodSeq == 1)    {
                    currSeq.printSequence(goodFile);
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + locationsString + "\n"); }
                    if (thisPStart != -1)   { locations.insert(thisPStart);  }
                    if (locationsString != "") { locationsFile << locationsString; }
                }
				else {  
                    badSeqNames.insert(currSeq.getName()); 
                    currSeq.setName(currSeq.getName() + '|' + trashCode);
                    currSeq.printSequence(badFile); 
                }
                count++;
			}
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = inFASTA.tellg();
            if ((pos == -1) || (pos >= filePos.end)) { break; }
#else
            if (inFASTA.eof()) { break; }
#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count)+"\n");		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count)+"\n"); 	}
		
        badFile.close();
		goodFile.close();
		inFASTA.close();
        locationsFile.close();
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: fileAligned = " + toString(fileAligned) +'\n'); }
    
        if (fileAligned && !keepdots) { //print out smallest start value and largest end value
            if (locations.size() > 1) { adjustNeeded = true; }
            if (primers.size() != 0)    {   set<int>::iterator it = locations.begin();  pstart = *it;  }
        }

		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "driverPcr");
		exit(1);
	}
}
//********************************************************************/
bool PcrSeqsCommand::isAligned(string seq, map<int, int>& aligned){
	try {
        aligned.clear();
        bool isAligned = false;
        
        int countBases = 0;
        for (int i = 0; i < seq.length(); i++) {
            if (!isalpha(seq[i])) { isAligned = true; }
            else { aligned[countBases] = i; countBases++; } //maps location in unaligned -> location in aligned.
        }                                                   //ie. the 3rd base may be at spot 10 in the alignment
                                                            //later when we trim we want to trim from spot 10.
        return isAligned;
    }
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "isAligned");
		exit(1);
	}
}
//**********************************************************************************************************************
int PcrSeqsCommand::adjustDots(string goodFasta, string locations, int pstart, int pend){
    try {
        ifstream inFasta;
        util.openInputFile(goodFasta, inFasta);
        
        ifstream inLocations;
        util.openInputFile(locations, inLocations);
        
        ofstream out;
        util.openOutputFile(goodFasta+".temp", out);
        
        set<int> lengths;
        //cout << pstart << '\t' << pend << endl;
        //if (pstart > pend) { //swap them
        
        while(!inFasta.eof()) {
            if(m->getControl_pressed()) { break; }
            
            Sequence seq(inFasta); util.gobble(inFasta);
            
            string name = "";
            int thisStart = -1; int thisEnd = -1;
            if (numFPrimers != 0)    { inLocations >> name >> thisStart; util.gobble(inLocations); }
            if (numRPrimers != 0)  { inLocations >> name >> thisEnd;   util.gobble(inLocations); }
            
            
            //cout << "'" << seq.getName() << "'" << '\t' << thisStart << '\t' << thisEnd << '\t' << seq.getAligned().length() << endl;
            //cout << "'" << name << "'" << '\t' << pstart << '\t' << pend << endl;
            
            if (name != seq.getName()) { m->mothurOut("[ERROR]: name mismatch in pcr.seqs.\n"); }
            else {
                if (pstart != -1) {
                    if (thisStart != -1) {
                        if (thisStart != pstart) {
                            string dots = "";
                            for (int i = pstart; i < thisStart; i++) { dots += "."; }
                            thisEnd += dots.length();
                            dots += seq.getAligned();
                            seq.setAligned(dots);
                        }
                    }
                }
                
                if (pend != -1) {
                    if (thisEnd != -1) {
                        if (thisEnd != pend) {
                            string dots = seq.getAligned();
                            for (int i = thisEnd; i < pend; i++) { dots += "."; }
                            seq.setAligned(dots);
                        }
                    }
                }
                lengths.insert(seq.getAligned().length());
            }
            
            seq.printSequence(out);
        }
        inFasta.close();
        inLocations.close();
        out.close();
        util.mothurRemove(locations);
        util.mothurRemove(goodFasta);
        util.renameFile(goodFasta+".temp", goodFasta);
        
        //cout << "final lengths = \n";
        //for (set<int>::iterator it = lengths.begin(); it != lengths.end(); it++) {
           //cout << *it << endl;
           // cout << lengths.count(*it) << endl;
       // }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "PcrSeqsCommand", "adjustDots");
        exit(1);
    }
}
//***************************************************************************************************************
bool PcrSeqsCommand::readEcoli(){
	try {
		ifstream in;
		util.openInputFile(ecolifile, in);
		
        //read seq
        if (!in.eof()){ 
            Sequence ecoli(in); 
            length = ecoli.getAligned().length();
            start = ecoli.getStartPos();
            end = ecoli.getEndPos();
        }else { in.close(); m->setControl_pressed(true); return false; }
        in.close();    
			
        return true;
    }
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readEcoli");
		exit(1);
	}
    
}
//***************************************************************************************************************
int PcrSeqsCommand::writeAccnos(set<string> badNames){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
		string outputFileName = getOutputFileName("accnos",variables);
        outputNames.push_back(outputFileName); outputTypes["accnos"].push_back(outputFileName);
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        for (set<string>::iterator it = badNames.begin(); it != badNames.end(); it++) {
            if (m->getControl_pressed()) { break; }
            out << (*it) << endl;
        }
        
        out.close();
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "writeAccnos");
		exit(1);
	}
    
}
//***************************************************************************************************************
int PcrSeqsCommand::readName(set<string>& names){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
        variables["[extension]"] = util.getExtension(namefile);
		string outputFileName = getOutputFileName("name", variables);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
			in >> firstCol;		util.gobble(in);		
			in >> secondCol;			
			
            string savedSecond = secondCol;
			vector<string> parsedNames;
			util.splitAtComma(secondCol, parsedNames);
			
			vector<string> validSecond;  validSecond.clear();
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
					validSecond.push_back(parsedNames[i]);
				}
			}
			
			if (validSecond.size() != parsedNames.size()) {  //we want to get rid of someone, so get rid of everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]);  }
				removedCount += parsedNames.size();
			}else {
                out << firstCol << '\t' << savedSecond << endl;
                wroteSomething = true;
            }
			util.gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["name"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your name file."); m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readName");
		exit(1);
	}
}
//**********************************************************************************************************************
int PcrSeqsCommand::readGroup(set<string> names){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(groupfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[extension]"] = util.getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
			in >> name;		util.gobble(in);		//read from first column
			in >> group;	util.gobble(in);		//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << group << endl;
			}else {  removedCount++;  }
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your group file."); m->mothurOutEndLine();
        
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int PcrSeqsCommand::readTax(set<string> names){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(taxfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
        variables["[extension]"] = util.getExtension(taxfile);
		string outputFileName = getOutputFileName("taxonomy", variables);

		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		ifstream in;
		util.openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
            in >> name; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << tax << endl;
			}else {  removedCount++;  }
            
        }
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["taxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your taxonomy file."); m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readTax");
		exit(1);
	}
}
//***************************************************************************************************************
int PcrSeqsCommand::readCount(set<string> badSeqNames){
	try {
		ifstream in;
		util.openInputFile(countfile, in);
		set<string>::iterator it;
		
		map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[extension]"] = util.getExtension(countfile);
		string goodCountFile = getOutputFileName("count", variables);

        outputNames.push_back(goodCountFile);  outputTypes["count"].push_back(goodCountFile);
		ofstream goodCountOut;	util.openOutputFile(goodCountFile, goodCountOut);
		
        string headers = util.getline(in); util.gobble(in);
        goodCountOut << headers << endl;
        string test = headers; vector<string> pieces = util.splitWhiteSpace(test);
        
        string name, rest; int thisTotal, removedCount; removedCount = 0; rest = "";
        bool wroteSomething = false;
        while (!in.eof()) {
            
			if (m->getControl_pressed()) { goodCountOut.close(); in.close(); util.mothurRemove(goodCountFile); return 0; }
            
			in >> name; util.gobble(in); 
            in >> thisTotal; util.gobble(in);
            if (pieces.size() > 2) {  rest = util.getline(in); util.gobble(in);  }
            
			if (badSeqNames.count(name) != 0) { removedCount+=thisTotal; }
			else{
                wroteSomething = true;
				goodCountOut << name << '\t' << thisTotal << '\t' << rest << endl;
			}
		}
		in.close();
		goodCountOut.close();
        
        if (m->getControl_pressed()) { util.mothurRemove(goodCountFile);   }
        
        if (wroteSomething == false) {  m->mothurOut("Your count file contains only sequences from the .accnos file."); m->mothurOutEndLine(); }
        
        //check for groups that have been eliminated
        CountTable ct;
        if (ct.testGroups(goodCountFile)) {
            ct.readTable(goodCountFile, true, false);
            ct.printTable(goodCountFile);
        }
		
		if (m->getControl_pressed()) { util.mothurRemove(goodCountFile);   }
        
        m->mothurOut("Removed " + toString(removedCount) + " sequences from your count file."); m->mothurOutEndLine();

		
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readCOunt");
		exit(1);
	}
}
//***************************************************************************************************************

int PcrSeqsCommand::readOligos(){
	try {
        oligos.read(oligosfile);
        
        if (m->getControl_pressed()) { return false; } //error in reading oligos
        
        if (oligos.hasPairedPrimers()) {
            pairedOligos = true;
            numFPrimers = oligos.getPairedPrimers().size();
            numRPrimers = numFPrimers;
        }else {
            pairedOligos = false;
            numFPrimers = oligos.getPrimers().size();
            numRPrimers = oligos.getReversePrimers().size();
        }
        
        if (oligos.getLinkers().size() != 0) { m->mothurOut("[WARNING]: pcr.seqs is not setup to remove linkers, ignoring.\n"); }
        if (oligos.getSpacers().size() != 0) { m->mothurOut("[WARNING]: pcr.seqs is not setup to remove spacers, ignoring.\n"); }

        return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readOligos");
		exit(1);
	}
}

/**************************************************************************************/


