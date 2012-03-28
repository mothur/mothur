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
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter poligos("oligos", "InputTypes", "", "", "ecolioligos", "none", "none",false,false); parameters.push_back(poligos);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
        CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
        CommandParameter ptax("taxonomy", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(ptax);
        CommandParameter pecoli("ecoli", "InputTypes", "", "", "ecolioligos", "none", "none",false,false); parameters.push_back(pecoli);
		CommandParameter pstart("start", "Number", "", "-1", "", "", "",false,false); parameters.push_back(pstart);
		CommandParameter pend("end", "Number", "", "-1", "", "", "",false,false); parameters.push_back(pend);
 		CommandParameter pnomatch("nomatch", "Multiple", "reject-keep", "reject", "", "", "",false,false); parameters.push_back(pnomatch);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ppdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pkeepprimer("keepprimer", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pkeepprimer);
        CommandParameter pkeepdots("keepdots", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pkeepdots);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
        
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
		helpString += "The pcr.seqs command reads a fasta file ...\n";
		
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

PcrSeqsCommand::PcrSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
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
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
                
                it = parameters.find("ecoli");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ecoli"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
			}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else if (fastafile == "not open") { fastafile = ""; abort = true; }	
			else { m->setFastaFile(fastafile); }
			
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "keepprimer", false);  if (temp == "not found")    {	temp = "f";	}
			keepprimer = m->isTrue(temp);	
            
            temp = validParameter.validFile(parameters, "keepdots", false);  if (temp == "not found")    {	temp = "t";	}
			keepdots = m->isTrue(temp);	
            
			temp = validParameter.validFile(parameters, "oligos", true);
			if (temp == "not found"){	oligosfile = "";		}
			else if(temp == "not open"){	oligosfile = ""; abort = true;	} 
			else					{	oligosfile = temp; m->setOligosFile(oligosfile);		}
			
            ecolifile = validParameter.validFile(parameters, "ecoli", true);
			if (ecolifile == "not found"){	ecolifile = "";		}
			else if(ecolifile == "not open"){	ecolifile = ""; abort = true;	} 
			
            namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not found"){	namefile = "";		}
			else if(namefile == "not open"){	namefile = ""; abort = true;	} 
            else { m->setNameFile(namefile); }
            
            groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not found"){	groupfile = "";		}
			else if(groupfile == "not open"){	groupfile = ""; abort = true;	} 
            else { m->setGroupFile(groupfile); }
            
            taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not found"){	taxfile = "";		}
			else if(taxfile == "not open"){	taxfile = ""; abort = true;	} 
            else { m->setTaxonomyFile(taxfile); }
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, pdiffs);
 			
			temp = validParameter.validFile(parameters, "start", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, start);
            
            temp = validParameter.validFile(parameters, "end", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, end);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
			
            nomatch = validParameter.validFile(parameters, "nomatch", false);	if (nomatch == "not found") { nomatch = "reject"; }
			
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
			if (namefile == "") {
				vector<string> files; files.push_back(fastafile);
				parser.getNameFile(files);
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
        
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        int start = time(NULL);
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string trimSeqFile = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "pcr.fasta";
		outputNames.push_back(trimSeqFile); outputTypes["fasta"].push_back(trimSeqFile);
        
        string badSeqFile = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "pcr.scrap.fasta";
		outputNames.push_back(badSeqFile); outputTypes["fasta"].push_back(badSeqFile);
		
        length = 0;
		if(oligosfile != ""){    readOligos();     }  if (m->control_pressed) {  return 0; }
        if(ecolifile != "") {    readEcoli();      }  if (m->control_pressed) {  return 0; }
        
        vector<unsigned long long> positions; 
        int numFastaSeqs = 0;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = m->divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        if (processors == 1) {
            lines.push_back(linePair(0, 1000));
        }else {
            positions = m->setFilePosFasta(fastafile, numFastaSeqs); 
            if (positions.size() < processors) { processors = positions.size(); }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif
        if (m->control_pressed) {  return 0; }

        set<string> badNames;
        if(processors == 1) {    numFastaSeqs = driverPcr(fastafile, trimSeqFile, badSeqFile, badNames, lines[0]);   }
        else                {    numFastaSeqs = createProcesses(fastafile, trimSeqFile, badSeqFile, badNames);       }	
		
		if (m->control_pressed) {  return 0; }		
        
        //don't write or keep if blank
        if (badNames.size() != 0)   { writeAccnos(badNames);        }   
        if (m->isBlank(badSeqFile)) { m->mothurRemove(badSeqFile);  }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        if (namefile != "")			{		readName(badNames);		}   
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        if (groupfile != "")		{		readGroup(badNames);    }
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
		if (taxfile != "")			{		readTax(badNames);		}
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
        
        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
		}
		
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTaxonomyFile(current); }
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
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverPcr(filename, goodFileName + toString(getpid()) + ".temp", badFileName + toString(getpid()) + ".temp", badSeqNames, lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
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
		
        num = driverPcr(filename, goodFileName, badFileName, badSeqNames, lines[0]);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
            int numBadNames = 0; string name = "";
			if (!in.eof()) { int tempNum = 0; in >> tempNum >> numBadNames; num += tempNum; m->gobble(in); }
            for (int j = 0; j < numBadNames; j++) {
                in >> name; m->gobble(in);
                badSeqNames.insert(name);
            }
			in.close(); m->mothurRemove(tempFile);
            
            m->appendFiles((goodFileName + toString(processIDS[i]) + ".temp"), goodFileName);
            m->mothurRemove((goodFileName + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((badFileName + toString(processIDS[i]) + ".temp"), badFileName);
            m->mothurRemove((badFileName + toString(processIDS[i]) + ".temp"));
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
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
            string extension = "";
            if (i!=0) {extension += toString(i) + ".temp"; processIDS.push_back(i); }
            
			// Allocate memory for thread data.
			pcrData* tempPcr = new pcrData(filename, goodFileName+extension, badFileName+extension, m, oligosfile, ecolifile, primers, revPrimer, nomatch, keepprimer, keepdots, start, end, length, lines[i].start, lines[i].end);
			pDataArray.push_back(tempPcr);
			
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyPcrThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        //do your part
        num = driverPcr(filename, (goodFileName+toString(processors-1)+".temp"), (badFileName+toString(processors-1)+".temp"),badSeqNames, lines[processors-1]);
        processIDS.push_back(processors-1);
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            for (set<string>::iterator it = pDataArray[i]->badSeqNames.begin(); it != pDataArray[i]->badSeqNames.end(); it++) {	badSeqNames.insert(*it);       }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            m->appendFiles((goodFileName + toString(processIDS[i]) + ".temp"), goodFileName);
            m->mothurRemove((goodFileName + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((badFileName + toString(processIDS[i]) + ".temp"), badFileName);
            m->mothurRemove((badFileName + toString(processIDS[i]) + ".temp"));
		}
        
#endif	
        
        return num;
        
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "createProcesses");
		exit(1);
	}
}

//**********************************************************************************************************************
int PcrSeqsCommand::driverPcr(string filename, string goodFasta, string badFasta, set<string>& badSeqNames, linePair filePos){
	try {
		ofstream goodFile;
		m->openOutputFile(goodFasta, goodFile);
        
        ofstream badFile;
		m->openOutputFile(badFasta, badFile);
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);
        
		inFASTA.seekg(filePos.start);
        
		bool done = false;
		int count = 0;
        set<int> lengths;
        
		while (!done) {
            
			if (m->control_pressed) {  break; }
			
			Sequence currSeq(inFASTA); m->gobble(inFASTA);
            
            string trashCode = "";
			if (currSeq.getName() != "") {
                
                bool goodSeq = true;
                if (oligosfile != "") {
                    map<int, int> mapAligned;
                    bool aligned = isAligned(currSeq.getAligned(), mapAligned);
                    
                    //process primers
                    if (primers.size() != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        bool good = findForward(currSeq, primerStart, primerEnd);
                        
                        if(!good){	if (nomatch == "reject") { goodSeq = false; } trashCode += "f";	}
                        else{
                            //are you aligned
                            if (aligned) { 
                                if (!keepprimer)    {  
                                    if (keepdots)   { currSeq.filterToPos(mapAligned[primerEnd]);   }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerEnd]));                                              }
                                } 
                                else                {  
                                    if (keepdots)   { currSeq.filterToPos(mapAligned[primerStart]);  }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerStart]));                                              }
                                }
                            }else { 
                                if (!keepprimer)    { currSeq.setAligned(currSeq.getUnaligned().substr(primerEnd)); } 
                                else                { currSeq.setAligned(currSeq.getUnaligned().substr(primerStart)); } 
                            }
                        }
                    }
                    
                    //process reverse primers
                    if (revPrimer.size() != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        bool good = findReverse(currSeq, primerStart, primerEnd);
                        if(!good){	if (nomatch == "reject") { goodSeq = false; } trashCode += "r";	}
                        else{ 
                            //are you aligned
                            if (aligned) { 
                                if (!keepprimer)    {  
                                    if (keepdots)   { currSeq.filterFromPos(mapAligned[primerStart]); }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerStart]));   }
                                } 
                                else                {  
                                    if (keepdots)   { currSeq.filterFromPos(mapAligned[primerEnd]); }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerEnd]));   }
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
                    lengths.insert(currSeq.getAligned().length());
                    if (lengths.size() > 1) { m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); m->control_pressed = true; break; }
                    else if (currSeq.getAligned().length() != length) {
                        m->mothurOut("[ERROR]: seqs are not the same length as ecoli seq. When using ecoli option your sequences must be aligned and the same length as the ecoli sequence.\n"); m->control_pressed = true; break; 
                    }else {
                        if (keepdots)   { 
                            currSeq.filterToPos(start); 
                            currSeq.filterFromPos(end);
                        }else {
                            string seqString = currSeq.getAligned().substr(0, end);
                            seqString = seqString.substr(start);
                            currSeq.setAligned(seqString); 
                        }
                    }
                }else{ //using start and end to trim
                    //make sure the seqs are aligned
                    lengths.insert(currSeq.getAligned().length());
                    if (lengths.size() > 1) { m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); m->control_pressed = true; break; }
                    else {
                        if (end != -1) {
                            if (end > currSeq.getAligned().length()) {  m->mothurOut("[ERROR]: end is longer than your sequence length, aborting.\n"); m->control_pressed = true; break; }
                            else {
                                if (keepdots)   { currSeq.filterFromPos(end); }
                                else {
                                    string seqString = currSeq.getAligned().substr(0, end);
                                    currSeq.setAligned(seqString); 
                                }
                            }
                        }
                        if (start != -1) { 
                            if (keepdots)   {  currSeq.filterToPos(start);  }
                            else {
                                string seqString = currSeq.getAligned().substr(start);
                                currSeq.setAligned(seqString); 
                            }
                        }
                    }
                }
                                    
				if(goodSeq == 1)    {   currSeq.printSequence(goodFile);        }
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
			if((count) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		
        badFile.close();
		goodFile.close();
		inFASTA.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "driverPcr");
		exit(1);
	}
}
//********************************************************************/
bool PcrSeqsCommand::findForward(Sequence& seq, int& primerStart, int& primerEnd){
	try {
		
		string rawSequence = seq.getUnaligned();
		
		for(int j=0;j<primers.size();j++){
			string oligo = primers[j];
			
			if(rawSequence.length() < oligo.length()) {  break;  }
			
			//search for primer
            int olength = oligo.length();
            for (int j = 0; j < rawSequence.length()-olength; j++){
                if (m->control_pressed) {  primerStart = 0; primerEnd = 0; return false; }
                string rawChunk = rawSequence.substr(j, olength);
                if(compareDNASeq(oligo, rawChunk)) {
                    primerStart = j;
                    primerEnd = primerStart + olength;
                    return true;
                }
                
            }
        }	
		
        primerStart = 0; primerEnd = 0;
		return false;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "stripForward");
		exit(1);
	}
}
//******************************************************************/
bool PcrSeqsCommand::findReverse(Sequence& seq, int& primerStart, int& primerEnd){
	try {
		
		string rawSequence = seq.getUnaligned();
		
		for(int i=0;i<revPrimer.size();i++){
			string oligo = revPrimer[i];
			if(rawSequence.length() < oligo.length()) {  break;  }
			
			//search for primer
            int olength = oligo.length();
            for (int j = rawSequence.length()-olength; j >= 0; j--){
                 if (m->control_pressed) {  primerStart = 0; primerEnd = 0; return false; }
                string rawChunk = rawSequence.substr(j, olength);
            
                if(compareDNASeq(oligo, rawChunk)) {
                    primerStart = j;
                    primerEnd = primerStart + olength;
                    return true;
                }
                
            }
		}	
		
        primerStart = 0; primerEnd = 0;
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "findReverse");
		exit(1);
	}
}
//********************************************************************/
bool PcrSeqsCommand::isAligned(string seq, map<int, int>& aligned){
	try {
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
//********************************************************************/
string PcrSeqsCommand::reverseOligo(string oligo){
	try {
        string reverse = "";
       
        for(int i=oligo.length()-1;i>=0;i--){
            
            if(oligo[i] == 'A')		{	reverse += 'T';	}
            else if(oligo[i] == 'T'){	reverse += 'A';	}
            else if(oligo[i] == 'U'){	reverse += 'A';	}
            
            else if(oligo[i] == 'G'){	reverse += 'C';	}
            else if(oligo[i] == 'C'){	reverse += 'G';	}
            
            else if(oligo[i] == 'R'){	reverse += 'Y';	}
            else if(oligo[i] == 'Y'){	reverse += 'R';	}
            
            else if(oligo[i] == 'M'){	reverse += 'K';	}
            else if(oligo[i] == 'K'){	reverse += 'M';	}
            
            else if(oligo[i] == 'W'){	reverse += 'W';	}
            else if(oligo[i] == 'S'){	reverse += 'S';	}
            
            else if(oligo[i] == 'B'){	reverse += 'V';	}
            else if(oligo[i] == 'V'){	reverse += 'B';	}
            
            else if(oligo[i] == 'D'){	reverse += 'H';	}
            else if(oligo[i] == 'H'){	reverse += 'D';	}

            else						{	reverse += 'N';	}
        }

        
        return reverse;
    }
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "reverseOligo");
		exit(1);
	}
}

//***************************************************************************************************************
bool PcrSeqsCommand::readOligos(){
	try {
		ifstream inOligos;
		m->openInputFile(oligosfile, inOligos);
		
		string type, oligo, group;
		
		while(!inOligos.eof()){
            
			inOligos >> type; 
            
			if(type[0] == '#'){ //ignore
				while (!inOligos.eof())	{	char c = inOligos.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				m->gobble(inOligos);
			}else{
				m->gobble(inOligos);
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				inOligos >> oligo;
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{	
                        char c = inOligos.get(); 
                        if (c == 10 || c == 13){	break;	} 
                        else if (c == 32 || c == 9){;} //space or tab
					} 
					primers.push_back(oligo);
                }else if(type == "REVERSE"){
                    string oligoRC = reverseOligo(oligo);
                    revPrimer.push_back(oligoRC);
                    //cout << "oligo = " << oligo << " reverse = " << oligoRC << endl;
				}else if(type == "BARCODE"){
					inOligos >> group;
				}else if((type == "LINKER")||(type == "SPACER")) {;}
				else{	m->mothurOut(type + " is not recognized as a valid type. Choices are forward, reverse, linker, spacer and barcode. Ignoring " + oligo + "."); m->mothurOutEndLine(); m->control_pressed = true; }
			}
			m->gobble(inOligos);
		}	
		inOligos.close();
		
		if ((primers.size() == 0) && (revPrimer.size() == 0)) {
			m->mothurOut("[ERROR]: your oligos file does not contain valid primers or reverse primers.  Please correct."); m->mothurOutEndLine();
            m->control_pressed = true;
			return false;
		}
        
        return true;
        
    }catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readOligos");
		exit(1);
	}
}
//***************************************************************************************************************
bool PcrSeqsCommand::readEcoli(){
	try {
		ifstream in;
		m->openInputFile(ecolifile, in);
		
        //read seq
        if (!in.eof()){ 
            Sequence ecoli(in); 
            length = ecoli.getAligned().length();
            start = ecoli.getStartPos();
            end = ecoli.getEndPos();
        }else { in.close(); m->control_pressed = true; return false; }
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "bad.accnos";
        outputNames.push_back(outputFileName); outputTypes["accnos"].push_back(outputFileName);
        
        ofstream out;
        m->openOutputFile(outputFileName, out);
        
        for (set<string>::iterator it = badNames.begin(); it != badNames.end(); it++) {
            if (m->control_pressed) { break; }
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
//******************************************************************/
bool PcrSeqsCommand::compareDNASeq(string oligo, string seq){
	try {
		bool success = 1;
		int length = oligo.length();
		
		for(int i=0;i<length;i++){
			
			if(oligo[i] != seq[i]){
				if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C')	{	success = 0; 	}
				else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))				{	success = 0;	}
				else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))					{	success = 0;	}
				else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))					{	success = 0;	}
				else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))					{	success = 0;	}
				else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
				else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
				else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))	{	success = 0;	}
				else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))	{	success = 0;	}			
				
				if(success == 0)	{	break;	 }
			}
			else{
				success = 1;
			}
		}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "compareDNASeq");
		exit(1);
	}
	
}
//***************************************************************************************************************
int PcrSeqsCommand::readName(set<string>& names){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(namefile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(namefile)) + "pcr" + m->getExtension(namefile);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			in >> firstCol;		m->gobble(in);		
			in >> secondCol;			
			
            string savedSecond = secondCol;
			vector<string> parsedNames;
			m->splitAtComma(secondCol, parsedNames);
			
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
			m->gobble(in);
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "pcr" + m->getExtension(groupfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << group << endl;
			}else {  removedCount++;  }
            
			m->gobble(in);
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
		if (outputDir == "") {  thisOutputDir += m->hasPath(taxfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(taxfile)) + "pcr" + m->getExtension(taxfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			in >> name;				//read from first column
			in >> tax;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << tax << endl;
			}else {  removedCount++;  }
            
			m->gobble(in);
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
/**************************************************************************************/


