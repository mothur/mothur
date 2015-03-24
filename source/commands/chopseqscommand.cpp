/*
 *  chopseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chopseqscommand.h"
#include "sequence.hpp"
#include "removeseqscommand.h"

//**********************************************************************************************************************
vector<string> ChopSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false,true); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pnumbases("numbases", "Number", "", "0", "", "", "","",false,true,true); parameters.push_back(pnumbases);
		CommandParameter pcountgaps("countgaps", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pcountgaps);
		CommandParameter pshort("short", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pshort);
		CommandParameter pkeep("keep", "Multiple", "front-back", "front", "", "", "","",false,false); parameters.push_back(pkeep);
        CommandParameter pkeepn("keepn", "Boolean", "", "f", "", "", "","",false,false); parameters.push_back(pkeepn);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChopSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chop.seqs command reads a fasta file and outputs a .chop.fasta containing the trimmed sequences. Note: If a sequence is completely 'chopped', an accnos file will be created with the names of the sequences removed. \n";
		helpString += "The chop.seqs command parameters are fasta, name, group, count, numbases, countgaps and keep. fasta is required unless you have a valid current fasta file. numbases is required.\n";
		helpString += "The chop.seqs command should be in the following format: chop.seqs(fasta=yourFasta, numbases=yourNum, keep=yourKeep).\n";
        helpString += "If you provide a name, group or count file any sequences removed from the fasta file will also be removed from those files.\n";
		helpString += "The numbases parameter allows you to specify the number of bases you want to keep.\n";
		helpString += "The keep parameter allows you to specify whether you want to keep the front or the back of your sequence, default=front.\n";
		helpString += "The countgaps parameter allows you to specify whether you want to count gaps as bases, default=false.\n";
		helpString += "The short parameter allows you to specify you want to keep sequences that are too short to chop, default=false.\n";
        helpString += "The keepn parameter allows you to specify you want to keep ambigous bases, default=false.\n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "For example, if you ran chop.seqs with numbases=200 and short=t, if a sequence had 100 bases mothur would keep the sequence rather than eliminate it.\n";
		helpString += "Example chop.seqs(fasta=amazon.fasta, numbases=200, keep=front).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChopSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],chop.fasta"; }
        else if (type == "name") {  pattern = "[filename],chop.names"; }
        else if (type == "group") {  pattern = "[filename],chop.groups"; }
        else if (type == "count") {  pattern = "[filename],chop.count_table"; } 
        else if (type == "accnos") {  pattern = "[filename],chop.accnos"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChopSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChopSeqsCommand::ChopSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "ChopSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ChopSeqsCommand::ChopSeqsCommand(string option)  {
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
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
            outputTypes["name"] = tempOutNames;
            outputTypes["group"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
		
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  				//if there is a current fasta file, use it
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); } 	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            			
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }
			else { m->setCountTableFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			string temp = validParameter.validFile(parameters, "numbases", false);	if (temp == "not found") { temp = "0"; } 
			m->mothurConvert(temp, numbases);   
			
            temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
			temp = validParameter.validFile(parameters, "countgaps", false);	if (temp == "not found") { temp = "f"; } 
			countGaps = m->isTrue(temp);  
			
			temp = validParameter.validFile(parameters, "short", false);	if (temp == "not found") { temp = "f"; } 
			Short = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "keepn", false);	if (temp == "not found") { temp = "f"; }
            keepN = m->isTrue(temp);
		
			keep = validParameter.validFile(parameters, "keep", false);		if (keep == "not found") { keep = "front"; } 
				
			if (numbases == 0)  { m->mothurOut("You must provide the number of bases you want to keep for the chops.seqs command."); m->mothurOutEndLine(); abort = true;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "ChopSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChopSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> variables;
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(fastafile));
        string outputFileName = getOutputFileName("fasta", variables);
        outputNames.push_back(outputFileName); outputTypes["fasta"].push_back(outputFileName);
        string outputFileNameAccnos = getOutputFileName("accnos", variables);        
        
        vector<unsigned long long> positions; 
        vector<linePair> lines;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = m->divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        int numSeqs = 0;
        positions = m->setFilePosFasta(fastafile, numSeqs); 
        if (positions.size() < processors) { processors = positions.size(); }
		
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        bool wroteAccnos = false;
        if(processors == 1) {   wroteAccnos = driver(lines[0], fastafile, outputFileName, outputFileNameAccnos);        }
        else                {   wroteAccnos = createProcesses(lines, fastafile, outputFileName, outputFileNameAccnos);  }
        
        if (m->control_pressed) {  return 0; }
		
        if (wroteAccnos) {
            outputNames.push_back(outputFileNameAccnos); outputTypes["accnos"].push_back(outputFileNameAccnos);
            
             //use remove.seqs to create new name, group and count file
            if ((countfile != "") || (namefile != "") || (groupfile != "")) {
                string inputString = "accnos=" + outputFileNameAccnos;
                
                if (countfile != "") {  inputString += ", count=" + countfile;  }
                else{
                    if (namefile != "") {  inputString += ", name=" + namefile;  }
                    if (groupfile != "") {  inputString += ", group=" + groupfile;  }
                }
                
                m->mothurOut("/******************************************/"); m->mothurOutEndLine();
                m->mothurOut("Running command: remove.seqs(" + inputString + ")"); m->mothurOutEndLine();
                m->mothurCalling = true;
                
                Command* removeCommand = new RemoveSeqsCommand(inputString);
                removeCommand->execute();
                
                map<string, vector<string> > filenames = removeCommand->getOutputFiles();
                
                delete removeCommand;
                m->mothurCalling = false;
                m->mothurOut("/******************************************/"); m->mothurOutEndLine();
                
                if (groupfile != "") {
                    thisOutputDir = outputDir;
                    if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
                    variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(groupfile));
                    string outGroup = getOutputFileName("group", variables);
                    m->renameFile(filenames["group"][0], outGroup);
                    outputNames.push_back(outGroup); outputTypes["group"].push_back(outGroup);
                }
                
                if (namefile != "") {
                    thisOutputDir = outputDir;
                    if (outputDir == "") {  thisOutputDir += m->hasPath(namefile);  }
                    variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(namefile));
                    string outName = getOutputFileName("name", variables);
                    m->renameFile(filenames["name"][0], outName);
                    outputNames.push_back(outName); outputTypes["name"].push_back(outName);
                }
                
                if (countfile != "") {
                    thisOutputDir = outputDir;
                    if (outputDir == "") {  thisOutputDir += m->hasPath(countfile);  }
                    variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(countfile));
                    string outCount = getOutputFileName("count", variables);
                    m->renameFile(filenames["count"][0], outCount);
                    outputNames.push_back(outCount); outputTypes["count"].push_back(outCount);
                }
            }
        }
		else {  m->mothurRemove(outputFileNameAccnos);  }
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
        
		if (wroteAccnos) { //set accnos file as new current accnosfile
			itTypes = outputTypes.find("accnos");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
			}
            
            itTypes = outputTypes.find("name");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
            }
            
            itTypes = outputTypes.find("group");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
            }
            
            itTypes = outputTypes.find("count");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
            }
		}
		
        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
bool ChopSeqsCommand::createProcesses(vector<linePair> lines, string filename, string outFasta, string outAccnos) {
	try {
		int process = 1;
		bool wroteAccnos = false;
		vector<int> processIDS;
        vector<string> nonBlankAccnosFiles;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				wroteAccnos = driver(lines[process], filename, outFasta + m->mothurGetpid(process) + ".temp", outAccnos + m->mothurGetpid(process) + ".temp");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = fastafile + m->mothurGetpid(process) + ".bool.temp";
				m->openOutputFile(tempFile, out);
				out << wroteAccnos << endl;				
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do your part
		wroteAccnos = driver(lines[0], filename, outFasta, outAccnos);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
        
		if (wroteAccnos) { nonBlankAccnosFiles.push_back(outAccnos); }
		else { m->mothurRemove(outAccnos); } //remove so other files can be renamed to it
        
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = fastafile + toString(processIDS[i]) + ".bool.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			bool temp;
			in >> temp; m->gobble(in); 
            if (temp) { wroteAccnos = temp; nonBlankAccnosFiles.push_back(outAccnos + toString(processIDS[i]) + ".temp");  }
			else { m->mothurRemove((outAccnos + toString(processIDS[i]) + ".temp"));  }
            
			in.close();
			m->mothurRemove(tempFilename);
		}
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<chopData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
            string extension = "";
            if (i != 0) { extension = toString(i) + ".temp"; processIDS.push_back(i); }
			// Allocate memory for thread data.
			chopData* tempChop = new chopData(filename, (outFasta+extension), (outAccnos+extension), m, lines[i].start, lines[i].end, keep, countGaps, numbases, Short, keepN);
			pDataArray.push_back(tempChop);
			
			//MyChopThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyChopThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        //do your part
		wroteAccnos = driver(lines[processors-1], filename, (outFasta + toString(processors-1) + ".temp"), (outAccnos + toString(processors-1) + ".temp"));
        processIDS.push_back(processors-1);
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
        if (wroteAccnos) { nonBlankAccnosFiles.push_back(outAccnos); }
		else { m->mothurRemove(outAccnos); } //remove so other files can be renamed to it

		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            if (pDataArray[i]->wroteAccnos) { wroteAccnos = pDataArray[i]->wroteAccnos; nonBlankAccnosFiles.push_back(outAccnos + toString(processIDS[i]) + ".temp");  }
			else { m->mothurRemove((outAccnos + toString(processIDS[i]) + ".temp"));  }
            //check to make sure the process finished
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif		
                
		for (int i = 0; i < processIDS.size(); i++) {
			m->appendFiles((outFasta + toString(processIDS[i]) + ".temp"), outFasta);
			m->mothurRemove((outFasta + toString(processIDS[i]) + ".temp"));
		}
		
        if (nonBlankAccnosFiles.size() != 0) { 
			m->renameFile(nonBlankAccnosFiles[0], outAccnos);
			
			for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
				m->appendFiles(nonBlankAccnosFiles[h], outAccnos);
				m->mothurRemove(nonBlankAccnosFiles[h]);
			}
		}else { //recreate the accnosfile if needed
			ofstream out;
			m->openOutputFile(outAccnos, out);
			out.close();
		}

		return wroteAccnos;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************/
bool ChopSeqsCommand::driver(linePair filePos, string filename, string outFasta, string outAccnos) {	
	try {
		
		ofstream out;
		m->openOutputFile(outFasta, out);
        
        ofstream outAcc;
		m->openOutputFile(outAccnos, outAcc);
        
		ifstream in;
		m->openInputFile(filename, in);
        
		in.seekg(filePos.start);
        
        //adjust
        if (filePos.start == 0) {
            m->zapGremlins(in); m->gobble(in);
        }
        
		bool done = false;
        bool wroteAccnos = false;
		int count = 0;
        
		while (!done) {
            
			if (m->control_pressed) { in.close(); out.close(); return 1; }
            
			Sequence seq(in); m->gobble(in);
			
			if (m->control_pressed) {  in.close(); out.close(); outAcc.close(); m->mothurRemove(outFasta); m->mothurRemove(outAccnos); return 0;  }
			
			if (seq.getName() != "") {
				string newSeqString = getChopped(seq);
				
				//output trimmed sequence
				if (newSeqString != "") {
					out << ">" << seq.getName() << endl << newSeqString << endl;
				}else{
					outAcc << seq.getName() << endl;
					wroteAccnos = true;
				}
                count++;
			}
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= filePos.end)) { break; }
#else
            if (in.eof()) { break; }
#endif
            //report progress
			if((count) % 1000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 1000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}

		
		in.close();
        out.close();
        outAcc.close();
		
		return wroteAccnos;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChopSeqsCommand::getChopped(Sequence seq) {
	try {
		string temp = seq.getAligned();
		string tempUnaligned = seq.getUnaligned();
		
		if (countGaps) {
			//if needed trim sequence
			if (keep == "front") {//you want to keep the beginning
				int tempLength = temp.length();

				if (tempLength > numbases) { //you have enough bases to remove some
				
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = 0; i < temp.length(); i++) {
						//eliminate N's
                        if (!keepN) { if (toupper(temp[i]) == 'N') { temp[i] = '.'; } }
						
						numBasesCounted++; 
						
						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
					
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(0, stopSpot+1);  }
							
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}
			}else { //you are keeping the back
				int tempLength = temp.length();
				if (tempLength > numbases) { //you have enough bases to remove some
					
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = (temp.length()-1); i >= 0; i--) {
						//eliminate N's
                        if (!keepN) { if (toupper(temp[i]) == 'N') { temp[i] = '.'; } }
						
						numBasesCounted++; 

						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
				
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(stopSpot+1);  }
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}
			}

		}else{
				
			//if needed trim sequence
			if (keep == "front") {//you want to keep the beginning
				int tempLength = tempUnaligned.length();

				if (tempLength > numbases) { //you have enough bases to remove some
					
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = 0; i < temp.length(); i++) {
						//eliminate N's
                        if (!keepN) {
                            if (toupper(temp[i]) == 'N') {
                                temp[i] = '.';
                                tempLength--;
                                if (tempLength < numbases) { stopSpot = 0; break; }
                            }
                        }
						if(isalpha(temp[i])) { numBasesCounted++; }
						
						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
					
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(0, stopSpot+1);  }
							
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}				
			}else { //you are keeping the back
				int tempLength = tempUnaligned.length();
				if (tempLength > numbases) { //you have enough bases to remove some
					
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = (temp.length()-1); i >= 0; i--) {
                        if (!keepN) {
                            //eliminate N's
                            if (toupper(temp[i]) == 'N') {
                                temp[i] = '.';
                                tempLength--;
                                if (tempLength < numbases) { stopSpot = 0; break; }
                            }
                        }
						if(isalpha(temp[i])) { numBasesCounted++; }

						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
				
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(stopSpot);  }
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}
			}
		}
		
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "getChopped");
		exit(1);
	}
}
//**********************************************************************************************************************


