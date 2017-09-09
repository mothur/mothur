/*
 *  filterseqscommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "filterseqscommand.h"
#include "sequence.hpp"


//**********************************************************************************************************************
vector<string> FilterSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-filter",false,true, true); parameters.push_back(pfasta);
		CommandParameter phard("hard", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(phard);
		CommandParameter ptrump("trump", "String", "", "*", "", "", "","",false,false, true); parameters.push_back(ptrump);
		CommandParameter psoft("soft", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psoft);
		CommandParameter pvertical("vertical", "Boolean", "", "T", "", "", "","",false,false, true); parameters.push_back(pvertical);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false, true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string FilterSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The filter.seqs command reads a file containing sequences and creates a .filter and .filter.fasta file.\n";
		helpString += "The filter.seqs command parameters are fasta, trump, soft, hard, processors and vertical. \n";
		helpString += "The fasta parameter is required, unless you have a valid current fasta file. You may enter several fasta files to build the filter from and filter, by separating their names with -'s.\n";
		helpString += "For example: fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The trump option will remove a column if the trump character is found at that position in any sequence of the alignment. Default=*, meaning no trump. \n";
		helpString += "A soft mask removes any column where the dominant base (i.e. A, T, G, C, or U) does not occur in at least a designated percentage of sequences. Default=0.\n";
		helpString += "The hard parameter allows you to enter a file containing the filter you want to use.\n";
		helpString += "The vertical parameter removes columns where all sequences contain a gap character. The default is T.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The filter.seqs command should be in the following format: \n";
		helpString += "filter.seqs(fasta=yourFastaFile, trump=yourTrump) \n";
		helpString += "Example filter.seqs(fasta=abrecovery.fasta, trump=.).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string FilterSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],filter.fasta"; } 
        else if (type == "filter") {  pattern =  "[filename],filter"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "FilterSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
FilterSeqsCommand::FilterSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["filter"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "FilterSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
FilterSeqsCommand::FilterSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;  recalced = false;
		filterFileName = "";
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("filter.seqs");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["filter"] = tempOutNames;
		
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
				
				it = parameters.find("hard");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["hard"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fasta = validParameter.validFile(parameters, "fasta", false);
			if (fasta == "not found") { 				
				fasta = m->getFastaFile(); 
				if (fasta != "") { 
                    fastafileNames.push_back(fasta);  
                    m->mothurOut("Using " + fasta + " as input file for the fasta parameter."); m->mothurOutEndLine();
                    string simpleName = m->getSimpleName(fasta);
                    filterFileName += simpleName.substr(0, simpleName.find_first_of('.'));
                }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else { 
				m->splitAtDash(fasta, fastafileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastafileNames.size(); i++) {
					
					bool ignore = false;
					if (fastafileNames[i] == "current") { 
						fastafileNames[i] = m->getFastaFile(); 
						if (fastafileNames[i] != "") {  m->mothurOut("Using " + fastafileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastafileNames.erase(fastafileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						if (inputDir != "") {
							string path = m->hasPath(fastafileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	fastafileNames[i] = inputDir + fastafileNames[i];		}
						}

						ifstream in;
						bool ableToOpen = m->openInputFile(fastafileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(fastafileNames[i]);
								m->mothurOut("Unable to open " + fastafileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastafileNames[i] = tryPath;
							}
						}
						
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(fastafileNames[i]);
								m->mothurOut("Unable to open " + fastafileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastafileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (!ableToOpen) { 
							m->mothurOut("Unable to open " + fastafileNames[i] + ". It will be disregarded."); m->mothurOutEndLine();
							//erase from file list
							fastafileNames.erase(fastafileNames.begin()+i);
							i--;
						}else{  
							string simpleName = m->getSimpleName(fastafileNames[i]);
							filterFileName += simpleName.substr(0, simpleName.find_first_of('.'));
							m->setFastaFile(fastafileNames[i]);
						}
						in.close();
					}
				}
				
				//make sure there is at least one valid file left
				if (fastafileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (!abort) {
				//if the user changes the output directory command factory will send this info to us in the output parameter 
				outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
					outputDir = "";	
					outputDir += m->hasPath(fastafileNames[0]); //if user entered a file with a path then preserve it	
				}
			}
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			hard = validParameter.validFile(parameters, "hard", true);				if (hard == "not found") { hard = ""; }
			else if (hard == "not open") { hard = ""; abort = true; }	

			temp = validParameter.validFile(parameters, "trump", false);			if (temp == "not found") { temp = "*"; }
			trump = temp[0];
			
			temp = validParameter.validFile(parameters, "soft", false);				if (temp == "not found") { soft = 0; }
			else {  soft = (float)atoi(temp.c_str()) / 100.0;  }
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
			
			vertical = validParameter.validFile(parameters, "vertical", false);		
			if (vertical == "not found") { 
				if ((hard == "") && (trump == '*') && (soft == 0)) { vertical = "T"; } //you have not given a hard file or set the trump char.
				else { vertical = "F";  }
			}
			
			numSeqs = 0;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "FilterSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/

int FilterSeqsCommand::execute() {	
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		ifstream inFASTA;
		m->openInputFile(fastafileNames[0], inFASTA);
		
		Sequence testSeq(inFASTA);
		alignmentLength = testSeq.getAlignLength();
		inFASTA.close();
		
		////////////create filter/////////////////
		m->mothurOut("Creating Filter... "); m->mothurOutEndLine();
		
		filter = createFilter();
		
		m->mothurOutEndLine();  m->mothurOutEndLine();
		
		if (m->control_pressed) { outputTypes.clear(); return 0; }
		
		ofstream outFilter;
		
		//prevent giantic file name
        map<string, string> variables; 
        variables["[filename]"] = outputDir + filterFileName + ".";
		if (fastafileNames.size() > 3) { variables["[filename]"] = outputDir + "merge."; }
		string filterFile = getOutputFileName("filter", variables);  
		
		m->openOutputFile(filterFile, outFilter);
		outFilter << filter << endl;
		outFilter.close();
		outputNames.push_back(filterFile); outputTypes["filter"].push_back(filterFile);
				
		////////////run filter/////////////////
		
		m->mothurOut("Running Filter... "); m->mothurOutEndLine();
		
		filterSequences();
		
		m->mothurOutEndLine();	m->mothurOutEndLine();
					
		int filteredLength = 0;
		for(int i=0;i<alignmentLength;i++){
			if(filter[i] == '1'){	filteredLength++;	}
		}
		
		if (m->control_pressed) {  outputTypes.clear(); for(int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }

		
		m->mothurOutEndLine();
		m->mothurOut("Length of filtered alignment: " + toString(filteredLength)); m->mothurOutEndLine();
		m->mothurOut("Number of columns removed: " + toString((alignmentLength-filteredLength))); m->mothurOutEndLine();
		m->mothurOut("Length of the original alignment: " + toString(alignmentLength)); m->mothurOutEndLine();
		m->mothurOut("Number of sequences used to construct filter: " + toString(numSeqs)); m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
int FilterSeqsCommand::filterSequences() {	
	try {
		
		numSeqs = 0;
		
		for (int s = 0; s < fastafileNames.size(); s++) {
			
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
				
                map<string, string> variables; 
                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafileNames[s]));
				string filteredFasta = getOutputFileName("fasta", variables);
            
            vector<unsigned long long> positions;
            if (savedPositions.size() != 0) { positions = savedPositions[s]; }
            else {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				positions = m->divideFile(fastafileNames[s], processors);
#else
                if(processors != 1){
                    int numFastaSeqs = 0;
                    positions = m->setFilePosFasta(fastafileNames[s], numFastaSeqs); 
                    if (numFastaSeqs < processors) { processors = numFastaSeqs; }
                }
#endif
            }
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			//vector<unsigned long long> positions = m->divideFile(fastafileNames[s], processors);
			
			for (int i = 0; i < (positions.size()-1); i++) {
				lines.push_back(new linePair(positions[i], positions[(i+1)]));
			}	
			
				if(processors == 1){
					int numFastaSeqs = driverRunFilter(filter, filteredFasta, fastafileNames[s], lines[0]);
					numSeqs += numFastaSeqs;
				}else{
					int numFastaSeqs = createProcessesRunFilter(filter, fastafileNames[s], filteredFasta); 
					numSeqs += numFastaSeqs;
				}
				
				if (m->control_pressed) {  return 1; }
		#else
            if(processors == 1){
                lines.push_back(new linePair(0, 1000));
				int numFastaSeqs = driverRunFilter(filter, filteredFasta, fastafileNames[s], lines[0]);
				numSeqs += numFastaSeqs;
            }else {
                int numFastaSeqs = positions.size()-1;
                //positions = m->setFilePosFasta(fastafileNames[s], numFastaSeqs); 
                
                //figure out how many sequences you have to process
                int numSeqsPerProcessor = numFastaSeqs / processors;
                for (int i = 0; i < processors; i++) {
                    int startIndex =  i * numSeqsPerProcessor;
                    if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                    lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
                }
                
                numFastaSeqs = createProcessesRunFilter(filter, fastafileNames[s], filteredFasta); 
                numSeqs += numFastaSeqs;
            }

				if (m->control_pressed) {  return 1; }
		#endif

			outputNames.push_back(filteredFasta); outputTypes["fasta"].push_back(filteredFasta);
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "filterSequences");
		exit(1);
	}
}
/**************************************************************************************/
int FilterSeqsCommand::driverRunFilter(string F, string outputFilename, string inputFilename, linePair* filePos) {	
	try {
		ofstream out;
		m->openOutputFile(outputFilename, out);
		
		ifstream in;
		m->openInputFile(inputFilename, in);
				
		in.seekg(filePos->start);
        
        //adjust start if null strings
        if (filePos->start == 0) {  m->zapGremlins(in); m->gobble(in);  }

		bool done = false;
		int count = 0;
	
		while (!done) {
				
				if (m->control_pressed) { in.close(); out.close(); return 0; }
				
				Sequence seq(in); m->gobble(in);
				if (seq.getName() != "") {
					string align = seq.getAligned();
					string filterSeq = "";
					
					for(int j=0;j<alignmentLength;j++){
						if(filter[j] == '1'){
							filterSeq += align[j];
						}
					}
					
					out << '>' << seq.getName() << endl << filterSeq << endl;
				count++;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = in.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (in.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen(toString(count)+"\n"); 	}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen(toString(count)+"\n"); 		}
		
		
		out.close();
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverRunFilter");
		exit(1);
	}
}
/**************************************************************************************************/

int FilterSeqsCommand::createProcessesRunFilter(string F, string filename, string filteredFastaName) {
	try {
        
        int process = 1;
		int num = 0;
		processIDS.clear();
        bool recalc = false;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				string filteredFasta = filename + m->mothurGetpid(process) + ".temp";
				num = driverRunFilter(F, filteredFasta, filename, lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename +  m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
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
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(filename + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(filename + (toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(filename + (toString(processIDS[i]) + ".temp"));m->mothurRemove(filename + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(new linePair(positions[i], positions[(i+1)]));  }
            
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
                    string filteredFasta = filename + m->mothurGetpid(process) + ".temp";
                    num = driverRunFilter(F, filteredFasta, filename, lines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = filename +  m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    out << num << endl;
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }

        }
		
        num = driverRunFilter(F, filteredFastaName, filename, lines[0]);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}	
					
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
            
            m->appendFiles((filename + toString(processIDS[i]) + ".temp"), filteredFastaName);
            m->mothurRemove((filename + toString(processIDS[i]) + ".temp"));
		}
               
#else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the filterData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to F.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<filterRunData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++){
			
            string extension = "";
			if (i != 0) { extension = toString(i) + ".temp"; }
            
			filterRunData* tempFilter = new filterRunData(filter, filename, (filteredFastaName + extension), m, lines[i]->start, lines[i]->end, alignmentLength, i);
			pDataArray.push_back(tempFilter);
			processIDS.push_back(i);
            
			hThreadArray[i] = CreateThread(NULL, 0, MyRunFilterThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
        
        num = driverRunFilter(F, (filteredFastaName + toString(processors-1) + ".temp"), filename, lines[processors-1]);
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        
        for (int i = 1; i < processors; i++) {
            m->appendFiles((filteredFastaName + toString(i) + ".temp"), filteredFastaName);
            m->mothurRemove((filteredFastaName + toString(i) + ".temp"));
		}
#endif	
        
        return num;
        
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createProcessesRunFilter");
		exit(1);
	}
}
/**************************************************************************************/
string FilterSeqsCommand::createFilter() {	
	try {
		string filterString = "";			
		Filters F;
		
		if (soft != 0)			{  F.setSoft(soft);		}
		if (trump != '*')		{  F.setTrump(trump);	}
		
		F.setLength(alignmentLength);
		
		if(trump != '*' || m->isTrue(vertical) || soft != 0){
			F.initialize();
		}
		
		if(hard.compare("") != 0)	{	F.doHard(hard);		}
		else						{	F.setFilter(string(alignmentLength, '1'));	}
		
		numSeqs = 0;
		if(trump != '*' || m->isTrue(vertical) || soft != 0){
			for (int s = 0; s < fastafileNames.size(); s++) {
			
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			
				
                vector<unsigned long long> positions;
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				positions = m->divideFile(fastafileNames[s], processors);
				for (int i = 0; i < (positions.size()-1); i++) {
					lines.push_back(new linePair(positions[i], positions[(i+1)]));
				}	
				
				if(processors == 1){
					int numFastaSeqs = driverCreateFilter(F, fastafileNames[s], lines[0]);
					numSeqs += numFastaSeqs;
				}else{
					int numFastaSeqs = createProcessesCreateFilter(F, fastafileNames[s]); 
					numSeqs += numFastaSeqs;
				}
		#else
                if(processors == 1){
                    lines.push_back(new linePair(0, 1000));
                    int numFastaSeqs = driverCreateFilter(F, fastafileNames[s], lines[0]);
                    numSeqs += numFastaSeqs;
				}else {
                    int numFastaSeqs = 0;
                    positions = m->setFilePosFasta(fastafileNames[s], numFastaSeqs); 
                    if (numFastaSeqs < processors) { processors = numFastaSeqs; }
                    
                    //figure out how many sequences you have to process
                    int numSeqsPerProcessor = numFastaSeqs / processors;
                    for (int i = 0; i < processors; i++) {
                        int startIndex =  i * numSeqsPerProcessor;
                        if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                        lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
                    }
                    
                    numFastaSeqs = createProcessesCreateFilter(F, fastafileNames[s]); 
					numSeqs += numFastaSeqs;
                }
		#endif
                //save the file positions so we can reuse them in the runFilter function
                if (!recalced) {  savedPositions[s] = positions; }
                
				if (m->control_pressed) {  return filterString; }
			
			}
		}

		F.setNumSeqs(numSeqs);
		if(m->isTrue(vertical) == 1)	{	F.doVertical();	}
		if(soft != 0)				{	F.doSoft();		}
		filterString = F.getFilter();
        
		return filterString;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createFilter");
		exit(1);
	}
}
/**************************************************************************************/
int FilterSeqsCommand::driverCreateFilter(Filters& F, string filename, linePair* filePos) {	
	try {
		
		ifstream in;
		m->openInputFile(filename, in);
				
		in.seekg(filePos->start);
        
        //adjust start if null strings
        if (filePos->start == 0) {  m->zapGremlins(in); m->gobble(in);  }

		bool done = false;
		int count = 0;
        bool error = false;
        
		while (!done) {
				
			if (m->control_pressed) { in.close(); return 1; }
					
			Sequence seq(in); m->gobble(in);
			if (seq.getName() != "") {
                    if (m->debug) { m->mothurOutJustToScreen("[DEBUG]: " + seq.getName() + " length = " + toString(seq.getAligned().length()) + '\n'); }
                if (seq.getAligned().length() != alignmentLength) { m->mothurOut("[ERROR]: Sequences are not all the same length, please correct.\n"); error = true; if (!m->debug) { m->control_pressed = true; }else{ m->mothurOutJustToLog("[DEBUG]: " + seq.getName() + " length = " + toString(seq.getAligned().length()) + '\n'); } }
					
					if(trump != '*')			{	F.doTrump(seq);		}
					if(m->isTrue(vertical) || soft != 0)	{	F.getFreqs(seq);	}
					cout.flush();
					count++;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = in.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (in.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen(toString(count)+"\n"); 		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen(toString(count)+"\n"); 	}
		in.close();
		
        if (error) { m->control_pressed = true; }
        
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverCreateFilter");
		exit(1);
	}
}
/**************************************************************************************************/

int FilterSeqsCommand::createProcessesCreateFilter(Filters& F, string filename) {
	try {
        int process = 1;
		int num = 0;
		processIDS.clear();
        bool recalc = false;

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				//reset child's filter counts to 0;
				F.a.clear(); F.a.resize(alignmentLength, 0);
				F.t.clear(); F.t.resize(alignmentLength, 0);
				F.g.clear(); F.g.resize(alignmentLength, 0);
				F.c.clear(); F.c.resize(alignmentLength, 0);
				F.gap.clear(); F.gap.resize(alignmentLength, 0);
				
				num = driverCreateFilter(F, filename, lines[process]);
				
				//write out filter counts to file
				filename += m->mothurGetpid(process) + "filterValues.temp";
				ofstream out;
				m->openOutputFile(filename, out);
				
				out << num << endl;
				out << F.getFilter() << endl;
				for (int k = 0; k < alignmentLength; k++) {		out << F.a[k] << '\t'; }  out << endl;
				for (int k = 0; k < alignmentLength; k++) {		out << F.t[k] << '\t'; }  out << endl;
				for (int k = 0; k < alignmentLength; k++) {		out << F.g[k] << '\t'; }  out << endl;
				for (int k = 0; k < alignmentLength; k++) {		out << F.c[k] << '\t'; }  out << endl;
				for (int k = 0; k < alignmentLength; k++) {		out << F.gap[k] << '\t'; }  out << endl;

				//cout << F.getFilter() << endl;
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
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(filename + (toString(processIDS[i]) + "filterValues.temp"));
                }
                recalc = true;
                break;

			}
		}
		
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(filename + (toString(processIDS[i]) + "filterValues.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(new linePair(positions[i], positions[(i+1)]));  }
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            recalced = true;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    //reset child's filter counts to 0;
                    F.a.clear(); F.a.resize(alignmentLength, 0);
                    F.t.clear(); F.t.resize(alignmentLength, 0);
                    F.g.clear(); F.g.resize(alignmentLength, 0);
                    F.c.clear(); F.c.resize(alignmentLength, 0);
                    F.gap.clear(); F.gap.resize(alignmentLength, 0);
                    
                    num = driverCreateFilter(F, filename, lines[process]);
                    
                    //write out filter counts to file
                    filename += m->mothurGetpid(process) + "filterValues.temp";
                    ofstream out;
                    m->openOutputFile(filename, out);
                    
                    out << num << endl;
                    out << F.getFilter() << endl;
                    for (int k = 0; k < alignmentLength; k++) {		out << F.a[k] << '\t'; }  out << endl;
                    for (int k = 0; k < alignmentLength; k++) {		out << F.t[k] << '\t'; }  out << endl;
                    for (int k = 0; k < alignmentLength; k++) {		out << F.g[k] << '\t'; }  out << endl;
                    for (int k = 0; k < alignmentLength; k++) {		out << F.c[k] << '\t'; }  out << endl;
                    for (int k = 0; k < alignmentLength; k++) {		out << F.gap[k] << '\t'; }  out << endl;
                    
                    //cout << F.getFilter() << endl;
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

        
		//parent do your part
		num = driverCreateFilter(F, filename, lines[0]);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<(processors-1);i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//parent reads in and combines Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = filename + toString(processIDS[i]) + "filterValues.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			int temp, tempNum;
			string tempFilterString;

			in >> tempNum; m->gobble(in); num += tempNum;

			in >> tempFilterString;
			F.mergeFilter(tempFilterString);

			for (int k = 0; k < alignmentLength; k++) {		in >> temp; F.a[k] += temp; }		m->gobble(in);
			for (int k = 0; k < alignmentLength; k++) {		in >> temp; F.t[k] += temp; }		m->gobble(in);
			for (int k = 0; k < alignmentLength; k++) {		in >> temp; F.g[k] += temp; }		m->gobble(in);
			for (int k = 0; k < alignmentLength; k++) {		in >> temp; F.c[k] += temp; }		m->gobble(in);
			for (int k = 0; k < alignmentLength; k++) {		in >> temp; F.gap[k] += temp; }	m->gobble(in);
				
			in.close();
			m->mothurRemove(tempFilename);
		}
		
		
#else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the filterData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to F.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<filterData*> pDataArray; 
		DWORD   dwThreadIdArray[processors];
		HANDLE  hThreadArray[processors]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors; i++ ){
			
			filterData* tempFilter = new filterData(filename, m, lines[i]->start, lines[i]->end, alignmentLength, trump, vertical, soft, hard, i);
			pDataArray.push_back(tempFilter);
			processIDS.push_back(i);
            
			hThreadArray[i] = CreateThread(NULL, 0, MyCreateFilterThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            F.mergeFilter(pDataArray[i]->F.getFilter());
            
			for (int k = 0; k < alignmentLength; k++) {	 F.a[k] += pDataArray[i]->F.a[k];       }
			for (int k = 0; k < alignmentLength; k++) {	 F.t[k] += pDataArray[i]->F.t[k];       }
			for (int k = 0; k < alignmentLength; k++) {	 F.g[k] += pDataArray[i]->F.g[k];       }
			for (int k = 0; k < alignmentLength; k++) {	 F.c[k] += pDataArray[i]->F.c[k];       }
			for (int k = 0; k < alignmentLength; k++) {	 F.gap[k] += pDataArray[i]->F.gap[k];   }

			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
#endif	
        return num;
        
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createProcessesCreateFilter");
		exit(1);
	}
}
/**************************************************************************************/
