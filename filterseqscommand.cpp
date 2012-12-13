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
		abort = false; calledHelp = false;   
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
						int ableToOpen = m->openInputFile(fastafileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (ableToOpen == 1) {
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
						if (ableToOpen == 1) {
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
						
						if (ableToOpen == 1) { 
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
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output the filter
		#endif
		
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
		
		#ifdef USE_MPI
			}
		#endif
		
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
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor, num; 
				int tag = 2001;
				vector<unsigned long long>MPIPos;
						
				MPI_Status status; 
				MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
				MPI_File outMPI;
				MPI_File inMPI;
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				int inMode=MPI_MODE_RDONLY; 
				
				char outFilename[1024];
				strcpy(outFilename, filteredFasta.c_str());
			
				char inFileName[1024];
				strcpy(inFileName, fastafileNames[s].c_str());
				
				MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPI);

				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  return 0;  }

				if (pid == 0) { //you are the root process 
					
					MPIPos = m->setFilePosFasta(fastafileNames[s], num); //fills MPIPos, returns numSeqs
					numSeqs += num;
					
					//send file positions to all processes
					for(int i = 1; i < processors; i++) { 
						MPI_Send(&num, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&MPIPos[0], (num+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
					}
					
					//figure out how many sequences you have to do
					numSeqsPerProcessor = num / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = num - pid * numSeqsPerProcessor; 	}
					
				
					//do your part
					driverMPIRun(startIndex, numSeqsPerProcessor, inMPI, outMPI, MPIPos);
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  return 0;  }
					
					//wait on chidren
					for(int i = 1; i < processors; i++) { 
						char buf[5];
						MPI_Recv(buf, 5, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
					}
					
				}else { //you are a child process
					MPI_Recv(&num, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(num+1);
					numSeqs += num;
					MPI_Recv(&MPIPos[0], (num+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = num / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = num - pid * numSeqsPerProcessor; 	}
					
					
					//align your part
					driverMPIRun(startIndex, numSeqsPerProcessor, inMPI, outMPI, MPIPos);		
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  return 0;  }
					
					char buf[5];
					strcpy(buf, "done"); 
					
					//tell parent you are done.
					MPI_Send(buf, 5, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
				}
				
				MPI_File_close(&outMPI);
				MPI_File_close(&inMPI);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
#else
            
            vector<unsigned long long> positions;
            if (savedPositions.size() != 0) { positions = savedPositions[s]; }
            else {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				positions = m->divideFile(fastafileNames[s], processors);
#else
                if(processors != 1){
                    int numFastaSeqs = 0;
                    positions = m->setFilePosFasta(fastafileNames[s], numFastaSeqs); 
                    if (positions.size() < processors) { processors = positions.size(); }
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
#ifdef USE_MPI
/**************************************************************************************/
int FilterSeqsCommand::driverMPIRun(int start, int num, MPI_File& inMPI, MPI_File& outMPI, vector<unsigned long long>& MPIPos) {	
	try {
		string outputString = "";
		int count = 0;
		MPI_Status status; 
		
		for(int i=0;i<num;i++){
		
			if (m->control_pressed) { return 0; }
		
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];
			char* buf4 = new char[length];
			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);
			delete buf4;
	
			Sequence seq(iss);  m->gobble(iss);
			
			if (seq.getName() != "") {
				string align = seq.getAligned();
				string filterSeq = "";
					
				for(int j=0;j<alignmentLength;j++){
					if(filter[j] == '1'){
						filterSeq += align[j];
					}
				}
				
				count++;
				outputString += ">" + seq.getName() + "\n" + filterSeq + "\n";
				
				if(count % 10 == 0){ //output to file 
					//send results to parent
					int length = outputString.length();
					char* buf = new char[length];
					memcpy(buf, outputString.c_str(), length);
				
					MPI_File_write_shared(outMPI, buf, length, MPI_CHAR, &status);
					outputString = "";
					delete buf;
				}

			}
			
			if((i+1) % 100 == 0){	cout << (i+1) << endl;	 m->mothurOutJustToLog(toString(i+1) + "\n");	}
		}
		
		if(outputString != ""){ //output to file 
			//send results to parent
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write_shared(outMPI, buf, length, MPI_CHAR, &status);
			outputString = "";
			delete buf;
		}
		
		if((num) % 100 != 0){	cout << (num) << endl;	 m->mothurOutJustToLog(toString(num) + "\n");	}
			
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverRunFilter");
		exit(1);
	}
}
#endif
/**************************************************************************************/
int FilterSeqsCommand::driverRunFilter(string F, string outputFilename, string inputFilename, linePair* filePos) {	
	try {
		ofstream out;
		m->openOutputFile(outputFilename, out);
		
		ifstream in;
		m->openInputFile(inputFilename, in);
				
		in.seekg(filePos->start);

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
			if((count) % 100 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		
		
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
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				string filteredFasta = filename + toString(getpid()) + ".temp";
				num = driverRunFilter(F, filteredFasta, filename, lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename +  toString(getpid()) + ".num.temp";
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
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor, num; 
				int tag = 2001;
				vector<unsigned long long> MPIPos;
				
				MPI_Status status; 
				MPI_File inMPI; 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
							
				//char* tempFileName = new char(fastafileNames[s].length());
				//tempFileName = &(fastafileNames[s][0]);
				
				char tempFileName[1024];
				strcpy(tempFileName, fastafileNames[s].c_str());
		
				MPI_File_open(MPI_COMM_WORLD, tempFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  return 0;  }
				
				if (pid == 0) { //you are the root process
						MPIPos = m->setFilePosFasta(fastafileNames[s], num); //fills MPIPos, returns numSeqs
						numSeqs += num;
						
						//send file positions to all processes
						for(int i = 1; i < processors; i++) { 
							MPI_Send(&num, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
							MPI_Send(&MPIPos[0], (num+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
						}
								
						//figure out how many sequences you have to do
						numSeqsPerProcessor = num / processors;
						int startIndex =  pid * numSeqsPerProcessor;
						if(pid == (processors - 1)){	numSeqsPerProcessor = num - pid * numSeqsPerProcessor; 	}
						
				
						//do your part
						MPICreateFilter(startIndex, numSeqsPerProcessor, F, inMPI, MPIPos);
						
						if (m->control_pressed) {  MPI_File_close(&inMPI);  return 0;  }
												
				}else { //i am the child process
					MPI_Recv(&num, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(num+1);
					numSeqs += num;
					MPI_Recv(&MPIPos[0], (num+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = num / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = num - pid * numSeqsPerProcessor; 	}
					
					
					//do your part
					MPICreateFilter(startIndex, numSeqsPerProcessor, F, inMPI,  MPIPos);
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  return 0;  }
				}
				
				MPI_File_close(&inMPI);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
#else
				
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
                    if (positions.size() < processors) { processors = positions.size(); }
                    
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
                savedPositions[s] = positions;
                
				if (m->control_pressed) {  return filterString; }
#endif
			
			}
		}


#ifdef USE_MPI	
		int pid;
		int Atag = 1; int Ttag = 2; int Ctag = 3; int Gtag = 4; int Gaptag = 5;
		MPI_Status status;
		
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if(trump != '*' || m->isTrue(vertical) || soft != 0){
			
			if (pid == 0) { //only one process should output the filter
			
				vector<int> temp; temp.resize(alignmentLength+1);
								
				//get the frequencies from the child processes
				for(int i = 1; i < processors; i++) { 
				
					for (int j = 0; j < 5; j++) {
					
						MPI_Recv(&temp[0], (alignmentLength+1), MPI_INT, i, 2001, MPI_COMM_WORLD, &status); 
						int receiveTag = temp[temp.size()-1];  //child process added a int to the end to indicate what letter count this is for
						
						if (receiveTag == Atag) { //you are recieveing the A frequencies
							for (int k = 0; k < alignmentLength; k++) {		F.a[k] += temp[k];	}
						}else if (receiveTag == Ttag) { //you are recieveing the T frequencies
							for (int k = 0; k < alignmentLength; k++) {		F.t[k] += temp[k];	}
						}else if (receiveTag == Ctag) { //you are recieveing the C frequencies
							for (int k = 0; k < alignmentLength; k++) {		F.c[k] += temp[k];	}
						}else if (receiveTag == Gtag) { //you are recieveing the G frequencies
							for (int k = 0; k < alignmentLength; k++) {		F.g[k] += temp[k];	}
						}else if (receiveTag == Gaptag) { //you are recieveing the gap frequencies
							for (int k = 0; k < alignmentLength; k++) {		F.gap[k] += temp[k];	}
						}
					}
				} 
			}else{
			
				//send my fequency counts
				F.a.push_back(Atag);
				int ierr = MPI_Send(&(F.a[0]), (alignmentLength+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
				F.t.push_back(Ttag);
				ierr = MPI_Send (&(F.t[0]), (alignmentLength+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
				F.c.push_back(Ctag);
				ierr = MPI_Send(&(F.c[0]), (alignmentLength+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
				F.g.push_back(Gtag);
				ierr = MPI_Send(&(F.g[0]), (alignmentLength+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
				F.gap.push_back(Gaptag);
				ierr = MPI_Send(&(F.gap[0]), (alignmentLength+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
			}
			
		}
		
		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
		
		if (pid == 0) { //only one process should output the filter
#endif

		F.setNumSeqs(numSeqs);
		if(m->isTrue(vertical) == 1)	{	F.doVertical();	}
		if(soft != 0)				{	F.doSoft();		}
		filterString = F.getFilter();
		
#ifdef USE_MPI
		//send filter string to kids
		//for(int i = 1; i < processors; i++) { 
		//	MPI_Send(&filterString[0], alignmentLength, MPI_CHAR, i, 2001, MPI_COMM_WORLD);
		//}
		MPI_Bcast(&filterString[0], alignmentLength, MPI_CHAR, 0, MPI_COMM_WORLD);
	}else{
		//recieve filterString
		char* tempBuf = new char[alignmentLength];
		//MPI_Recv(&tempBuf[0], alignmentLength, MPI_CHAR, 0, 2001, MPI_COMM_WORLD, &status);
		MPI_Bcast(tempBuf, alignmentLength, MPI_CHAR, 0, MPI_COMM_WORLD);
		
		filterString = tempBuf;
		if (filterString.length() > alignmentLength) { filterString = filterString.substr(0, alignmentLength);  }
		delete tempBuf;	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
#endif
            
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

		bool done = false;
		int count = 0;
	
		while (!done) {
				
			if (m->control_pressed) { in.close(); return 1; }
					
			Sequence seq(in); m->gobble(in);
			if (seq.getName() != "") {
					if (seq.getAligned().length() != alignmentLength) { m->mothurOut("Sequences are not all the same length, please correct."); m->mothurOutEndLine(); m->control_pressed = true;  }
					
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
			if((count) % 100 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverCreateFilter");
		exit(1);
	}
}
#ifdef USE_MPI
/**************************************************************************************/
int FilterSeqsCommand::MPICreateFilter(int start, int num, Filters& F, MPI_File& inMPI, vector<unsigned long long>& MPIPos) {	
	try {
		
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		for(int i=0;i<num;i++){
			
			if (m->control_pressed) { return 0; }
			
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];
	
			char* buf4 = new char[length];
			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);
			delete buf4;

			Sequence seq(iss);  

			if (seq.getAligned().length() != alignmentLength) {  cout << "Alignment length is " << alignmentLength << " and sequence " << seq.getName() << " has length " << seq.getAligned().length() << ", please correct." << endl; exit(1);  }
			
			if(trump != '*'){	F.doTrump(seq);	}
			if(m->isTrue(vertical) || soft != 0){	F.getFreqs(seq);	}
			cout.flush();
						
			//report progress
			if((i+1) % 100 == 0){	cout << (i+1) << endl;	 m->mothurOutJustToLog(toString(i+1) + "\n");	}
		}
		
		//report progress
		if((num) % 100 != 0){	cout << num << endl; m->mothurOutJustToLog(toString(num) + "\n"); 	}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "MPICreateFilter");
		exit(1);
	}
}
#endif
/**************************************************************************************************/

int FilterSeqsCommand::createProcessesCreateFilter(Filters& F, string filename) {
	try {
        int process = 1;
		int num = 0;
		processIDS.clear();

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
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
				filename += toString(getpid()) + "filterValues.temp";
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
