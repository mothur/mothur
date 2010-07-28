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

/**************************************************************************************/

FilterSeqsCommand::FilterSeqsCommand(string option)  {
	try {
		abort = false;
		filterFileName = "";
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "trump", "soft", "hard", "vertical", "outputdir","inputdir", "processors"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("filter.seqs");
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
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("hard");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["hard"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fasta = validParameter.validFile(parameters, "fasta", false);
			if (fasta == "not found") { m->mothurOut("fasta is a required parameter for the filter.seqs command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				splitAtDash(fasta, fastafileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastafileNames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(fastafileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastafileNames[i] = inputDir + fastafileNames[i];		}
					}

					ifstream in;
					int ableToOpen = openInputFile(fastafileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + getSimpleName(fastafileNames[i]);
							m->mothurOut("Unable to open " + fastafileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ableToOpen = openInputFile(tryPath, in, "noerror");
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
						string simpleName = getSimpleName(fastafileNames[i]);
						filterFileName += simpleName.substr(0, simpleName.find_first_of('.'));
					}
					in.close();
				}
				
				//make sure there is at least one valid file left
				if (fastafileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (!abort) {
				//if the user changes the output directory command factory will send this info to us in the output parameter 
				outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
					outputDir = "";	
					outputDir += hasPath(fastafileNames[0]); //if user entered a file with a path then preserve it	
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
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
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

//**********************************************************************************************************************

void FilterSeqsCommand::help(){
	try {
				
		m->mothurOut("The filter.seqs command reads a file containing sequences and creates a .filter and .filter.fasta file.\n");
		m->mothurOut("The filter.seqs command parameters are fasta, trump, soft, hard and vertical. \n");
		m->mothurOut("The fasta parameter is required. You may enter several fasta files to build the filter from and filter, by separating their names with -'s.\n");
		m->mothurOut("For example: fasta=abrecovery.fasta-amazon.fasta \n");
		m->mothurOut("The trump parameter .... The default is ...\n");
		m->mothurOut("The soft parameter .... The default is ....\n");
		m->mothurOut("The hard parameter allows you to enter a file containing the filter you want to use.\n");
		m->mothurOut("The vertical parameter removes columns where all sequences contain a gap character. The default is T.\n");
		m->mothurOut("The filter.seqs command should be in the following format: \n");
		m->mothurOut("filter.seqs(fasta=yourFastaFile, trump=yourTrump) \n");
		m->mothurOut("Example filter.seqs(fasta=abrecovery.fasta, trump=.).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "help");
		exit(1);
	}
}

/**************************************************************************************/

int FilterSeqsCommand::execute() {	
	try {
	
		if (abort == true) { return 0; }
		
		ifstream inFASTA;
		openInputFile(fastafileNames[0], inFASTA);
		
		Sequence testSeq(inFASTA);
		alignmentLength = testSeq.getAlignLength();
		inFASTA.close();
		
		////////////create filter/////////////////
		m->mothurOut("Creating Filter... "); m->mothurOutEndLine();
		
		filter = createFilter();
		
		m->mothurOutEndLine();  m->mothurOutEndLine();
		
		if (m->control_pressed) { return 0; }
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output the filter
		#endif
		
		ofstream outFilter;
		
		string filterFile = outputDir + filterFileName + ".filter";
		openOutputFile(filterFile, outFilter);
		outFilter << filter << endl;
		outFilter.close();
		outputNames.push_back(filterFile);
		
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
		
		if (m->control_pressed) {  for(int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); }  return 0; }

		
		m->mothurOutEndLine();
		m->mothurOut("Length of filtered alignment: " + toString(filteredLength)); m->mothurOutEndLine();
		m->mothurOut("Number of columns removed: " + toString((alignmentLength-filteredLength))); m->mothurOutEndLine();
		m->mothurOut("Length of the original alignment: " + toString(alignmentLength)); m->mothurOutEndLine();
		m->mothurOut("Number of sequences used to construct filter: " + toString(numSeqs)); m->mothurOutEndLine();
		
		
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
				
				string filteredFasta = outputDir + getRootName(getSimpleName(fastafileNames[s])) + "filter.fasta";
#ifdef USE_MPI	
				int pid, start, end, numSeqsPerProcessor, num; 
				int tag = 2001;
				vector<unsigned long int>MPIPos;
						
				MPI_Status status; 
				MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
				MPI_File outMPI;
				MPI_File tempMPI;
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
					
					MPIPos = setFilePosFasta(fastafileNames[s], num); //fills MPIPos, returns numSeqs
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
						char buf[4];
						MPI_Recv(buf, 4, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
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
					
					char buf[4];
					strcpy(buf, "done"); 
					
					//tell parent you are done.
					MPI_Send(buf, 4, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
				}
				
				MPI_File_close(&outMPI);
				MPI_File_close(&inMPI);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
#else
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					ifstream inFASTA;
					int numFastaSeqs;
					openInputFile(fastafileNames[s], inFASTA);
					getNumSeqs(inFASTA, numFastaSeqs);
					inFASTA.close();
					
					lines.push_back(new linePair(0, numFastaSeqs));
					
					numSeqs += numFastaSeqs;
					
					driverRunFilter(filter, filteredFasta, fastafileNames[s], lines[0]);
				}else{
					setLines(fastafileNames[s]);
					createProcessesRunFilter(filter, fastafileNames[s]); 
				
					rename((fastafileNames[s] + toString(processIDS[0]) + ".temp").c_str(), filteredFasta.c_str());
				
					//append fasta files
					for(int i=1;i<processors;i++){
						appendFiles((fastafileNames[s] + toString(processIDS[i]) + ".temp"), filteredFasta);
						remove((fastafileNames[s] + toString(processIDS[i]) + ".temp").c_str());
					}
				}
				
				if (m->control_pressed) {  return 1; }
		#else
				ifstream inFASTA;
				int numFastaSeqs;
				openInputFile(fastafileNames[s], inFASTA);
				getNumSeqs(inFASTA, numFastaSeqs);
				inFASTA.close();
					
				lines.push_back(new linePair(0, numFastaSeqs));
				
				numSeqs += numFastaSeqs;
				
				driverRunFilter(filter, filteredFasta, fastafileNames[s], lines[0]);

				if (m->control_pressed) {  return 1; }
		#endif
#endif
			outputNames.push_back(filteredFasta);
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
int FilterSeqsCommand::driverMPIRun(int start, int num, MPI_File& inMPI, MPI_File& outMPI, vector<unsigned long int>& MPIPos) {	
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
	
			Sequence seq(iss);  gobble(iss);
			
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
int FilterSeqsCommand::driverRunFilter(string F, string outputFilename, string inputFilename, linePair* line) {	
	try {
		ofstream out;
		openOutputFile(outputFilename, out);
		
		ifstream in;
		openInputFile(inputFilename, in);
				
		in.seekg(line->start);
		
		for(int i=0;i<line->num;i++){
				
				if (m->control_pressed) { in.close(); out.close(); return 0; }
				
				Sequence seq(in);
				if (seq.getName() != "") {
					string align = seq.getAligned();
					string filterSeq = "";
					
					for(int j=0;j<alignmentLength;j++){
						if(filter[j] == '1'){
							filterSeq += align[j];
						}
					}
					
					out << '>' << seq.getName() << endl << filterSeq << endl;
				}
				gobble(in);
		}
		out.close();
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverRunFilter");
		exit(1);
	}
}
/**************************************************************************************************/

int FilterSeqsCommand::createProcessesRunFilter(string F, string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		int exitCommand = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				string filteredFasta = filename + toString(getpid()) + ".temp";
				driverRunFilter(F, filteredFasta, filename, lines[process]);
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		return exitCommand;
#endif		
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
		
		if(trump != '*' || isTrue(vertical) || soft != 0){
			F.initialize();
		}
		
		if(hard.compare("") != 0)	{	F.doHard(hard);		}
		else						{	F.setFilter(string(alignmentLength, '1'));	}
		
		numSeqs = 0;
		if(trump != '*' || isTrue(vertical) || soft != 0){
			for (int s = 0; s < fastafileNames.size(); s++) {
			
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor, num; 
				int tag = 2001;
				vector<unsigned long int> MPIPos;
				
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
						MPIPos = setFilePosFasta(fastafileNames[s], num); //fills MPIPos, returns numSeqs
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
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					ifstream inFASTA;
					int numFastaSeqs;
					openInputFile(fastafileNames[s], inFASTA);
					getNumSeqs(inFASTA, numFastaSeqs);
					inFASTA.close();
					
					numSeqs += numFastaSeqs;
					
					lines.push_back(new linePair(0, numFastaSeqs));
					
					driverCreateFilter(F, fastafileNames[s], lines[0]);
				}else{
					setLines(fastafileNames[s]);					
					createProcessesCreateFilter(F, fastafileNames[s]); 
				}
				
				if (m->control_pressed) {  return filterString; }
		#else
				ifstream inFASTA;
				int numFastaSeqs;
				openInputFile(fastafileNames[s], inFASTA);
				getNumSeqs(inFASTA, numFastaSeqs);
				inFASTA.close();
				
				numSeqs += numFastaSeqs;
				
				lines.push_back(new linePair(0, numFastaSeqs));
				
				driverCreateFilter(F, fastafileNames[s], lines[0]);
				if (m->control_pressed) {  return filterString; }
		#endif
#endif
			
			}
		}


#ifdef USE_MPI	
		int pid;
		int Atag = 1; int Ttag = 2; int Ctag = 3; int Gtag = 4; int Gaptag = 5;
		MPI_Status status;
		
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if(trump != '*' || isTrue(vertical) || soft != 0){
			
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
				
		if(isTrue(vertical) == 1)	{	F.doVertical();	}
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
int FilterSeqsCommand::driverCreateFilter(Filters& F, string filename, linePair* line) {	
	try {
		
		ifstream in;
		openInputFile(filename, in);
				
		in.seekg(line->start);
		
		for(int i=0;i<line->num;i++){
				
			if (m->control_pressed) { in.close(); return 1; }
					
			Sequence seq(in);
			if (seq.getName() != "") {
					if (seq.getAligned().length() != alignmentLength) { m->mothurOut("Sequences are not all the same length, please correct."); m->mothurOutEndLine(); m->control_pressed = true;  }
					
					if(trump != '*'){	F.doTrump(seq);	}
					if(isTrue(vertical) || soft != 0){	F.getFreqs(seq);	}
					cout.flush();
			}
			
			//report progress
			if((i+1) % 100 == 0){	m->mothurOut(toString(i+1)); m->mothurOutEndLine();		}
		}
		
		//report progress
		if((line->num) % 100 != 0){	m->mothurOut(toString(line->num)); m->mothurOutEndLine();		}
		
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverCreateFilter");
		exit(1);
	}
}
#ifdef USE_MPI
/**************************************************************************************/
int FilterSeqsCommand::MPICreateFilter(int start, int num, Filters& F, MPI_File& inMPI, vector<unsigned long int>& MPIPos) {	
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
			if(isTrue(vertical) || soft != 0){	F.getFreqs(seq);	}
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
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		int exitCommand = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = vfork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driverCreateFilter(F, filename, lines[process]);
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		return exitCommand;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createProcessesCreateFilter");
		exit(1);
	}
}
/**************************************************************************************************/

int FilterSeqsCommand::setLines(string filename) {
	try {
		
		vector<unsigned long int> positions;
		bufferSizes.clear();
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);
			
		string input;
		while(!inFASTA.eof()){	
			input = getline(inFASTA);

			if (input.length() != 0) {
				if(input[0] == '>'){ unsigned long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
			}
		}
		inFASTA.close();
		
		int numFastaSeqs = positions.size();
	
		FILE * pFile;
		unsigned long int size;
		
		//get num bytes in file
		pFile = fopen (filename.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
	
		numSeqs += numFastaSeqs;
		
		int numSeqsPerProcessor = numFastaSeqs / processors;
		
		for (int i = 0; i < processors; i++) {

			unsigned long int startPos = positions[ i * numSeqsPerProcessor ];
			if(i == processors - 1){
				numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;
				bufferSizes.push_back(size - startPos);
			}else{  
				unsigned long int myEnd = positions[ (i+1) * numSeqsPerProcessor ];
				bufferSizes.push_back(myEnd-startPos);
			}
			lines.push_back(new linePair(startPos, numSeqsPerProcessor));
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "setLines");
		exit(1);
	}
}
/**************************************************************************************/
