/*
 *  chimeracheckcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeracheckcommand.h"

//***************************************************************************************************************

ChimeraCheckCommand::ChimeraCheckCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","processors","increment","template","ksize","svg", "name","outputdir","inputdir" };
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.check");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				it = parameters.find("template");
				//user has given a template file
				if(it != parameters.end()){ 
					string path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["template"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.check command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(fastaFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
					}
	
					int ableToOpen;
					ifstream in;
					
					#ifdef USE_MPI	
						int pid;
						MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
						if (pid == 0) {
					#endif

					ableToOpen = openInputFile(fastaFileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ableToOpen = openInputFile(tryPath, in, "noerror");
							fastaFileNames[i] = tryPath;
						}
					}
					in.close();
					
					#ifdef USE_MPI	
							for (int j = 1; j < processors; j++) {
								MPI_Send(&ableToOpen, 1, MPI_INT, j, 2001, MPI_COMM_WORLD); 
							}
						}else{
							MPI_Status status;
							MPI_Recv(&ableToOpen, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
						}
						
					#endif

					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + fastaFileNames[i] +". It will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

			templatefile = validParameter.validFile(parameters, "template", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = "";  m->mothurOut("template is a required parameter for the chimera.check command."); m->mothurOutEndLine(); abort = true;  }	
			
			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = ""; }
			else { 
				splitAtDash(namefile, nameFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < nameFileNames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(nameFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	nameFileNames[i] = inputDir + nameFileNames[i];		}
					}
	
					int ableToOpen;
					ifstream in;
					
					#ifdef USE_MPI	
						int pid;
						MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
						if (pid == 0) {
					#endif

					ableToOpen = openInputFile(nameFileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + getSimpleName(nameFileNames[i]);
							m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ableToOpen = openInputFile(tryPath, in, "noerror");
							nameFileNames[i] = tryPath;
						}
					}
					in.close();
					
					#ifdef USE_MPI	
							for (int j = 1; j < processors; j++) {
								MPI_Send(&ableToOpen, 1, MPI_INT, j, 2001, MPI_COMM_WORLD); 
							}
						}else{
							MPI_Status status;
							MPI_Recv(&ableToOpen, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
						}
						
					#endif

					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + nameFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						nameFileNames.erase(nameFileNames.begin()+i);
						i--;
					}
					
				}
				
				//make sure there is at least one valid file left
				if (nameFileNames.size() != 0) {
					if (nameFileNames.size() != fastaFileNames.size()) { 
						 m->mothurOut("Different number of valid name files and fasta files, aborting command."); m->mothurOutEndLine(); 
						 abort = true;
					}
				}
			}

			string temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			convert(temp, ksize);
			
			temp = validParameter.validFile(parameters, "svg", false);				if (temp == "not found") { temp = "F"; }
			svg = isTrue(temp);
			if (nameFileNames.size() != 0) { svg = true; }
			
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, increment);			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "ChimeraCheckCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraCheckCommand::help(){
	try {
	
		m->mothurOut("The chimera.check command reads a fastafile and templatefile and outputs potentially chimeric sequences.\n");
		m->mothurOut("This command was created using the algorythms described in CHIMERA_CHECK version 2.7 written by Niels Larsen. \n");
		m->mothurOut("The chimera.check command parameters are fasta, template, processors, ksize, increment, svg and name.\n");
		m->mothurOut("The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required. \n");
		m->mothurOut("You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences, and is required. \n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
		m->mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default is 10.\n");
		m->mothurOut("The ksize parameter allows you to input kmersize, default is 7. \n");
		m->mothurOut("The svg parameter allows you to specify whether or not you would like a svg file outputted for each query sequence, default is False.\n");
		m->mothurOut("The name parameter allows you to enter a file containing names of sequences you would like .svg files for.\n");
		m->mothurOut("You may enter multiple name files by separating their names with dashes. ie. fasta=abrecovery.svg.names-amzon.svg.names \n");
		m->mothurOut("The chimera.check command should be in the following format: \n");
		m->mothurOut("chimera.check(fasta=yourFastaFile, template=yourTemplateFile, processors=yourProcessors, ksize=yourKmerSize) \n");
		m->mothurOut("Example: chimera.check(fasta=AD.fasta, template=core_set_aligned,imputed.fasta, processors=4, ksize=8) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraCheckCommand::~ChimeraCheckCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraCheckCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[i] + " ..." ); m->mothurOutEndLine();
			
			int start = time(NULL);	
			
			string thisNameFile = "";
			if (nameFileNames.size() != 0) { thisNameFile = nameFileNames[i]; }
			
			chimera = new ChimeraCheckRDP(fastaFileNames[i], templatefile, thisNameFile, svg, increment, ksize, outputDir);			

			if (m->control_pressed) { delete chimera;	return 0;	}
			
			string outputFileName = outputDir + getRootName(getSimpleName(fastaFileNames[i]))  + "chimeracheck.chimeras";
			outputNames.push_back(outputFileName);
			
		#ifdef USE_MPI
		
				int pid, end, numSeqsPerProcessor; 
				int tag = 2001;
				vector<unsigned long int> MPIPos;
				
				MPI_Status status; 
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				MPI_Comm_size(MPI_COMM_WORLD, &processors); 

				MPI_File inMPI;
				MPI_File outMPI;
							
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				int inMode=MPI_MODE_RDONLY; 
							
				char outFilename[1024];
				strcpy(outFilename, outputFileName.c_str());
			
				char inFileName[1024];
				strcpy(inFileName, fastaFileNames[i].c_str());

				MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPI);
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} delete chimera; return 0;  }
				
				if (pid == 0) { //you are the root process 
					MPIPos = setFilePosFasta(fastaFileNames[i], numSeqs); //fills MPIPos, returns numSeqs
					
					//send file positions to all processes
					for(int j = 1; j < processors; j++) { 
						MPI_Send(&numSeqs, 1, MPI_INT, j, tag, MPI_COMM_WORLD);
						MPI_Send(&MPIPos[0], (numSeqs+1), MPI_LONG, j, tag, MPI_COMM_WORLD);
					}	
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
					
				
					//align your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, MPIPos);
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}   delete chimera; return 0;  }
					
					//wait on chidren
					for(int j = 1; j < processors; j++) { 
						char buf[4];
						MPI_Recv(buf, 4, MPI_CHAR, j, tag, MPI_COMM_WORLD, &status); 
					}
				}else{ //you are a child process
					MPI_Recv(&numSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(numSeqs+1);
					MPI_Recv(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
					
					//align your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, MPIPos);
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  delete chimera; return 0;  }
					
					//tell parent you are done.
					char buf[4];
					strcpy(buf, "done"); 
					MPI_Send(buf, 4, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
				}
				
				//close files 
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPI);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
		#else
			
			//break up file
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					ifstream inFASTA;
					openInputFile(fastaFileNames[i], inFASTA);
					getNumSeqs(inFASTA, numSeqs);
					inFASTA.close();
					
					lines.push_back(new linePair(0, numSeqs));
					
					driver(lines[0], outputFileName, fastaFileNames[i]);
					
					if (m->control_pressed) { 
						for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} 
						for (int j = 0; j < lines.size(); j++) {  delete lines[j];  }  lines.clear();
						delete chimera;
						return 0;
					}
									
				}else{
					vector<unsigned long int> positions;
					processIDS.resize(0);
					
					ifstream inFASTA;
					openInputFile(fastaFileNames[i], inFASTA);
					
					string input;
					while(!inFASTA.eof()){
						input = getline(inFASTA);
						if (input.length() != 0) {
							if(input[0] == '>'){	unsigned long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
						}
					}
					inFASTA.close();
					
					numSeqs = positions.size();
					
					int numSeqsPerProcessor = numSeqs / processors;
					
					for (int j = 0; j < processors; j++) {
						unsigned long int startPos = positions[ j * numSeqsPerProcessor ];
						if(j == processors - 1){
							numSeqsPerProcessor = numSeqs - j * numSeqsPerProcessor;
						}
						lines.push_back(new linePair(startPos, numSeqsPerProcessor));
					}
					
					
					createProcesses(outputFileName, fastaFileNames[i]); 
				
					rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
						
					//append output files
					for(int j=1;j<processors;j++){
						appendFiles((outputFileName + toString(processIDS[j]) + ".temp"), outputFileName);
						remove((outputFileName + toString(processIDS[j]) + ".temp").c_str());
					}
					
					if (m->control_pressed) { 
						for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} 
						for (int j = 0; j < lines.size(); j++) {  delete lines[j];  }  lines.clear();
						delete chimera;
						return 0;
					}
				}

			#else
				ifstream inFASTA;
				openInputFile(fastaFileNames[i], inFASTA);
				getNumSeqs(inFASTA, numSeqs);
				inFASTA.close();
				lines.push_back(new linePair(0, numSeqs));
				
				driver(lines[0], outputFileName, fastaFileNames[i]);
				
				if (m->control_pressed) { 
						for (int j = 0; j < lines.size(); j++) {  delete lines[j];  }  lines.clear();
						for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} 
						delete chimera;
						return 0;
				}
			#endif
		#endif		
			delete chimera;
			for (int j = 0; j < lines.size(); j++) {  delete lines[j];  }  lines.clear();
			
			m->mothurOutEndLine(); m->mothurOut("This method does not determine if a sequence is chimeric, but allows you to make that determination based on the IS values."); m->mothurOutEndLine(); 
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine(); m->mothurOutEndLine();

		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
	
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraCheckCommand::driver(linePair* line, string outputFName, string filename){
	try {
		ofstream out;
		openOutputFile(outputFName, out);
		
		ofstream out2;
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);

		inFASTA.seekg(line->start);
		
		for(int i=0;i<line->numSeqs;i++){
		
			if (m->control_pressed) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				//find chimeras
				chimera->getChimeras(candidateSeq);
				
				if (m->control_pressed) {	delete candidateSeq; return 1;	}
	
				//print results
				chimera->print(out, out2);
			}
			delete candidateSeq;
			
			//report progress
			if((i+1) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(i+1)); m->mothurOutEndLine();		}
		}
		//report progress
		if((line->numSeqs) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(line->numSeqs)); m->mothurOutEndLine();		}
		
		out.close();
		inFASTA.close();
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ChimeraCheckCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& outMPI, vector<unsigned long int>& MPIPos){
	try {
		MPI_File outAccMPI;
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

			Sequence* candidateSeq = new Sequence(iss);  gobble(iss);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				//find chimeras
				chimera->getChimeras(candidateSeq);
					
				//print results
				chimera->print(outMPI, outAccMPI);
			}
			delete candidateSeq;
			
			//report progress
			if((i+1) % 100 == 0){  cout << "Processing sequence: " << (i+1) << endl;	m->mothurOutJustToLog("Processing sequence: " + toString(i+1) + "\n");		}
		}
		//report progress
		if(num % 100 != 0){		cout << "Processing sequence: " << num << endl;	m->mothurOutJustToLog("Processing sequence: " + toString(num) + "\n"); 	}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "driverMPI");
		exit(1);
	}
}
#endif

/**************************************************************************************************/

int ChimeraCheckCommand::createProcesses(string outputFileName, string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		//		processIDS.resize(0);
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driver(lines[process], outputFileName + toString(getpid()) + ".temp", filename);
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		return 0;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/


