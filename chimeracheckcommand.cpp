/*
 *  chimeracheckcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeracheckcommand.h"

//**********************************************************************************************************************
vector<string> ChimeraCheckCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter psvg("svg", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psvg);
		CommandParameter pincrement("increment", "Number", "", "10", "", "", "",false,false); parameters.push_back(pincrement);
		CommandParameter pksize("ksize", "Number", "", "7", "", "", "",false,false); parameters.push_back(pksize);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraCheckCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.check command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command was created using the algorythms described in CHIMERA_CHECK version 2.7 written by Niels Larsen. \n";
		helpString += "The chimera.check command parameters are fasta, reference, processors, ksize, increment, svg and name.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required unless you have a valid current fasta file. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
#ifdef USE_MPI
		helpString += "When using MPI, the processors parameter is set to the number of MPI processes running. \n";
#endif
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default is 10.\n";
		helpString += "The ksize parameter allows you to input kmersize, default is 7. \n";
		helpString += "The svg parameter allows you to specify whether or not you would like a svg file outputted for each query sequence, default is False.\n";
		helpString += "The name parameter allows you to enter a file containing names of sequences you would like .svg files for.\n";
		helpString += "You may enter multiple name files by separating their names with dashes. ie. fasta=abrecovery.svg.names-amzon.svg.names \n";
		helpString += "The chimera.check command should be in the following format: \n";
		helpString += "chimera.check(fasta=yourFastaFile, reference=yourTemplateFile, processors=yourProcessors, ksize=yourKmerSize) \n";
		helpString += "Example: chimera.check(fasta=AD.fasta, reference=core_set_aligned,imputed.fasta, processors=4, ksize=8) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
ChimeraCheckCommand::ChimeraCheckCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "ChimeraCheckCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraCheckCommand::ChimeraCheckCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.check");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					string path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
				string filename = m->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				m->splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = m->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
					
					
						if (inputDir != "") {
							string path = m->hasPath(fastaFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
						}
		
						int ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
								m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastaFileNames[i] = tryPath;
							}
						}
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(fastaFileNames[i]);
								m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastaFileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + fastaFileNames[i] +". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

			templatefile = validParameter.validFile(parameters, "reference", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = "";  m->mothurOut("reference is a required parameter for the chimera.check command."); m->mothurOutEndLine(); abort = true;  }	
			
			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = ""; }
			else { 
				m->splitAtDash(namefile, nameFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < nameFileNames.size(); i++) {
					
					bool ignore = false;
					if (nameFileNames[i] == "current") { 
						nameFileNames[i] = m->getNameFile(); 
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + nameFileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
					
						if (inputDir != "") {
							string path = m->hasPath(nameFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	nameFileNames[i] = inputDir + nameFileNames[i];		}
						}
		
						int ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(nameFileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + nameFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
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

			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			convert(temp, ksize);
			
			temp = validParameter.validFile(parameters, "svg", false);				if (temp == "not found") { temp = "F"; }
			svg = m->isTrue(temp);
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
//***************************************************************************************************************

int ChimeraCheckCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[i] + " ..." ); m->mothurOutEndLine();
			
			int start = time(NULL);	
			
			string thisNameFile = "";
			if (nameFileNames.size() != 0) { thisNameFile = nameFileNames[i]; }
			
			chimera = new ChimeraCheckRDP(fastaFileNames[i], templatefile, thisNameFile, svg, increment, ksize, outputDir);			

			if (m->control_pressed) { delete chimera;	return 0;	}
			
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[i]);  }//if user entered a file with a path then preserve it
			string outputFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[i]))  + "chimeracheck.chimeras";
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			
		#ifdef USE_MPI
		
				int pid, numSeqsPerProcessor; 
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
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} outputTypes.clear(); delete chimera; return 0;  }
				
				if (pid == 0) { //you are the root process 
					MPIPos = m->setFilePosFasta(fastaFileNames[i], numSeqs); //fills MPIPos, returns numSeqs
					
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
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}   outputTypes.clear(); delete chimera; return 0;  }
					
					//wait on chidren
					for(int j = 1; j < processors; j++) { 
						char buf[5];
						MPI_Recv(buf, 5, MPI_CHAR, j, tag, MPI_COMM_WORLD, &status); 
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
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  outputTypes.clear(); delete chimera; return 0;  }
					
					//tell parent you are done.
					char buf[5];
					strcpy(buf, "done"); 
					MPI_Send(buf, 5, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
				}
				
				//close files 
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPI);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
		#else
			
			vector<unsigned long int> positions = m->divideFile(fastaFileNames[i], processors);
				
			for (int s = 0; s < (positions.size()-1); s++) {
				lines.push_back(new linePair(positions[s], positions[(s+1)]));
			}	
			
			//break up file
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					numSeqs = driver(lines[0], outputFileName, fastaFileNames[i]);
					
					if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} for (int j = 0; j < lines.size(); j++) {  delete lines[j];  } outputTypes.clear();  lines.clear(); delete chimera; return 0; }
									
				}else{
					processIDS.resize(0);
					
					numSeqs = createProcesses(outputFileName, fastaFileNames[i]); 
				
					rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
						
					//append output files
					for(int j=1;j<processors;j++){
						m->appendFiles((outputFileName + toString(processIDS[j]) + ".temp"), outputFileName);
						remove((outputFileName + toString(processIDS[j]) + ".temp").c_str());
					}
					
					if (m->control_pressed) { 
						for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} outputTypes.clear();
						for (int j = 0; j < lines.size(); j++) {  delete lines[j];  }  lines.clear();
						delete chimera;
						return 0;
					}
				}

			#else
				numSeqs = driver(lines[0], outputFileName, fastaFileNames[i]);
				
				if (m->control_pressed) { for (int j = 0; j < lines.size(); j++) {  delete lines[j];  }  lines.clear(); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} outputTypes.clear(); delete chimera; return 0; }
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

int ChimeraCheckCommand::driver(linePair* filePos, string outputFName, string filename){
	try {
		ofstream out;
		m->openOutputFile(outputFName, out);
		
		ofstream out2;
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
	
		while (!done) {

			if (m->control_pressed) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  m->gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				//find chimeras
				chimera->getChimeras(candidateSeq);
				
				if (m->control_pressed) {	delete candidateSeq; return 1;	}
	
				//print results
				chimera->print(out, out2);
			}
			delete candidateSeq;
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				unsigned long int pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		
		out.close();
		inFASTA.close();
				
		return count;
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

			Sequence* candidateSeq = new Sequence(iss);  m->gobble(iss);
				
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
		int num = 0;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], outputFileName + toString(getpid()) + ".temp", filename);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(getpid()) + ".num.temp";
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
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); remove(tempFile.c_str());
		}
		
		return num;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/


