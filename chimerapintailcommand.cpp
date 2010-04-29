/*
 *  chimerapintailcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerapintailcommand.h"
#include "pintail.h"

//***************************************************************************************************************

ChimeraPintailCommand::ChimeraPintailCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","filter","processors","window" ,"increment","template","conservation","quantile","mask","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
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
				
				it = parameters.find("template");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["template"] = inputDir + it->second;		}
				}
				
				it = parameters.find("conservation");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["conservation"] = inputDir + it->second;		}
				}
				
				it = parameters.find("quantile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["quantile"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.pintail command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			templatefile = validParameter.validFile(parameters, "template", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = "";  m->mothurOut("template is a required parameter for the chimera.pintail command."); m->mothurOutEndLine(); abort = true;  }
			
			consfile = validParameter.validFile(parameters, "conservation", true);
			if (consfile == "not open") { abort = true; }
			else if (consfile == "not found") { consfile = "";  }	
			
			quanfile = validParameter.validFile(parameters, "quantile", true);
			if (quanfile == "not open") { abort = true; }
			else if (quanfile == "not found") { quanfile = "";  }
			
			maskfile = validParameter.validFile(parameters, "mask", false);
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  { 
				if (inputDir != "") {
					string path = hasPath(maskfile);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	maskfile = inputDir + maskfile;		}
				}

				ifstream in;
				int	ableToOpen = openInputFile(maskfile, in);
				if (ableToOpen == 1) { abort = true; }
				in.close();
			}
						
			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, window);
			
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "25"; }
			convert(temp, increment);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "ChimeraPintailCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraPintailCommand::help(){
	try {
	
		m->mothurOut("The chimera.pintail command reads a fastafile and templatefile and outputs potentially chimeric sequences.\n");
		m->mothurOut("This command was created using the algorythms described in the 'At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies' paper by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.\n");
		m->mothurOut("The chimera.pintail command parameters are fasta, template, filter, mask, processors, window, increment, conservation and quantile.\n");
		m->mothurOut("The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required. \n");
		m->mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences, and is required. \n");
		m->mothurOut("The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n");
		m->mothurOut("The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences, by default no mask is applied.  You can apply an ecoli mask by typing, mask=default. \n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
		m->mothurOut("The window parameter allows you to specify the window size for searching for chimeras, default=300. \n");
		m->mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default=25.\n");
		m->mothurOut("The conservation parameter allows you to enter a frequency file containing the highest bases frequency at each place in the alignment.\n");
		m->mothurOut("The quantile parameter allows you to enter a file containing quantiles for a template files sequences, if you use the filter the quantile file generated becomes unique to the fasta file you used.\n");
		m->mothurOut("The chimera.pintail command should be in the following format: \n");
		m->mothurOut("chimera.seqs(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors, method=bellerophon) \n");
		m->mothurOut("Example: chimera.seqs(fasta=AD.align, filter=True, correction=true, method=bellerophon, window=200) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraPintailCommand::~ChimeraPintailCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraPintailCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		int start = time(NULL);	
		
		//set user options
		if (maskfile == "default") { m->mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); m->mothurOutEndLine();  }
		
		chimera = new Pintail(fastafile, templatefile, filter, processors, maskfile, consfile, quanfile, window, increment, outputDir);
		
		string outputFileName, accnosFileName;
		if (maskfile != "") {
			outputFileName = outputDir + getRootName(getSimpleName(fastafile)) + maskfile + ".pintail.chimeras";
			accnosFileName = outputDir + getRootName(getSimpleName(fastafile)) + maskfile + ".pintail.accnos";
		}else {
			outputFileName = outputDir + getRootName(getSimpleName(fastafile))  + "pintail.chimeras";
			accnosFileName = outputDir + getRootName(getSimpleName(fastafile))  + "pintail.accnos";
		}
		bool hasAccnos = true;
		
		if (m->control_pressed) { delete chimera;	return 0;	}
		
		if (chimera->getUnaligned()) { 
			m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); 
			delete chimera;
			return 0; 
		}
		templateSeqsLength = chimera->getLength();
	
	#ifdef USE_MPI
		int pid, end, numSeqsPerProcessor; 
			int tag = 2001;
			vector<long> MPIPos;
			MPIWroteAccnos = false;
			
			MPI_Status status; 
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
			MPI_Comm_size(MPI_COMM_WORLD, &processors); 

			MPI_File inMPI;
			MPI_File outMPI;
			MPI_File outMPIAccnos;
			
			int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
			int inMode=MPI_MODE_RDONLY; 
			
			//char* outFilename = new char[outputFileName.length()];
			//memcpy(outFilename, outputFileName.c_str(), outputFileName.length());
			
			char outFilename[1024];
			strcpy(outFilename, outputFileName.c_str());
			
			//char* outAccnosFilename = new char[accnosFileName.length()];
			//memcpy(outAccnosFilename, accnosFileName.c_str(), accnosFileName.length());
			
			char outAccnosFilename[1024];
			strcpy(outAccnosFilename, accnosFileName.c_str());

			//char* inFileName = new char[fastafile.length()];
			//memcpy(inFileName, fastafile.c_str(), fastafile.length());
			
			char inFileName[1024];
			strcpy(inFileName, fastafile.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
			MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPI);
			MPI_File_open(MPI_COMM_WORLD, outAccnosFilename, outMode, MPI_INFO_NULL, &outMPIAccnos);
			
			//delete inFileName;
			//delete outFilename;
			//delete outAccnosFilename;

			if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos);  delete chimera; return 0;  }

			if (pid == 0) { //you are the root process 
							
				MPIPos = setFilePosFasta(fastafile, numSeqs); //fills MPIPos, returns numSeqs
				
				//send file positions to all processes
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD);  //send numSeqs
				MPI_Bcast(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //send file pos	
				
				//figure out how many sequences you have to align
				numSeqsPerProcessor = numSeqs / processors;
				int startIndex =  pid * numSeqsPerProcessor;
				if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				
			
				//align your part
				driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, outMPIAccnos, MPIPos);
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos);  remove(outputFileName.c_str());  remove(accnosFileName.c_str());  delete chimera; return 0;  }
				
				for (int i = 1; i < processors; i++) {
					bool tempResult;
					MPI_Recv(&tempResult, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
					if (tempResult != 0) { MPIWroteAccnos = true; }
				}
			}else{ //you are a child process
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD); //get numSeqs
				MPIPos.resize(numSeqs+1);
				MPI_Bcast(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //get file positions
				
				//figure out how many sequences you have to align
				numSeqsPerProcessor = numSeqs / processors;
				int startIndex =  pid * numSeqsPerProcessor;
				if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				
				
				//align your part
				driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, outMPIAccnos, MPIPos);
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos);  delete chimera; return 0;  }

				MPI_Send(&MPIWroteAccnos, 1, MPI_INT, 0, tag, MPI_COMM_WORLD); 
			}
			
			//close files 
			MPI_File_close(&inMPI);
			MPI_File_close(&outMPI);
			MPI_File_close(&outMPIAccnos);
			
			//delete accnos file if blank
			if (pid == 0) {
				if (!MPIWroteAccnos) { 
					//MPI_Info info;
					//MPI_File_delete(outAccnosFilename, info);
					hasAccnos = false;	
					remove(accnosFileName.c_str()); 
				}
			}

	#else
	
		//break up file
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				ifstream inFASTA;
				openInputFile(fastafile, inFASTA);
				numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
				inFASTA.close();
				
				lines.push_back(new linePair(0, numSeqs));
				
				driver(lines[0], outputFileName, fastafile, accnosFileName);
				
				if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
				}
				
				//delete accnos file if its blank 
				if (isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
								
			}else{
				vector<int> positions;
				processIDS.resize(0);
				
				ifstream inFASTA;
				openInputFile(fastafile, inFASTA);
				
				string input;
				while(!inFASTA.eof()){
					input = getline(inFASTA);
					if (input.length() != 0) {
						if(input[0] == '>'){	long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
					}
				}
				inFASTA.close();
				
				numSeqs = positions.size();
				
				int numSeqsPerProcessor = numSeqs / processors;
				
				for (int i = 0; i < processors; i++) {
					long int startPos = positions[ i * numSeqsPerProcessor ];
					if(i == processors - 1){
						numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor;
					}
					lines.push_back(new linePair(startPos, numSeqsPerProcessor));
				}
				
				
				createProcesses(outputFileName, fastafile, accnosFileName); 
			
				rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
					
				//append output files
				for(int i=1;i<processors;i++){
					appendFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
					remove((outputFileName + toString(processIDS[i]) + ".temp").c_str());
				}
				
				vector<string> nonBlankAccnosFiles;
				//delete blank accnos files generated with multiple processes
				for(int i=0;i<processors;i++){  
					if (!(isBlank(accnosFileName + toString(processIDS[i]) + ".temp"))) {
						nonBlankAccnosFiles.push_back(accnosFileName + toString(processIDS[i]) + ".temp");
					}else { remove((accnosFileName + toString(processIDS[i]) + ".temp").c_str());  }
				}
				
				//append accnos files
				if (nonBlankAccnosFiles.size() != 0) { 
					rename(nonBlankAccnosFiles[0].c_str(), accnosFileName.c_str());
					
					for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
						appendFiles(nonBlankAccnosFiles[h], accnosFileName);
						remove(nonBlankAccnosFiles[h].c_str());
					}
				}else{ hasAccnos = false;  }
				
				if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
				}
			}

		#else
			ifstream inFASTA;
			openInputFile(fastafile, inFASTA);
			numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
			inFASTA.close();
			lines.push_back(new linePair(0, numSeqs));
			
			driver(lines[0], outputFileName, fastafile, accnosFileName);
			
			if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
			}
			
			//delete accnos file if its blank 
			if (isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
		#endif
		
	#endif	
	
		delete chimera;
		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		if (hasAccnos) {  m->mothurOut(accnosFileName); m->mothurOutEndLine();  }
		m->mothurOutEndLine();
		m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraPintailCommand::driver(linePair* line, string outputFName, string filename, string accnos){
	try {
		ofstream out;
		openOutputFile(outputFName, out);
		
		ofstream out2;
		openOutputFile(accnos, out2);
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);

		inFASTA.seekg(line->start);
		
		for(int i=0;i<line->numSeqs;i++){
		
			if (m->control_pressed) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength)  {  //chimeracheck does not require seqs to be aligned
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->control_pressed) {	delete candidateSeq; return 1;	}
		
					//print results
					chimera->print(out, out2);
				}
			}
			delete candidateSeq;
			
			//report progress
			if((i+1) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(i+1)); m->mothurOutEndLine();		}
		}
		//report progress
		if((line->numSeqs) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(line->numSeqs)); m->mothurOutEndLine();		}
		
		out.close();
		out2.close();
		inFASTA.close();
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ChimeraPintailCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& outMPI, MPI_File& outAccMPI, vector<long>& MPIPos){
	try {
				
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		for(int i=0;i<num;i++){
			
			if (m->control_pressed) {	return 1;	}
			
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
				
				if	(candidateSeq->getAligned().length() != templateSeqsLength) {  //chimeracheck does not require seqs to be aligned
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->control_pressed) {	delete candidateSeq; return 1;	}
		
					//print results
					bool isChimeric = chimera->print(outMPI, outAccMPI);
					if (isChimeric) { MPIWroteAccnos = true;  }
				}
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
		m->errorOut(e, "ChimeraPintailCommand", "driverMPI");
		exit(1);
	}
}
#endif

/**************************************************************************************************/

int ChimeraPintailCommand::createProcesses(string outputFileName, string filename, string accnos) {
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
				driver(lines[process], outputFileName + toString(getpid()) + ".temp", filename, accnos + toString(getpid()) + ".temp");
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
		m->errorOut(e, "ChimeraPintailCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/


