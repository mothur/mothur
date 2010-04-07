/*
 *  chimeraccodecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraccodecommand.h"
#include "ccode.h"

//***************************************************************************************************************

ChimeraCcodeCommand::ChimeraCcodeCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "filter", "processors", "window", "template", "mask", "numwanted", "outputdir","inputdir", };
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
			}

			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.ccode command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			templatefile = validParameter.validFile(parameters, "template", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = ""; m->mothurOut("template is a required parameter for the chimera.ccode command."); m->mothurOutEndLine(); abort = true;  }		
			
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
			
			temp = validParameter.validFile(parameters, "numwanted", false);		if (temp == "not found") { temp = "20"; }
			convert(temp, numwanted);

		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "ChimeraCcodeCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraCcodeCommand::help(){
	try {
	
		m->mothurOut("The chimera.ccode command reads a fastafile and templatefile and outputs potentially chimeric sequences.\n");
		m->mothurOut("This command was created using the algorythms described in the 'Evaluating putative chimeric sequences from PCR-amplified products' paper by Juan M. Gonzalez, Johannes Zimmerman and Cesareo Saiz-Jimenez.\n");
		m->mothurOut("The chimera.ccode command parameters are fasta, template, filter, mask, processors, window and numwanted.\n");
		m->mothurOut("The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required. \n");
		m->mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences, and is required. \n");
		m->mothurOut("The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
		m->mothurOut("The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences. \n");
		m->mothurOut("The window parameter allows you to specify the window size for searching for chimeras. \n");
		m->mothurOut("The numwanted parameter allows you to specify how many sequences you would each query sequence compared with.\n");
		m->mothurOut("The chimera.ccode command should be in the following format: \n");
		m->mothurOut("chimera.ccode(fasta=yourFastaFile, template=yourTemplate) \n");
		m->mothurOut("Example: chimera.seqs(fasta=AD.align, template=core_set_aligned.imputed.fasta) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraCcodeCommand::~ChimeraCcodeCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraCcodeCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		int start = time(NULL);	
		
		//set user options
		if (maskfile == "default") { m->mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); m->mothurOutEndLine();  }

		chimera = new Ccode(fastafile, templatefile, filter, maskfile, window, numwanted, outputDir);	
		
		//is your template aligned?
		if (chimera->getUnaligned()) { m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); delete chimera; return 0; }
		templateSeqsLength = chimera->getLength();
		
		string outputFileName, accnosFileName;
		if (maskfile != "") {
			outputFileName = outputDir + getRootName(getSimpleName(fastafile)) + maskfile + ".ccode.chimeras";
			accnosFileName = outputDir + getRootName(getSimpleName(fastafile)) + maskfile + ".ccode.accnos";
		}else {
			outputFileName = outputDir + getRootName(getSimpleName(fastafile))  + "ccode.chimeras";
			accnosFileName = outputDir + getRootName(getSimpleName(fastafile))  + "ccode.accnos";
		}

		string mapInfo = outputDir + getRootName(getSimpleName(fastafile)) + "mapinfo";
		bool hasAccnos = true;
		
		if (m->control_pressed) { delete chimera;	return 0;	}
		
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
							
			char outFilename[outputFileName.length()];
			strcpy(outFilename, outputFileName.c_str());
			
			char outAccnosFilename[accnosFileName.length()];
			strcpy(outAccnosFilename, accnosFileName.c_str());
			
			char inFileName[fastafile.length()];
			strcpy(inFileName, fastafile.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
			MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPI);
			MPI_File_open(MPI_COMM_WORLD, outAccnosFilename, outMode, MPI_INFO_NULL, &outMPIAccnos);
		
			if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos);  delete chimera; return 0;  }
		
			if (pid == 0) { //you are the root process 
				string outTemp = "For full window mapping info refer to " + mapInfo + "\n\n";
				
				//print header
				int length = outTemp.length();
				char buf2[length];
				strcpy(buf2, outTemp.c_str()); 
				MPI_File_write_shared(outMPI, buf2, length, MPI_CHAR, &status);
				
				MPIPos = setFilePosFasta(fastafile, numSeqs); //fills MPIPos, returns numSeqs
				
				//send file positions to all processes
				MPI_Bcast(&numSeqs, 1, MPI_INT, 0, MPI_COMM_WORLD);  //send numSeqs
				MPI_Bcast(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, MPI_COMM_WORLD); //send file pos	
				
				//figure out how many sequences you have to align
				numSeqsPerProcessor = numSeqs / processors;
				if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				int startIndex =  pid * numSeqsPerProcessor;
			
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
				if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				int startIndex =  pid * numSeqsPerProcessor;
				
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
		ofstream outHeader;
		string tempHeader = outputDir + getRootName(getSimpleName(fastafile)) + maskfile + "ccode.chimeras.tempHeader";
		openOutputFile(tempHeader, outHeader);
		
		outHeader << "For full window mapping info refer to " << mapInfo << endl << endl;

		outHeader.close();
		
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
					remove(tempHeader.c_str()); 
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
			openInputFile(candidateFileNames[s], inFASTA);
			numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
			inFASTA.close();
			lines.push_back(new linePair(0, numSeqs));
			
			driver(lines[0], outputFileName, fastafile, accnosFileName);
			
			if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(tempHeader.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
			}
			
			//delete accnos file if its blank 
			if (isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
		#endif
		
		//m->mothurOut("Output File Names: ");
		//if ((filter) && (method == "bellerophon")) { m->mothurOut(
		//if (outputDir == "") { fastafile = getRootName(fastafile) + "filter.fasta"; }
		//	else				 { fastafile = outputDir + getRootName(getSimpleName(fastafile)) + "filter.fasta"; }
	
		appendFiles(outputFileName, tempHeader);
	
		remove(outputFileName.c_str());
		rename(tempHeader.c_str(), outputFileName.c_str());
	#endif
	
		delete chimera;
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		if (hasAccnos) {  m->mothurOut(accnosFileName); m->mothurOutEndLine();  }
		m->mothurOutEndLine();

		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		
		m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraCcodeCommand::driver(linePair* line, string outputFName, string filename, string accnos){
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
				
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
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
		m->errorOut(e, "ChimeraCcodeCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ChimeraCcodeCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& outMPI, MPI_File& outAccMPI, vector<long>& MPIPos){
	try {
				
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		for(int i=0;i<num;i++){
		
			if (m->control_pressed) { return 0; }
			
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];
	
			char buf4[length];
			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);

			Sequence* candidateSeq = new Sequence(iss);  gobble(iss);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
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
		m->errorOut(e, "ChimeraCcodeCommand", "driverMPI");
		exit(1);
	}
}
#endif

/**************************************************************************************************/

int ChimeraCcodeCommand::createProcesses(string outputFileName, string filename, string accnos) {
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
		m->errorOut(e, "ChimeraCcodeCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************

