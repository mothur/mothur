/*
 *  chimeraslayercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraslayercommand.h"
#include "bellerophon.h"
#include "pintail.h"
#include "ccode.h"
#include "chimeracheckrdp.h"
#include "chimeraslayer.h"


//***************************************************************************************************************

ChimeraSlayerCommand::ChimeraSlayerCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "processors", "window", "template","numwanted", "ksize", "match","mismatch", 
			"divergence", "minsim","mincov","minbs", "minsnp","parents", "iters","outputdir","inputdir", "search","realign" };
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
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.slayer command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			templatefile = validParameter.validFile(parameters, "template", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = "";  m->mothurOut("template is a required parameter for the chimera.slayer command."); m->mothurOutEndLine(); abort = true;  }	
						
			string temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			convert(temp, ksize);
						
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "50"; }			
			convert(temp, window);
			
			temp = validParameter.validFile(parameters, "match", false);			if (temp == "not found") { temp = "5"; }
			convert(temp, match);
			
			temp = validParameter.validFile(parameters, "mismatch", false);			if (temp == "not found") { temp = "-4"; }
			convert(temp, mismatch);
			
			temp = validParameter.validFile(parameters, "divergence", false);		if (temp == "not found") { temp = "1.007"; }
			convert(temp, divR);
			
			temp = validParameter.validFile(parameters, "minsim", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minSimilarity);
			
			temp = validParameter.validFile(parameters, "mincov", false);			if (temp == "not found") { temp = "70"; }
			convert(temp, minCoverage);
			
			temp = validParameter.validFile(parameters, "minbs", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minBS);
			
			temp = validParameter.validFile(parameters, "minsnp", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, minSNP);

			temp = validParameter.validFile(parameters, "parents", false);			if (temp == "not found") { temp = "3"; }
			convert(temp, parents); 
			
			temp = validParameter.validFile(parameters, "realign", false);			if (temp == "not found") { temp = "f"; }
			realign = isTrue(temp); 
			
			search = validParameter.validFile(parameters, "search", false);			if (search == "not found") { search = "distance"; }
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "100"; }		
			convert(temp, iters); 
			 
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "5"; }
			convert(temp, increment);
			
			temp = validParameter.validFile(parameters, "numwanted", false);		if (temp == "not found") { temp = "15"; }		
			convert(temp, numwanted);

			if ((search != "distance") && (search != "blast") && (search != "kmer")) { m->mothurOut(search + " is not a valid search."); m->mothurOutEndLine(); abort = true;  }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "ChimeraSlayerCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraSlayerCommand::help(){
	try {
	
		m->mothurOut("The chimera.slayer command reads a fastafile and templatefile and outputs potentially chimeric sequences.\n");
		m->mothurOut("This command was modeled after the chimeraSlayer written by the Broad Institute.\n");
		m->mothurOut("The chimera.slayer command parameters are fasta, template, filter, mask, processors, ksize, window, match, mismatch, divergence. minsim, mincov, minbs, minsnp, parents, search, iters, increment and numwanted.\n"); //realign,
		m->mothurOut("The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required. \n");
		m->mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences, and is required. \n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
		m->mothurOut("The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences. \n");
		m->mothurOut("The window parameter allows you to specify the window size for searching for chimeras, default=50. \n");
		m->mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default=5.\n");
		m->mothurOut("The numwanted parameter allows you to specify how many sequences you would each query sequence compared with, default=15.\n");
		m->mothurOut("The ksize parameter allows you to input kmersize, default is 7, used if search is kmer. \n");
		m->mothurOut("The match parameter allows you to reward matched bases in blast search, default is 5. \n");
		m->mothurOut("The parents parameter allows you to select the number of potential parents to investigate from the numwanted best matches after rating them, default is 3. \n");
		m->mothurOut("The mismatch parameter allows you to penalize mismatched bases in blast search, default is -4. \n");
		m->mothurOut("The divergence parameter allows you to set a cutoff for chimera determination, default is 1.007. \n");
		m->mothurOut("The iters parameter allows you to specify the number of bootstrap iters to do with the chimeraslayer method, default=100.\n");
		m->mothurOut("The minsim parameter allows you to specify a minimum similarity with the parent fragments, default=90. \n");
		m->mothurOut("The mincov parameter allows you to specify minimum coverage by closest matches found in template. Default is 70, meaning 70%. \n");
		m->mothurOut("The minbs parameter allows you to specify minimum bootstrap support for calling a sequence chimeric. Default is 90, meaning 90%. \n");
		m->mothurOut("The minsnp parameter allows you to specify percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default: 10) \n");
		m->mothurOut("The search parameter allows you to specify search method for finding the closest parent. Choices are distance, blast, and kmer, default distance. \n");
		//m->mothurOut("The realign parameter allows you to realign the query to the potential parents. Choices are true or false, default false. Found to make results worse. \n");
		m->mothurOut("NOT ALL PARAMETERS ARE USED BY ALL METHODS. Please look below for method specifics.\n\n");
		m->mothurOut("The chimera.slayer command should be in the following format: \n");
		m->mothurOut("chimera.slayer(fasta=yourFastaFile, template=yourTemplate, search=yourSearch) \n");
		m->mothurOut("Example: chimera.slayer(fasta=AD.align, template=core_set_aligned.imputed.fasta, search=kmer) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraSlayerCommand::~ChimeraSlayerCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraSlayerCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		int start = time(NULL);	
		
		chimera = new ChimeraSlayer(fastafile, templatefile, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign);	
						
		string outputFileName = outputDir + getRootName(getSimpleName(fastafile)) + "slayer.chimeras";
		string accnosFileName = outputDir + getRootName(getSimpleName(fastafile))  + "slayer.accnos";
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
				m->mothurOutEndLine();
				m->mothurOut("Only reporting sequence supported by " + toString(minBS) + "% of bootstrapped results.");
				m->mothurOutEndLine();
	
				string outTemp = "Name\tLeftParent\tRightParent\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
				
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
		string tempHeader = outputDir + getRootName(getSimpleName(fastafile)) + "slayer.chimeras.tempHeader";
		openOutputFile(tempHeader, outHeader);
		
		chimera->printHeader(outHeader);
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
		
		appendFiles(tempHeader, outputFileName);
	
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
		m->errorOut(e, "ChimeraSlayerCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraSlayerCommand::driver(linePair* line, string outputFName, string filename, string accnos){
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
		m->errorOut(e, "ChimeraSlayerCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ChimeraSlayerCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& outMPI, MPI_File& outAccMPI, vector<long>& MPIPos){
	try {
				
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		for(int i=0;i<num;i++){
			
			if (m->control_pressed) {	return 1;	}
			
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
		m->errorOut(e, "ChimeraSlayerCommand", "driverMPI");
		exit(1);
	}
}
#endif

/**************************************************************************************************/

int ChimeraSlayerCommand::createProcesses(string outputFileName, string filename, string accnos) {
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
		m->errorOut(e, "ChimeraSlayerCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/


