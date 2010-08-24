/*
 *  chimeraslayercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraslayercommand.h"
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
			
			ValidParameters validParameter("chimera.slayer");
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
				it = parameters.find("template");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["template"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.slayer command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
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
							ableToOpen = m->openInputFile(tryPath, in, "noerror");
							fastaFileNames[i] = tryPath;
						}
					}
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + fastaFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
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
			
			temp = validParameter.validFile(parameters, "minsnp", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, minSNP);

			temp = validParameter.validFile(parameters, "parents", false);			if (temp == "not found") { temp = "3"; }
			convert(temp, parents); 
			
			temp = validParameter.validFile(parameters, "realign", false);			if (temp == "not found") { temp = "f"; }
			realign = m->isTrue(temp); 
			
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
		m->mothurOut("The chimera.slayer command parameters are fasta, template, processors, ksize, window, match, mismatch, divergence. minsim, mincov, minbs, minsnp, parents, search, iters, increment and numwanted.\n"); //realign,
		m->mothurOut("The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required. \n");
		m->mothurOut("You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n");
		m->mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences, and is required. \n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
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
		m->mothurOut("The minsnp parameter allows you to specify percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default: 100) \n");
		m->mothurOut("The search parameter allows you to specify search method for finding the closest parent. Choices are distance, blast, and kmer, default distance. \n");
		m->mothurOut("The realign parameter allows you to realign the query to the potential parents. Choices are true or false, default false.  \n");
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
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
		
			int start = time(NULL);	
			
			chimera = new ChimeraSlayer(fastaFileNames[s], templatefile, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign);	
							
			string outputFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "slayer.chimeras";
			string accnosFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]))  + "slayer.accnos";
			
			if (m->control_pressed) { delete chimera; for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0;	}
			
			if (chimera->getUnaligned()) { 
				m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); 
				delete chimera;
				return 0; 
			}
			templateSeqsLength = chimera->getLength();
			
		#ifdef USE_MPI	
			int pid, end, numSeqsPerProcessor; 
				int tag = 2001;
				vector<unsigned long int> MPIPos;
				
				MPI_Status status; 
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				MPI_Comm_size(MPI_COMM_WORLD, &processors); 

				MPI_File inMPI;
				MPI_File outMPI;
				MPI_File outMPIAccnos;
				
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				int inMode=MPI_MODE_RDONLY; 
				
				char outFilename[1024];
				strcpy(outFilename, outputFileName.c_str());
				
				char outAccnosFilename[1024];
				strcpy(outAccnosFilename, accnosFileName.c_str());
				
				char inFileName[1024];
				strcpy(inFileName, fastaFileNames[s].c_str());

				MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPI);
				MPI_File_open(MPI_COMM_WORLD, outAccnosFilename, outMode, MPI_INFO_NULL, &outMPIAccnos);

				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}   delete chimera; return 0;  }
			
				if (pid == 0) { //you are the root process 
					m->mothurOutEndLine();
					m->mothurOut("Only reporting sequence supported by " + toString(minBS) + "% of bootstrapped results.");
					m->mothurOutEndLine();
		
					string outTemp = "Name\tLeftParent\tRightParent\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
					
					//print header
					int length = outTemp.length();
					char* buf2 = new char[length];
					memcpy(buf2, outTemp.c_str(), length);

					MPI_File_write_shared(outMPI, buf2, length, MPI_CHAR, &status);
					delete buf2;

					MPIPos = m->setFilePosFasta(fastaFileNames[s], numSeqs); //fills MPIPos, returns numSeqs
					
					//send file positions to all processes
					for(int i = 1; i < processors; i++) { 
						MPI_Send(&numSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&MPIPos[0], (numSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
					}
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				
					//do your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, outMPIAccnos, MPIPos);
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos);  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  remove(outputFileName.c_str());  remove(accnosFileName.c_str());  delete chimera; return 0;  }

				}else{ //you are a child process
					MPI_Recv(&numSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(numSeqs+1);
					MPI_Recv(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
					
					//do your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, outMPIAccnos, MPIPos);
					
					if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI);   MPI_File_close(&outMPIAccnos);  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  delete chimera; return 0;  }
				}
				
				//close files 
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPI);
				MPI_File_close(&outMPIAccnos);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
		#else
			ofstream outHeader;
			string tempHeader = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "slayer.chimeras.tempHeader";
			m->openOutputFile(tempHeader, outHeader);
			
			chimera->printHeader(outHeader);
			outHeader.close();
			
			vector<unsigned long int> positions = m->divideFile(fastaFileNames[s], processors);
				
			for (int i = 0; i < (positions.size()-1); i++) {
				lines.push_back(new linePair(positions[i], positions[(i+1)]));
			}	

			//break up file
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					numSeqs = driver(lines[0], outputFileName, fastaFileNames[s], accnosFileName);
					
					if (m->control_pressed) { remove(outputFileName.c_str()); remove(tempHeader.c_str()); remove(accnosFileName.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear(); delete chimera; return 0; }
					
				}else{
					processIDS.resize(0);
					
					numSeqs = createProcesses(outputFileName, fastaFileNames[s], accnosFileName); 
				
					rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
					rename((accnosFileName + toString(processIDS[0]) + ".temp").c_str(), accnosFileName.c_str());
						
					//append output files
					for(int i=1;i<processors;i++){
						m->appendFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
						remove((outputFileName + toString(processIDS[i]) + ".temp").c_str());
					}
					
					//append output files
					for(int i=1;i<processors;i++){
						m->appendFiles((accnosFileName + toString(processIDS[i]) + ".temp"), accnosFileName);
						remove((accnosFileName + toString(processIDS[i]) + ".temp").c_str());
					}
					
					if (m->control_pressed) { remove(outputFileName.c_str()); remove(accnosFileName.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear(); delete chimera; return 0; }
				}

			#else
				numSeqs = driver(lines[0], outputFileName, fastaFileNames[s], accnosFileName);
				
				if (m->control_pressed) { remove(outputFileName.c_str()); remove(tempHeader.c_str()); remove(accnosFileName.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear(); delete chimera; return 0; }
				
			#endif
			
			m->appendFiles(outputFileName, tempHeader);
		
			remove(outputFileName.c_str());
			rename(tempHeader.c_str(), outputFileName.c_str());
			
		#endif
			delete chimera;
			
			
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			
			outputNames.push_back(outputFileName);
			outputNames.push_back(accnosFileName); 
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraSlayerCommand::driver(linePair* filePos, string outputFName, string filename, string accnos){
	try {
		ofstream out;
		m->openOutputFile(outputFName, out);
		
		ofstream out2;
		m->openOutputFile(accnos, out2);
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
	
		while (!done) {
		
			if (m->control_pressed) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  m->gobble(inFASTA);
				
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
			count++;
			}
			delete candidateSeq;
			
			unsigned long int pos = inFASTA.tellg();
			if ((pos == -1) || (pos >= filePos->end)) { break; }
			
			//report progress
			if((count) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		
		out.close();
		out2.close();
		inFASTA.close();
				
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ChimeraSlayerCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& outMPI, MPI_File& outAccMPI, vector<unsigned long int>& MPIPos){
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

			Sequence* candidateSeq = new Sequence(iss);  m->gobble(iss);
		
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
		
					//find chimeras
					chimera->getChimeras(candidateSeq);
			
					if (m->control_pressed) {	delete candidateSeq; return 1;	}
		//cout << "about to print" << endl;
					//print results
					bool isChimeric = chimera->print(outMPI, outAccMPI);
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
		int num = 0;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], outputFileName + toString(getpid()) + ".temp", filename, accnos + toString(getpid()) + ".temp");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();
				
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
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
		m->errorOut(e, "ChimeraSlayerCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/


