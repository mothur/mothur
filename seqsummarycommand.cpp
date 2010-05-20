/*
 *  seqcoordcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "seqsummarycommand.h"
#include "sequence.hpp"

//***************************************************************************************************************

SeqSummaryCommand::SeqSummaryCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","processors","outputdir","inputdir"};
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
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the summary.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 


		}
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SeqSummaryCommand::help(){
	try {
		m->mothurOut("The summary.seqs command reads a fastafile and summarizes the sequences.\n");
		m->mothurOut("The summary.seqs command parameters are fasta and processors, fasta is required.\n");
		m->mothurOut("The summary.seqs command should be in the following format: \n");
		m->mothurOut("summary.seqs(fasta=yourFastaFile, processors=2) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

SeqSummaryCommand::~SeqSummaryCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int SeqSummaryCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		string summaryFile = outputDir + getSimpleName(fastafile) + ".summary";
				
		int numSeqs = 0;
		
		vector<int> startPosition;
		vector<int> endPosition;
		vector<int> seqLength;
		vector<int> ambigBases;
		vector<int> longHomoPolymer;
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor; 
				int tag = 2001;
				int startTag = 1; int endTag = 2; int lengthTag = 3; int baseTag = 4; int lhomoTag = 5;
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				vector<long> MPIPos;
				
				MPI_Status status; 
				MPI_Status statusOut;
				MPI_File inMPI; 
				MPI_File outMPI; 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
							
				char tempFileName[1024];
				strcpy(tempFileName, fastafile.c_str());
				
				char sumFileName[1024];
				strcpy(sumFileName, summaryFile.c_str());
		
				MPI_File_open(MPI_COMM_WORLD, tempFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, sumFileName, outMode, MPI_INFO_NULL, &outMPI);
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI); return 0;  }
				
				if (pid == 0) { //you are the root process
						//print header
						string outputString = "seqname\tstart\tend\tnbases\tambigs\tpolymer\n";	
						int length = outputString.length();
						char* buf2 = new char[length];
						memcpy(buf2, outputString.c_str(), length);
					
						MPI_File_write_shared(outMPI, buf2, length, MPI_CHAR, &statusOut);
						delete buf2;
						
						MPIPos = setFilePosFasta(fastafile, numSeqs); //fills MPIPos, returns numSeqs
					
						for(int i = 1; i < processors; i++) { 
							MPI_Send(&numSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
							MPI_Send(&MPIPos[0], (numSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
						}
						
						//figure out how many sequences you have to do
						numSeqsPerProcessor = numSeqs / processors;
						int startIndex =  pid * numSeqsPerProcessor;
						if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
						
						//do your part
						MPICreateSummary(startIndex, numSeqsPerProcessor, startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, inMPI, outMPI, MPIPos);
						
				}else { //i am the child process
			
					MPI_Recv(&numSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(numSeqs+1);
					MPI_Recv(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
				
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				
					//do your part
					MPICreateSummary(startIndex, numSeqsPerProcessor, startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, inMPI, outMPI, MPIPos);
				}
				
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPI);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
				if (pid == 0) {
					//get the info from the child processes
					for(int i = 1; i < processors; i++) { 
						int size;
						MPI_Recv(&size, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

						vector<int> temp; temp.resize(size+1);
						
						for(int j = 0; j < 5; j++) { 
						
							MPI_Recv(&temp[0], (size+1), MPI_INT, i, 2001, MPI_COMM_WORLD, &status); 
							int receiveTag = temp[temp.size()-1];  //child process added a int to the end to indicate what count this is for
							
							if (receiveTag == startTag) { 
								for (int k = 0; k < size; k++) {		startPosition.push_back(temp[k]);	}
							}else if (receiveTag == endTag) { 
								for (int k = 0; k < size; k++) {		endPosition.push_back(temp[k]);	}
							}else if (receiveTag == lengthTag) { 
								for (int k = 0; k < size; k++) {		seqLength.push_back(temp[k]);	}
							}else if (receiveTag == baseTag) { 
								for (int k = 0; k < size; k++) {		ambigBases.push_back(temp[k]);	}
							}else if (receiveTag == lhomoTag) { 
								for (int k = 0; k < size; k++) {		longHomoPolymer.push_back(temp[k]);	}
							}
						} 
					}

				}else{
				
					//send my counts
					int size = startPosition.size();
					MPI_Send(&size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
					
					startPosition.push_back(startTag);
					int ierr = MPI_Send(&(startPosition[0]), (size+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
					endPosition.push_back(endTag);
					ierr = MPI_Send (&(endPosition[0]), (size+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
					seqLength.push_back(lengthTag);
					ierr = MPI_Send(&(seqLength[0]), (size+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
					ambigBases.push_back(baseTag);
					ierr = MPI_Send(&(ambigBases[0]), (size+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
					longHomoPolymer.push_back(lhomoTag);
					ierr = MPI_Send(&(longHomoPolymer[0]), (size+1), MPI_INT, 0, 2001, MPI_COMM_WORLD);
				}
				
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#else
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					ifstream inFASTA;
					openInputFile(fastafile, inFASTA);
					numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
					inFASTA.close();
					
					lines.push_back(new linePair(0, numSeqs));
					
					driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile, lines[0]);
				}else{
					numSeqs = setLines(fastafile);					
					createProcessesCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile); 
					
					rename((summaryFile + toString(processIDS[0]) + ".temp").c_str(), summaryFile.c_str());
					//append files
					for(int i=1;i<processors;i++){
						appendFiles((summaryFile + toString(processIDS[i]) + ".temp"), summaryFile);
						remove((summaryFile + toString(processIDS[i]) + ".temp").c_str());
					}
				}
				
				if (m->control_pressed) {  return 0; }
		#else
				ifstream inFASTA;
				openInputFile(fastafile, inFASTA);
				numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
				inFASTA.close();
				
				lines.push_back(new linePair(0, numSeqs));
				
				driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile, lines[0]);
				if (m->control_pressed) {  return 0; }
		#endif
#endif
			
		#ifdef USE_MPI
			if (pid == 0) { 
		#endif
		
		sort(startPosition.begin(), startPosition.end());
		sort(endPosition.begin(), endPosition.end());
		sort(seqLength.begin(), seqLength.end());
		sort(ambigBases.begin(), ambigBases.end());
		sort(longHomoPolymer.begin(), longHomoPolymer.end());
		
		int ptile0_25	= int(numSeqs * 0.025);
		int ptile25		= int(numSeqs * 0.250);
		int ptile50		= int(numSeqs * 0.500);
		int ptile75		= int(numSeqs * 0.750);
		int ptile97_5	= int(numSeqs * 0.975);
		int ptile100	= numSeqs - 1;
		
		//to compensate for blank sequences that would result in startPosition and endPostion equalling -1
		if (startPosition[0] == -1) {  startPosition[0] = 0;	}
		if (endPosition[0] == -1)	{  endPosition[0] = 0;		}
		
		if (m->control_pressed) {  remove(summaryFile.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("\t\tStart\tEnd\tNBases\tAmbigs\tPolymer"); m->mothurOutEndLine();
		m->mothurOut("Minimum:\t" + toString(startPosition[0]) + "\t" + toString(endPosition[0]) + "\t" + toString(seqLength[0]) + "\t" + toString(ambigBases[0]) + "\t" + toString(longHomoPolymer[0])); m->mothurOutEndLine();
		m->mothurOut("2.5%-tile:\t" + toString(startPosition[ptile0_25]) + "\t" + toString(endPosition[ptile0_25]) + "\t" + toString(seqLength[ptile0_25]) + "\t" + toString(ambigBases[ptile0_25]) + "\t"+ toString(longHomoPolymer[ptile0_25])); m->mothurOutEndLine();
		m->mothurOut("25%-tile:\t" + toString(startPosition[ptile25]) + "\t" + toString(endPosition[ptile25]) + "\t" + toString(seqLength[ptile25]) + "\t" + toString(ambigBases[ptile25]) + "\t" + toString(longHomoPolymer[ptile25])); m->mothurOutEndLine();
		m->mothurOut("Median: \t" + toString(startPosition[ptile50]) + "\t" + toString(endPosition[ptile50]) + "\t" + toString(seqLength[ptile50]) + "\t" + toString(ambigBases[ptile50]) + "\t" + toString(longHomoPolymer[ptile50])); m->mothurOutEndLine();
		m->mothurOut("75%-tile:\t" + toString(startPosition[ptile75]) + "\t" + toString(endPosition[ptile75]) + "\t" + toString(seqLength[ptile75]) + "\t" + toString(ambigBases[ptile75]) + "\t" + toString(longHomoPolymer[ptile75])); m->mothurOutEndLine();
		m->mothurOut("97.5%-tile:\t" + toString(startPosition[ptile97_5]) + "\t" + toString(endPosition[ptile97_5]) + "\t" + toString(seqLength[ptile97_5]) + "\t" + toString(ambigBases[ptile97_5]) + "\t" + toString(longHomoPolymer[ptile97_5])); m->mothurOutEndLine();
		m->mothurOut("Maximum:\t" + toString(startPosition[ptile100]) + "\t" + toString(endPosition[ptile100]) + "\t" + toString(seqLength[ptile100]) + "\t" + toString(ambigBases[ptile100]) + "\t" + toString(longHomoPolymer[ptile100])); m->mothurOutEndLine();
		m->mothurOut("# of Seqs:\t" + toString(numSeqs)); m->mothurOutEndLine();
		
		if (m->control_pressed) {  remove(summaryFile.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(summaryFile); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		#ifdef USE_MPI
			}
		#endif

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
int SeqSummaryCommand::driverCreateSummary(vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, string filename, string sumFile, linePair* line) {	
	try {
		
		ofstream outSummary;
		openOutputFile(sumFile, outSummary);
		
		//print header if you are process 0
		if (line->start == 0) {
			outSummary << "seqname\tstart\tend\tnbases\tambigs\tpolymer" << endl;	
		}
				
		ifstream in;
		openInputFile(filename, in);
				
		in.seekg(line->start);
		
		for(int i=0;i<line->num;i++){
				
			if (m->control_pressed) { in.close(); outSummary.close(); return 1; }
					
			Sequence current(in);
			if (current.getName() != "") {
				startPosition.push_back(current.getStartPos());
				endPosition.push_back(current.getEndPos());
				seqLength.push_back(current.getNumBases());
				ambigBases.push_back(current.getAmbigBases());
				longHomoPolymer.push_back(current.getLongHomoPolymer());
				
				outSummary << current.getName() << '\t';
				outSummary << current.getStartPos() << '\t' << current.getEndPos() << '\t';
				outSummary << current.getNumBases() << '\t' << current.getAmbigBases() << '\t';
				outSummary << current.getLongHomoPolymer() << endl;
			}
			gobble(in);
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "driverCreateSummary");
		exit(1);
	}
}
#ifdef USE_MPI
/**************************************************************************************/
int SeqSummaryCommand::MPICreateSummary(int start, int num, vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, MPI_File& inMPI, MPI_File& outMPI, vector<long>& MPIPos) {	
	try {
		
		int pid;
		MPI_Status status; 
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 

		
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

			Sequence current(iss);  

			if (current.getName() != "") {
				startPosition.push_back(current.getStartPos());
				endPosition.push_back(current.getEndPos());
				seqLength.push_back(current.getNumBases());
				ambigBases.push_back(current.getAmbigBases());
				longHomoPolymer.push_back(current.getLongHomoPolymer());
				
				string outputString = current.getName() + "\t" + toString(current.getStartPos()) + "\t" + toString(current.getEndPos()) + "\t";
				outputString += toString(current.getNumBases()) + "\t" + toString(current.getAmbigBases()) + "\t" + toString(current.getLongHomoPolymer()) + "\n";
				
				//output to file
				length = outputString.length();
				char* buf3 = new char[length];
				memcpy(buf3, outputString.c_str(), length);
					
				MPI_File_write_shared(outMPI, buf3, length, MPI_CHAR, &status);
				delete buf3;
			}	
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "MPICreateSummary");
		exit(1);
	}
}
#endif
/**************************************************************************************************/
int SeqSummaryCommand::createProcessesCreateSummary(vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, string filename, string sumFile) {
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
				driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, sumFile + toString(getpid()) + ".temp", lines[process]);
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
		m->errorOut(e, "SeqSummaryCommand", "createProcessesCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/

int SeqSummaryCommand::setLines(string filename) {
	try {
		
		vector<long int> positions;
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);
			
		string input;
		while(!inFASTA.eof()){	
			input = getline(inFASTA);

			if (input.length() != 0) {
				if(input[0] == '>'){ long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
			}
		}
		inFASTA.close();
		
		int numFastaSeqs = positions.size();
	
		FILE * pFile;
		long size;
		
		//get num bytes in file
		pFile = fopen (filename.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
		
		int numSeqsPerProcessor = numFastaSeqs / processors;
		
		for (int i = 0; i < processors; i++) {

			long int startPos = positions[ i * numSeqsPerProcessor ];
			if(i == processors - 1){
				numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;
			}else{  
				long int myEnd = positions[ (i+1) * numSeqsPerProcessor ];
			}
			lines.push_back(new linePair(startPos, numSeqsPerProcessor));
		}
		
		return numFastaSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "setLines");
		exit(1);
	}
}
//***************************************************************************************************************


