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

//**********************************************************************************************************************
vector<string> SeqSummaryCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta","name","processors","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SeqSummaryCommand::SeqSummaryCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SeqSummaryCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SeqSummaryCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "getRequiredFiles");
		exit(1);
	}
}
//***************************************************************************************************************

SeqSummaryCommand::SeqSummaryCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","name","processors","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("summary.seqs");
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
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the summary.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
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
		m->mothurOut("The summary.seqs command parameters are fasta, name and processors, fasta is required.\n");
		m->mothurOut("The name parameter allows you to enter a name file associated with your fasta file. \n");
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
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		string summaryFile = outputDir + m->getSimpleName(fastafile) + ".summary";
				
		int numSeqs = 0;
		
		vector<int> startPosition;
		vector<int> endPosition;
		vector<int> seqLength;
		vector<int> ambigBases;
		vector<int> longHomoPolymer;
		
		if (namefile != "") { readNames(); }
		
		if (m->control_pressed) { return 0; }
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor; 
				int tag = 2001;
				int startTag = 1; int endTag = 2; int lengthTag = 3; int baseTag = 4; int lhomoTag = 5;
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				vector<unsigned long int> MPIPos;
				
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
						string outputString = "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs\n";	
						int length = outputString.length();
						char* buf2 = new char[length];
						memcpy(buf2, outputString.c_str(), length);
					
						MPI_File_write_shared(outMPI, buf2, length, MPI_CHAR, &statusOut);
						delete buf2;
						
						MPIPos = m->setFilePosFasta(fastafile, numSeqs); //fills MPIPos, returns numSeqs
					
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
			vector<unsigned long int> positions = m->divideFile(fastafile, processors);
				
			for (int i = 0; i < (positions.size()-1); i++) {
				lines.push_back(new linePair(positions[i], positions[(i+1)]));
			}	

		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					numSeqs = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile, lines[0]);
				}else{
					numSeqs = createProcessesCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile); 
					
					rename((summaryFile + toString(processIDS[0]) + ".temp").c_str(), summaryFile.c_str());
					//append files
					for(int i=1;i<processors;i++){
						m->appendFiles((summaryFile + toString(processIDS[i]) + ".temp"), summaryFile);
						remove((summaryFile + toString(processIDS[i]) + ".temp").c_str());
					}
				}
				
				if (m->control_pressed) {  return 0; }
		#else
				numSeqs = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile, lines[0]);
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
		if (namefile == "") {  m->mothurOut("# of Seqs:\t" + toString(numSeqs)); m->mothurOutEndLine(); }
		else { m->mothurOut("# of unique seqs:\t" + toString(numSeqs)); m->mothurOutEndLine(); m->mothurOut("total # of seqs:\t" + toString(startPosition.size())); m->mothurOutEndLine(); }
		
		if (m->control_pressed) {  remove(summaryFile.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(summaryFile); m->mothurOutEndLine();	outputNames.push_back(summaryFile); outputTypes["summary"].push_back(summaryFile);
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
int SeqSummaryCommand::driverCreateSummary(vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, string filename, string sumFile, linePair* filePos) {	
	try {
		
		ofstream outSummary;
		m->openOutputFile(sumFile, outSummary);
		
		//print header if you are process 0
		if (filePos->start == 0) {
			outSummary << "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs" << endl;	
		}
				
		ifstream in;
		m->openInputFile(filename, in);
				
		in.seekg(filePos->start);

		bool done = false;
		int count = 0;
	
		while (!done) {
				
			if (m->control_pressed) { in.close(); outSummary.close(); return 1; }
					
			Sequence current(in); m->gobble(in);
	
			if (current.getName() != "") {
				
				int num = 1;
				if (namefile != "") {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(current.getName());
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				//for each sequence this sequence represents
				for (int i = 0; i < num; i++) {
					startPosition.push_back(current.getStartPos());
					endPosition.push_back(current.getEndPos());
					seqLength.push_back(current.getNumBases());
					ambigBases.push_back(current.getAmbigBases());
					longHomoPolymer.push_back(current.getLongHomoPolymer());
				}
				
				count++;
				outSummary << current.getName() << '\t';
				outSummary << current.getStartPos() << '\t' << current.getEndPos() << '\t';
				outSummary << current.getNumBases() << '\t' << current.getAmbigBases() << '\t';
				outSummary << current.getLongHomoPolymer() << '\t' << num << endl;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				unsigned long int pos = in.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (in.eof()) { break; }
			#endif
			
			//report progress
			//if((count) % 100 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		//if((count) % 100 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "driverCreateSummary");
		exit(1);
	}
}
#ifdef USE_MPI
/**************************************************************************************/
int SeqSummaryCommand::MPICreateSummary(int start, int num, vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, MPI_File& inMPI, MPI_File& outMPI, vector<unsigned long int>& MPIPos) {	
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
				
				int num = 1;
				if (namefile != "") {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(current.getName());
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				//for each sequence this sequence represents
				for (int i = 0; i < num; i++) {
					startPosition.push_back(current.getStartPos());
					endPosition.push_back(current.getEndPos());
					seqLength.push_back(current.getNumBases());
					ambigBases.push_back(current.getAmbigBases());
					longHomoPolymer.push_back(current.getLongHomoPolymer());
				}
				
				string outputString = current.getName() + "\t" + toString(current.getStartPos()) + "\t" + toString(current.getEndPos()) + "\t";
				outputString += toString(current.getNumBases()) + "\t" + toString(current.getAmbigBases()) + "\t" + toString(current.getLongHomoPolymer()) + "\t" + toString(num) + "\n";
				
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
		int num = 0;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, sumFile + toString(getpid()) + ".temp", lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = fastafile + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				
				out << num << endl;
				out << startPosition.size() << endl;
				for (int k = 0; k < startPosition.size(); k++)		{		out << startPosition[k] << '\t'; }  out << endl;
				for (int k = 0; k < endPosition.size(); k++)		{		out << endPosition[k] << '\t'; }  out << endl;
				for (int k = 0; k < seqLength.size(); k++)			{		out << seqLength[k] << '\t'; }  out << endl;
				for (int k = 0; k < ambigBases.size(); k++)			{		out << ambigBases[k] << '\t'; }  out << endl;
				for (int k = 0; k < longHomoPolymer.size(); k++)	{		out << longHomoPolymer[k] << '\t'; }  out << endl;
				
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
		
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = fastafile + toString(processIDS[i]) + ".num.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			int temp, tempNum;
			in >> tempNum; m->gobble(in); num += tempNum;
			in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; startPosition.push_back(temp);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; endPosition.push_back(temp);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; seqLength.push_back(temp);			}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; ambigBases.push_back(temp);			}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; longHomoPolymer.push_back(temp);	}		m->gobble(in);
				
			in.close();
			remove(tempFilename.c_str());
		}
		
		return num;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "createProcessesCreateSummary");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SeqSummaryCommand::readNames() { 
	try {
		//open input file
		ifstream in;
		m->openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol >> secondCol; m->gobble(in);
			
			int num = m->getNumNames(secondCol);
			
			nameMap[firstCol] = num;
		}
		in.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "readNames");
		exit(1);
	}
}

/**********************************************************************************************************************/



