/*
 *  seqcoordcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "seqsummarycommand.h"


//**********************************************************************************************************************
vector<string> SeqSummaryCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SeqSummaryCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The summary.seqs command reads a fastafile and summarizes the sequences.\n";
		helpString += "The summary.seqs command parameters are fasta, name and processors, fasta is required, unless you have a valid current fasta file.\n";
		helpString += "The name parameter allows you to enter a name file associated with your fasta file. \n";
		helpString += "The summary.seqs command should be in the following format: \n";
		helpString += "summary.seqs(fasta=yourFastaFile, processors=2) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
SeqSummaryCommand::SeqSummaryCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SeqSummaryCommand::SeqSummaryCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
			else if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { m->setNameFile(namefile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);


		}
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SeqSummaryCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//set current fasta to fastafile
		m->setFastaFile(fastafile);
		
		string summaryFile = outputDir + m->getSimpleName(fastafile) + ".summary";
				
		int numSeqs = 0;
		
		vector<int> startPosition;
		vector<int> endPosition;
		vector<int> seqLength;
		vector<int> ambigBases;
		vector<int> longHomoPolymer;
		
		if (namefile != "") { nameMap = m->readNames(namefile); }
		
		if (m->control_pressed) { return 0; }
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor; 
				int tag = 2001;
				int startTag = 1; int endTag = 2; int lengthTag = 3; int baseTag = 4; int lhomoTag = 5;
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				vector<unsigned long long> MPIPos;
				
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
			vector<unsigned long long> positions; 
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				positions = m->divideFile(fastafile, processors);
				for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
			#else
				positions = m->setFilePosFasta(fastafile, numSeqs); 
		
				//figure out how many sequences you have to process
				int numSeqsPerProcessor = numSeqs / processors;
				for (int i = 0; i < processors; i++) {
					int startIndex =  i * numSeqsPerProcessor;
					if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
					lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
				}
			#endif
			

			if(processors == 1){
				numSeqs = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile, lines[0]);
			}else{
				numSeqs = createProcessesCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, summaryFile); 
			}
			
			if (m->control_pressed) {  return 0; }
#endif
			
		#ifdef USE_MPI
			if (pid == 0) { 
		#endif
		
		sort(startPosition.begin(), startPosition.end());
		sort(endPosition.begin(), endPosition.end());
		sort(seqLength.begin(), seqLength.end());
		sort(ambigBases.begin(), ambigBases.end());
		sort(longHomoPolymer.begin(), longHomoPolymer.end());
		int size = startPosition.size();
		
		//find means
		float meanStartPosition, meanEndPosition, meanSeqLength, meanAmbigBases, meanLongHomoPolymer;
		meanStartPosition = 0; meanEndPosition = 0; meanSeqLength = 0; meanAmbigBases = 0; meanLongHomoPolymer = 0;
		for (int i = 0; i < size; i++) {
			meanStartPosition += startPosition[i];
			meanEndPosition += endPosition[i];
			meanSeqLength += seqLength[i];
			meanAmbigBases += ambigBases[i];
			meanLongHomoPolymer += longHomoPolymer[i];
		}
		//this is an int divide so the remainder is lost
		meanStartPosition /= (float) size; meanEndPosition /= (float) size; meanLongHomoPolymer /= (float) size; meanSeqLength /= (float) size; meanAmbigBases /= (float) size;
				
		int ptile0_25	= int(size * 0.025);
		int ptile25		= int(size * 0.250);
		int ptile50		= int(size * 0.500);
		int ptile75		= int(size * 0.750);
		int ptile97_5	= int(size * 0.975);
		int ptile100	= size - 1;
		
		//to compensate for blank sequences that would result in startPosition and endPostion equalling -1
		if (startPosition[0] == -1) {  startPosition[0] = 0;	}
		if (endPosition[0] == -1)	{  endPosition[0] = 0;		}
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("\t\tStart\tEnd\tNBases\tAmbigs\tPolymer\tNumSeqs"); m->mothurOutEndLine();
		m->mothurOut("Minimum:\t" + toString(startPosition[0]) + "\t" + toString(endPosition[0]) + "\t" + toString(seqLength[0]) + "\t" + toString(ambigBases[0]) + "\t" + toString(longHomoPolymer[0]) + "\t" + toString(1)); m->mothurOutEndLine();
		m->mothurOut("2.5%-tile:\t" + toString(startPosition[ptile0_25]) + "\t" + toString(endPosition[ptile0_25]) + "\t" + toString(seqLength[ptile0_25]) + "\t" + toString(ambigBases[ptile0_25]) + "\t"+ toString(longHomoPolymer[ptile0_25]) + "\t" + toString(ptile0_25+1)); m->mothurOutEndLine();
		m->mothurOut("25%-tile:\t" + toString(startPosition[ptile25]) + "\t" + toString(endPosition[ptile25]) + "\t" + toString(seqLength[ptile25]) + "\t" + toString(ambigBases[ptile25]) + "\t" + toString(longHomoPolymer[ptile25]) + "\t" + toString(ptile25+1)); m->mothurOutEndLine();
		m->mothurOut("Median: \t" + toString(startPosition[ptile50]) + "\t" + toString(endPosition[ptile50]) + "\t" + toString(seqLength[ptile50]) + "\t" + toString(ambigBases[ptile50]) + "\t" + toString(longHomoPolymer[ptile50]) + "\t" + toString(ptile50+1)); m->mothurOutEndLine();
		m->mothurOut("75%-tile:\t" + toString(startPosition[ptile75]) + "\t" + toString(endPosition[ptile75]) + "\t" + toString(seqLength[ptile75]) + "\t" + toString(ambigBases[ptile75]) + "\t" + toString(longHomoPolymer[ptile75]) + "\t" + toString(ptile75+1)); m->mothurOutEndLine();
		m->mothurOut("97.5%-tile:\t" + toString(startPosition[ptile97_5]) + "\t" + toString(endPosition[ptile97_5]) + "\t" + toString(seqLength[ptile97_5]) + "\t" + toString(ambigBases[ptile97_5]) + "\t" + toString(longHomoPolymer[ptile97_5]) + "\t" + toString(ptile97_5+1)); m->mothurOutEndLine();
		m->mothurOut("Maximum:\t" + toString(startPosition[ptile100]) + "\t" + toString(endPosition[ptile100]) + "\t" + toString(seqLength[ptile100]) + "\t" + toString(ambigBases[ptile100]) + "\t" + toString(longHomoPolymer[ptile100]) + "\t" + toString(ptile100+1)); m->mothurOutEndLine();
		m->mothurOut("Mean:\t" + toString(meanStartPosition) + "\t" + toString(meanEndPosition) + "\t" + toString(meanSeqLength) + "\t" + toString(meanAmbigBases) + "\t" + toString(meanLongHomoPolymer)); m->mothurOutEndLine();

		if (namefile == "") {  m->mothurOut("# of Seqs:\t" + toString(numSeqs)); m->mothurOutEndLine(); }
		else { m->mothurOut("# of unique seqs:\t" + toString(numSeqs)); m->mothurOutEndLine(); m->mothurOut("total # of seqs:\t" + toString(startPosition.size())); m->mothurOutEndLine(); }
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
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
				unsigned long long pos = in.tellg();
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
int SeqSummaryCommand::MPICreateSummary(int start, int num, vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, MPI_File& inMPI, MPI_File& outMPI, vector<unsigned long long>& MPIPos) {	
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
					
					if (it == nameMap.end()) { cout << "[ERROR]: " << current.getName() << " is not in your namefile, please correct." << endl; m->control_pressed = true; }
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
		int process = 1;
		int num = 0;
		processIDS.clear();
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		
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
		
		//do your part
		num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, fastafile, sumFile, lines[0]);

		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
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
			m->mothurRemove(tempFilename);
			
			m->appendFiles((sumFile + toString(processIDS[i]) + ".temp"), sumFile);
			m->mothurRemove((sumFile + toString(processIDS[i]) + ".temp"));
		}
		
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<seqSumData*> pDataArray; 
		DWORD   dwThreadIdArray[processors];
		HANDLE  hThreadArray[processors]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors; i++ ){
			
			//cout << i << '\t' << lines[i]->start << '\t' << lines[i]->end << endl;
			// Allocate memory for thread data.
			seqSumData* tempSum = new seqSumData(&startPosition, &endPosition, &seqLength, &ambigBases, &longHomoPolymer, filename, (sumFile + toString(i) + ".temp"), m, lines[i]->start, lines[i]->end, namefile, nameMap);
			pDataArray.push_back(tempSum);
			processIDS.push_back(i);
				
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MySeqSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
			
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
		//rename((sumFile + toString(processIDS[0]) + ".temp").c_str(), sumFile.c_str());
		//append files
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((sumFile + toString(processIDS[i]) + ".temp"), sumFile);
			m->mothurRemove((sumFile + toString(processIDS[i]) + ".temp"));
		}
#endif		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqSummaryCommand", "createProcessesCreateSummary");
		exit(1);
	}
}
/**********************************************************************************************************************/



