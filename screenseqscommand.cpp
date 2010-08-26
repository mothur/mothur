/*
 *  screenseqscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/3/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "screenseqscommand.h"
#include "sequence.hpp"

//***************************************************************************************************************

ScreenSeqsCommand::ScreenSeqsCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta", "start", "end", "maxambig", "maxhomop", "minlength", "maxlength",
									"name", "group", "alignreport","processors","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("screen.seqs");
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("alignreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["alignreport"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { m->mothurOut("fasta is a required parameter for the screen.seqs command."); m->mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
	
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	

			alignreport = validParameter.validFile(parameters, "alignreport", true);
			if (alignreport == "not open") { abort = true; }
			else if (alignreport == "not found") { alignreport = ""; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "start", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, startPos); 
		
			temp = validParameter.validFile(parameters, "end", false);			if (temp == "not found") { temp = "-1"; }
			convert(temp, endPos);  

			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, maxAmbig);  

			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "-1"; }
			convert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "-1"; }
			convert(temp, maxLength); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 

		}

	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "ScreenSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ScreenSeqsCommand::help(){
	try {
		m->mothurOut("The screen.seqs command reads a fastafile and creates .....\n");
		m->mothurOut("The screen.seqs command parameters are fasta, start, end, maxambig, maxhomop, minlength, maxlength, name, group and processors.\n");
		m->mothurOut("The fasta parameter is required.\n");
		m->mothurOut("The start parameter .... The default is -1.\n");
		m->mothurOut("The end parameter .... The default is -1.\n");
		m->mothurOut("The maxambig parameter .... The default is -1.\n");
		m->mothurOut("The maxhomop parameter .... The default is -1.\n");
		m->mothurOut("The minlength parameter .... The default is -1.\n");
		m->mothurOut("The maxlength parameter .... The default is -1.\n");
		m->mothurOut("The processors parameter allows you to specify the number of processors to use while running the command. The default is 1.\n");
		m->mothurOut("The name parameter allows you to provide a namesfile, and the group parameter allows you to provide a groupfile.\n");
		m->mothurOut("The screen.seqs command should be in the following format: \n");
		m->mothurOut("screen.seqs(fasta=yourFastaFile, name=youNameFile, group=yourGroupFIle, start=yourStart, end=yourEnd, maxambig=yourMaxambig,  \n");
		m->mothurOut("maxhomop=yourMaxhomop, minlength=youMinlength, maxlength=yourMaxlength)  \n");	
		m->mothurOut("Example screen.seqs(fasta=abrecovery.fasta, name=abrecovery.names, group=abrecovery.groups, start=..., end=..., maxambig=..., maxhomop=..., minlength=..., maxlength=...).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ScreenSeqsCommand::~ScreenSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ScreenSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
				
		string goodSeqFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "good" + m->getExtension(fastafile);
		string badAccnosFile =  outputDir + m->getRootName(m->getSimpleName(fastafile)) + "bad.accnos";
		
		int numFastaSeqs = 0;
		set<string> badSeqNames;
		int start = time(NULL);
		
#ifdef USE_MPI	
			int pid, end, numSeqsPerProcessor; 
			int tag = 2001;
			vector<unsigned long int> MPIPos;
			
			MPI_Status status; 
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
			MPI_Comm_size(MPI_COMM_WORLD, &processors); 

			MPI_File inMPI;
			MPI_File outMPIGood;
			MPI_File outMPIBadAccnos;
			
			int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
			int inMode=MPI_MODE_RDONLY; 
			
			char outGoodFilename[1024];
			strcpy(outGoodFilename, goodSeqFile.c_str());

			char outBadAccnosFilename[1024];
			strcpy(outBadAccnosFilename, badAccnosFile.c_str());

			char inFileName[1024];
			strcpy(inFileName, fastafile.c_str());
			
			MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
			MPI_File_open(MPI_COMM_WORLD, outGoodFilename, outMode, MPI_INFO_NULL, &outMPIGood);
			MPI_File_open(MPI_COMM_WORLD, outBadAccnosFilename, outMode, MPI_INFO_NULL, &outMPIBadAccnos);
			
			if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIGood); MPI_File_close(&outMPIBadAccnos); return 0; }
			
			if (pid == 0) { //you are the root process 
				
				MPIPos = m->setFilePosFasta(fastafile, numFastaSeqs); //fills MPIPos, returns numSeqs
				
				//send file positions to all processes
				for(int i = 1; i < processors; i++) { 
					MPI_Send(&numFastaSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
				}
				
				//figure out how many sequences you have to align
				numSeqsPerProcessor = numFastaSeqs / processors;
				int startIndex =  pid * numSeqsPerProcessor;
				if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
				
				//align your part
				driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPIGood, outMPIBadAccnos, MPIPos, badSeqNames);

				if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIGood);  MPI_File_close(&outMPIBadAccnos);  return 0; }

				for (int i = 1; i < processors; i++) {
				
					//get bad lists
					int badSize;
					MPI_Recv(&badSize, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
					/*for (int j = 0; j < badSize; j++) {
						int length;
						MPI_Recv(&length, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);  //recv the length of the name
						char* buf2 = new char[length];										//make space to recieve it
						MPI_Recv(buf2, length, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);  //get name
						
						string tempBuf = buf2;
						if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
						delete buf2;
						
						badSeqNames.insert(tempBuf);
					}*/
				}
			}else{ //you are a child process
				MPI_Recv(&numFastaSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPIPos.resize(numFastaSeqs+1);
				MPI_Recv(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);

				//figure out how many sequences you have to align
				numSeqsPerProcessor = numFastaSeqs / processors;
				int startIndex =  pid * numSeqsPerProcessor;
				if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
				
				//align your part
				driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPIGood, outMPIBadAccnos, MPIPos, badSeqNames);

				if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIGood);  MPI_File_close(&outMPIBadAccnos); return 0; }
				
				//send bad list	
				int badSize = badSeqNames.size();
				MPI_Send(&badSize, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				
				/*
				set<string>::iterator it;
				for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {
					string name = *it;
					int length = name.length();
					char* buf2 = new char[length];
					memcpy(buf2, name.c_str(), length);
					
					MPI_Send(&length, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
					MPI_Send(buf2, length, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
				}*/
			}
			
			//close files 
			MPI_File_close(&inMPI);
			MPI_File_close(&outMPIGood);
			MPI_File_close(&outMPIBadAccnos);
			MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
					
#else
			vector<unsigned long int> positions = m->divideFile(fastafile, processors);
				
			for (int i = 0; i < (positions.size()-1); i++) {
				lines.push_back(new linePair(positions[i], positions[(i+1)]));
			}	
						
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				numFastaSeqs = driver(lines[0], goodSeqFile, badAccnosFile, fastafile, badSeqNames);
				
				if (m->control_pressed) { remove(goodSeqFile.c_str()); return 0; }
				
			}else{
				processIDS.resize(0);
				
				numFastaSeqs = createProcesses(goodSeqFile, badAccnosFile, fastafile, badSeqNames); 
				
				rename((goodSeqFile + toString(processIDS[0]) + ".temp").c_str(), goodSeqFile.c_str());
				rename((badAccnosFile + toString(processIDS[0]) + ".temp").c_str(), badAccnosFile.c_str());
				
				//append alignment and report files
				for(int i=1;i<processors;i++){
					m->appendFiles((goodSeqFile + toString(processIDS[i]) + ".temp"), goodSeqFile);
					remove((goodSeqFile + toString(processIDS[i]) + ".temp").c_str());
			
					m->appendFiles((badAccnosFile + toString(processIDS[i]) + ".temp"), badAccnosFile);
					remove((badAccnosFile + toString(processIDS[i]) + ".temp").c_str());
				}
				
				if (m->control_pressed) { remove(goodSeqFile.c_str()); return 0; }
				
				//read badSeqs in because root process doesnt know what other "bad" seqs the children found
				ifstream inBad;
				int ableToOpen = m->openInputFile(badAccnosFile, inBad, "no error");
				
				if (ableToOpen == 0) {
					badSeqNames.clear();
					string tempName;
					while (!inBad.eof()) {
						inBad >> tempName; m->gobble(inBad);
						badSeqNames.insert(tempName);
					}
					inBad.close();
				}
			}
	#else
			numFastaSeqs = driver(lines[0], goodSeqFile, badAccnosFile, fastafile, badSeqNames);
			
			if (m->control_pressed) { remove(goodSeqFile.c_str()); return 0; }
			
	#endif

#endif		

		#ifdef USE_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should fix files
			
				//read accnos file with all names in it, process 0 just has its names
				MPI_File inMPIAccnos;
				MPI_Offset size;
			
				char inFileName[1024];
				strcpy(inFileName, badAccnosFile.c_str());
			
				MPI_File_open(MPI_COMM_SELF, inFileName, inMode, MPI_INFO_NULL, &inMPIAccnos);  //comm, filename, mode, info, filepointer
				MPI_File_get_size(inMPIAccnos, &size);
			
				char* buffer = new char[size];
				MPI_File_read(inMPIAccnos, buffer, size, MPI_CHAR, &status);
			
				string tempBuf = buffer;
				if (tempBuf.length() > size) { tempBuf = tempBuf.substr(0, size);  }
				istringstream iss (tempBuf,istringstream::in);

				delete buffer;
				MPI_File_close(&inMPIAccnos);
				
				badSeqNames.clear();
				string tempName;
				while (!iss.eof()) {
					iss >> tempName; m->gobble(iss);
					badSeqNames.insert(tempName);
				}
		#endif
																					
		if(namefile != "" && groupfile != "")	{	
			screenNameGroupFile(badSeqNames);	
			if (m->control_pressed) {  remove(goodSeqFile.c_str()); return 0; }
		}else if(namefile != "")	{	
			screenNameGroupFile(badSeqNames);
			if (m->control_pressed) {  remove(goodSeqFile.c_str());  return 0; }	
		}else if(groupfile != "")				{	screenGroupFile(badSeqNames);		}	// this screens just the group
		
		if (m->control_pressed) { remove(goodSeqFile.c_str());  return 0; }

		if(alignreport != "")					{	screenAlignReport(badSeqNames);		}
		
		if (m->control_pressed) { remove(goodSeqFile.c_str());  return 0; }
		
		#ifdef USE_MPI
			}
		#endif

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(goodSeqFile); m->mothurOutEndLine();	
		m->mothurOut(badAccnosFile); m->mothurOutEndLine();	
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
		m->mothurOutEndLine();

		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to screen " + toString(numFastaSeqs) + " sequences.");
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************

int ScreenSeqsCommand::screenNameGroupFile(set<string> badSeqNames){
	try {
		ifstream inputNames;
		m->openInputFile(namefile, inputNames);
		set<string> badSeqGroups;
		string seqName, seqList, group;
		set<string>::iterator it;

		string goodNameFile = outputDir + m->getRootName(m->getSimpleName(namefile)) + "good" + m->getExtension(namefile);
		outputNames.push_back(goodNameFile); 
		
		ofstream goodNameOut;	m->openOutputFile(goodNameFile, goodNameOut);
	
		while(!inputNames.eof()){
			if (m->control_pressed) { goodNameOut.close();  inputNames.close(); remove(goodNameFile.c_str());  return 0; }

			inputNames >> seqName >> seqList;
			it = badSeqNames.find(seqName);
			
			if(it != badSeqNames.end()){
				badSeqNames.erase(it);
				
				if(namefile != ""){
					int start = 0;
					for(int i=0;i<seqList.length();i++){
						if(seqList[i] == ','){
							badSeqGroups.insert(seqList.substr(start,i-start));
							start = i+1;
						}					
					}
					badSeqGroups.insert(seqList.substr(start,seqList.length()-start));
				}
			}
			else{
				goodNameOut << seqName << '\t' << seqList << endl;
			}
			m->gobble(inputNames);
		}
		inputNames.close();
		goodNameOut.close();
	
		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your namefile does not include the sequence " + *it + " please correct."); 
				m->mothurOutEndLine();
			}
		}

		if(groupfile != ""){
			
			ifstream inputGroups;
			m->openInputFile(groupfile, inputGroups);

			string goodGroupFile = outputDir + m->getRootName(m->getSimpleName(groupfile)) + "good" + m->getExtension(groupfile);
			outputNames.push_back(goodGroupFile);  
			
			ofstream goodGroupOut;	m->openOutputFile(goodGroupFile, goodGroupOut);
			
			while(!inputGroups.eof()){
				if (m->control_pressed) { goodGroupOut.close(); inputGroups.close(); remove(goodNameFile.c_str());  remove(goodGroupFile.c_str()); return 0; }

				inputGroups >> seqName >> group;

				it = badSeqGroups.find(seqName);
				
				if(it != badSeqGroups.end()){
					badSeqGroups.erase(it);
				}
				else{
					goodGroupOut << seqName << '\t' << group << endl;
				}
				m->gobble(inputGroups);
			}
			inputGroups.close();
			goodGroupOut.close();
			
			//we were unable to remove some of the bad sequences
			if (badSeqGroups.size() != 0) {
				for (it = badSeqGroups.begin(); it != badSeqGroups.end(); it++) {  
					m->mothurOut("Your groupfile does not include the sequence " + *it + " please correct."); 
					m->mothurOutEndLine();
				}
			}
		}
			
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenNameGroupFile");
		exit(1);
	}
}

//***************************************************************************************************************

int ScreenSeqsCommand::screenGroupFile(set<string> badSeqNames){
	try {
		ifstream inputGroups;
		m->openInputFile(groupfile, inputGroups);
		string seqName, group;
		set<string>::iterator it;
		
		string goodGroupFile = outputDir + m->getRootName(m->getSimpleName(groupfile)) + "good" + m->getExtension(groupfile);
		outputNames.push_back(goodGroupFile); 
		ofstream goodGroupOut;	m->openOutputFile(goodGroupFile, goodGroupOut);
		
		while(!inputGroups.eof()){
			if (m->control_pressed) { goodGroupOut.close(); inputGroups.close(); remove(goodGroupFile.c_str()); return 0; }

			inputGroups >> seqName >> group;
			it = badSeqNames.find(seqName);
			
			if(it != badSeqNames.end()){
				badSeqNames.erase(it);
			}
			else{
				goodGroupOut << seqName << '\t' << group << endl;
			}
			m->gobble(inputGroups);
		}
		
		if (m->control_pressed) { goodGroupOut.close();  inputGroups.close(); remove(goodGroupFile.c_str());  return 0; }

		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your groupfile does not include the sequence " + *it + " please correct."); 
				m->mothurOutEndLine();
			}
		}
		
		inputGroups.close();
		goodGroupOut.close();
		
		if (m->control_pressed) { remove(goodGroupFile.c_str());   }
		
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenGroupFile");
		exit(1);
	}
}

//***************************************************************************************************************

int ScreenSeqsCommand::screenAlignReport(set<string> badSeqNames){
	try {
		ifstream inputAlignReport;
		m->openInputFile(alignreport, inputAlignReport);
		string seqName, group;
		set<string>::iterator it;
		
		string goodAlignReportFile = outputDir + m->getRootName(m->getSimpleName(alignreport)) + "good" + m->getExtension(alignreport);
		outputNames.push_back(goodAlignReportFile);  
		ofstream goodAlignReportOut;	m->openOutputFile(goodAlignReportFile, goodAlignReportOut);

		while (!inputAlignReport.eof())	{		//	need to copy header
			char c = inputAlignReport.get();
			goodAlignReportOut << c;
			if (c == 10 || c == 13){	break;	}	
		}

		while(!inputAlignReport.eof()){
			if (m->control_pressed) { goodAlignReportOut.close(); inputAlignReport.close(); remove(goodAlignReportFile.c_str()); return 0; }

			inputAlignReport >> seqName;
			it = badSeqNames.find(seqName);
			string line;		
			while (!inputAlignReport.eof())	{		//	need to copy header
				char c = inputAlignReport.get();
				line += c;
				if (c == 10 || c == 13){	break;	}	
			}
			
			if(it != badSeqNames.end()){
				badSeqNames.erase(it);
			}
			else{
				goodAlignReportOut << seqName << '\t' << line;
			}
			m->gobble(inputAlignReport);
		}
		
		if (m->control_pressed) { goodAlignReportOut.close();  inputAlignReport.close(); remove(goodAlignReportFile.c_str());  return 0; }

		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your alignreport file does not include the sequence " + *it + " please correct."); 
				m->mothurOutEndLine();
			}
		}

		inputAlignReport.close();
		goodAlignReportOut.close();
				
		if (m->control_pressed) {  remove(goodAlignReportFile.c_str());  return 0; }
		
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenAlignReport");
		exit(1);
	}
	
}
//**********************************************************************************************************************

int ScreenSeqsCommand::driver(linePair* filePos, string goodFName, string badAccnosFName, string filename, set<string>& badSeqNames){
	try {
		ofstream goodFile;
		m->openOutputFile(goodFName, goodFile);
		
		ofstream badAccnosFile;
		m->openOutputFile(badAccnosFName, badAccnosFile);
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
	
		while (!done) {
		
			if (m->control_pressed) {  return 0; }
			
			Sequence currSeq(inFASTA); m->gobble(inFASTA);
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
				if(goodSeq == 1 && startPos != -1 && startPos < currSeq.getStartPos())			{	goodSeq = 0;	}
				if(goodSeq == 1 && endPos != -1 && endPos > currSeq.getEndPos())				{	goodSeq = 0;	}
				if(goodSeq == 1 && maxAmbig != -1 && maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	}
				if(goodSeq == 1 && maxHomoP != -1 && maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	}
				if(goodSeq == 1 && minLength != -1 && minLength > currSeq.getNumBases())		{	goodSeq = 0;	}
				if(goodSeq == 1 && maxLength != -1 && maxLength < currSeq.getNumBases())		{	goodSeq = 0;	}
				
				if(goodSeq == 1){
					currSeq.printSequence(goodFile);	
				}
				else{
					badAccnosFile << currSeq.getName() << endl;
					badSeqNames.insert(currSeq.getName());
				}
			count++;
			}
			
			unsigned long int pos = inFASTA.tellg();
			if ((pos == -1) || (pos >= filePos->end)) { break; }
			
			//report progress
			if((count) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		
			
		goodFile.close();
		inFASTA.close();
		badAccnosFile.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ScreenSeqsCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& goodFile, MPI_File& badAccnosFile, vector<unsigned long int>& MPIPos, set<string>& badSeqNames){
	try {
		string outputString = "";
		MPI_Status statusGood; 
		MPI_Status statusBadAccnos; 
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

		for(int i=0;i<num;i++){
		
			if (m->control_pressed) {  return 0; }
			
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];

			char* buf4 = new char[length];
			memcpy(buf4, outputString.c_str(), length);

			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;	delete buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);
			
			Sequence currSeq(iss);			
			
			//process seq
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
				if(goodSeq == 1 && startPos != -1 && startPos < currSeq.getStartPos())			{	goodSeq = 0;	}
				if(goodSeq == 1 && endPos != -1 && endPos > currSeq.getEndPos())				{	goodSeq = 0;	}
				if(goodSeq == 1 && maxAmbig != -1 && maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	}
				if(goodSeq == 1 && maxHomoP != -1 && maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	}
				if(goodSeq == 1 && minLength != -1 && minLength > currSeq.getNumBases())		{	goodSeq = 0;	}
				if(goodSeq == 1 && maxLength != -1 && maxLength < currSeq.getNumBases())		{	goodSeq = 0;	}
				
				if(goodSeq == 1){
					outputString =  ">" + currSeq.getName() + "\n" + currSeq.getAligned() + "\n";
				
					//print good seq
					length = outputString.length();
					char* buf2 = new char[length];
					memcpy(buf2, outputString.c_str(), length);
					
					MPI_File_write_shared(goodFile, buf2, length, MPI_CHAR, &statusGood);
					delete buf2;
				}
				else{

					badSeqNames.insert(currSeq.getName());
					
					//write to bad accnos file
					outputString = currSeq.getName() + "\n";
				
					length = outputString.length();
					char* buf3 = new char[length];
					memcpy(buf3, outputString.c_str(), length);
					
					MPI_File_write_shared(badAccnosFile, buf3, length, MPI_CHAR, &statusBadAccnos);
					delete buf3;
				}
			}
		}
				
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driverMPI");
		exit(1);
	}
}
#endif
/**************************************************************************************************/

int ScreenSeqsCommand::createProcesses(string goodFileName, string badAccnos, string filename, set<string>& badSeqNames) {
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
				num = driver(lines[process], goodFileName + toString(getpid()) + ".temp", badAccnos + toString(getpid()) + ".temp", filename, badSeqNames);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename + toString(getpid()) + ".num.temp";
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
			string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); remove(tempFile.c_str());
		}
		
		return num;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "createProcesses");
		exit(1);
	}
}

//***************************************************************************************************************


