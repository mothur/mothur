/*
 *  distancecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "distancecommand.h"

//**********************************************************************************************************************
vector<string> DistanceCommand::setParameters(){	
	try {
		CommandParameter pcolumn("column", "InputTypes", "", "", "none", "none", "OldFastaColumn",false,false); parameters.push_back(pcolumn);
		CommandParameter poldfasta("oldfasta", "InputTypes", "", "", "none", "none", "OldFastaColumn",false,false); parameters.push_back(poldfasta);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter poutput("output", "Multiple", "column-lt-square", "column", "", "", "",false,false); parameters.push_back(poutput);
		CommandParameter pcalc("calc", "Multiple", "nogaps-eachgap-onegap", "onegap", "", "", "",false,false); parameters.push_back(pcalc);
		CommandParameter pcountends("countends", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pcountends);
		CommandParameter pcompress("compress", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pcompress);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pcutoff("cutoff", "Number", "", "1.0", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DistanceCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The dist.seqs command reads a file containing sequences and creates a distance file.\n";
		helpString += "The dist.seqs command parameters are fasta, oldfasta, column, calc, countends, output, compress, cutoff and processors.  \n";
		helpString += "The fasta parameter is required, unless you have a valid current fasta file.\n";
		helpString += "The oldfasta and column parameters allow you to append the distances calculated to the column file.\n";
		helpString += "The calc parameter allows you to specify the method of calculating the distances.  Your options are: nogaps, onegap or eachgap. The default is onegap.\n";
		helpString += "The countends parameter allows you to specify whether to include terminal gaps in distance.  Your options are: T or F. The default is T.\n";
		helpString += "The cutoff parameter allows you to specify maximum distance to keep. The default is 1.0.\n";
		helpString += "The output parameter allows you to specify format of your distance matrix. Options are column, lt, and square. The default is column.\n";
		helpString += "The processors parameter allows you to specify number of processors to use.  The default is 1.\n";
		helpString += "The compress parameter allows you to indicate that you want the resulting distance file compressed.  The default is false.\n";
		helpString += "The dist.seqs command should be in the following format: \n";
		helpString += "dist.seqs(fasta=yourFastaFile, calc=yourCalc, countends=yourEnds, cutoff= yourCutOff, processors=yourProcessors) \n";
		helpString += "Example dist.seqs(fasta=amazon.fasta, calc=eachgap, countends=F, cutoff= 2.0, processors=3).\n";
		helpString += "Note: No spaces between parameter labels (i.e. calc), '=' and parameters (i.e.yourCalc).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
DistanceCommand::DistanceCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["phylip"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "DistanceCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
DistanceCommand::DistanceCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		Estimators.clear();
				
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter("dist.seqs");
			map<string, string>::iterator it2;
		
			//check to make sure all parameters are valid for command
			for (it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["phylip"] = tempOutNames;
			outputTypes["column"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it2 = parameters.find("fasta");
				//user has given a template file
				if(it2 != parameters.end()){ 
					path = m->hasPath(it2->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it2->second;		}
				}
				
				it2 = parameters.find("oldfasta");
				//user has given a template file
				if(it2 != parameters.end()){ 
					path = m->hasPath(it2->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oldfasta"] = inputDir + it2->second;		}
				}
				
				it2 = parameters.find("column");
				//user has given a template file
				if(it2 != parameters.end()){ 
					path = m->hasPath(it2->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it2->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); 
					ifstream inFASTA;
					m->openInputFile(fastafile, inFASTA);
					alignDB = SequenceDB(inFASTA); 
					inFASTA.close();
				}else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else if (fastafile == "not open") { abort = true; }	
			else{
				ifstream inFASTA;
				m->openInputFile(fastafile, inFASTA);
				alignDB = SequenceDB(inFASTA); 
				inFASTA.close();
				m->setFastaFile(fastafile);
			}
			
			oldfastafile = validParameter.validFile(parameters, "oldfasta", true);
			if (oldfastafile == "not found") { oldfastafile = ""; }
			else if (oldfastafile == "not open") { abort = true; }	
			
			column = validParameter.validFile(parameters, "column", true);
			if (column == "not found") { column = ""; }
			else if (column == "not open") { abort = true; }	
			else { m->setColumnFile(column); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "onegap";  }
			else { 
				 if (calc == "default")  {  calc = "onegap";  }
			}
			m->splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "countends", false);	if(temp == "not found"){	temp = "T";	}
			convert(temp, countends); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if(temp == "not found"){	temp = "1.0"; }
			convert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "compress", false);		if(temp == "not found"){  temp = "F"; }
			convert(temp, compress);

			output = validParameter.validFile(parameters, "output", false);		if(output == "not found"){	output = "column"; }
			
			if (((column != "") && (oldfastafile == "")) || ((column == "") && (oldfastafile != ""))) { m->mothurOut("If you provide column or oldfasta, you must provide both."); m->mothurOutEndLine(); abort=true; }
			
			if ((column != "") && (oldfastafile != "") && (output != "column")) { m->mothurOut("You have provided column and oldfasta, indicating you want to append distances to your column file. Your output must be in column format to do so."); m->mothurOutEndLine(); abort=true; }
			
			if ((output != "column") && (output != "lt") && (output != "square")) { m->mothurOut(output + " is not a valid output form. Options are column, lt and square. I will use column."); m->mothurOutEndLine(); output = "column"; }

		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "DistanceCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int DistanceCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		int startTime = time(NULL);
		
		//save number of new sequence
		numNewFasta = alignDB.getNumSeqs();
		
		//sanity check the oldfasta and column file as well as add oldfasta sequences to alignDB
		if ((oldfastafile != "") && (column != ""))  {	if (!(sanityCheck())) { return 0; }  }
		
		if (m->control_pressed) { return 0; }
		
		int numSeqs = alignDB.getNumSeqs();
		cutoff += 0.005;
		
		if (!alignDB.sameLength()) {  m->mothurOut("[ERROR]: your sequences are not the same length, aborting."); m->mothurOutEndLine(); return 0; }
		
		string outputFile;
				
		if (output == "lt") { //does the user want lower triangle phylip formatted file 
			outputFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "phylip.dist";
			m->mothurRemove(outputFile); outputTypes["phylip"].push_back(outputFile);
			
			//output numSeqs to phylip formatted dist file
		}else if (output == "column") { //user wants column format
			outputFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "dist";
			outputTypes["column"].push_back(outputFile);
			
			//so we don't accidentally overwrite
			if (outputFile == column) { 
				string tempcolumn = column + ".old"; 
				rename(column.c_str(), tempcolumn.c_str());
			}
			
			m->mothurRemove(outputFile);
		}else { //assume square
			outputFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "square.dist";
			m->mothurRemove(outputFile);
			outputTypes["phylip"].push_back(outputFile);
		}
		

#ifdef USE_MPI
		
		int pid, start, end; 
		int tag = 2001;
				
		MPI_Status status; 
		MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		//each process gets where it should start and stop in the file
		if (output != "square") {
			start = int (sqrt(float(pid)/float(processors)) * numSeqs);
			end = int (sqrt(float(pid+1)/float(processors)) * numSeqs);
		}else{
			start = int ((float(pid)/float(processors)) * numSeqs);
			end = int ((float(pid+1)/float(processors)) * numSeqs);
		}
		
		if (output == "column") {
			MPI_File outMPI;
			int amode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 

			//char* filename = new char[outputFile.length()];
			//memcpy(filename, outputFile.c_str(), outputFile.length());
			
			char filename[1024];
			strcpy(filename, outputFile.c_str());
			
			MPI_File_open(MPI_COMM_WORLD, filename, amode, MPI_INFO_NULL, &outMPI);
			//delete filename;

			if (pid == 0) { //you are the root process 
				
				//do your part
				string outputMyPart;
				
				driverMPI(start, end, outMPI, cutoff); 
				
				if (m->control_pressed) { outputTypes.clear(); MPI_File_close(&outMPI);   return 0; }
			
				//wait on chidren
				for(int i = 1; i < processors; i++) { 
					if (m->control_pressed) { outputTypes.clear();  MPI_File_close(&outMPI);    return 0; }
					
					char buf[5];
					MPI_Recv(buf, 5, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
				}
			}else { //you are a child process
				//do your part
				driverMPI(start, end, outMPI, cutoff); 
				
				if (m->control_pressed) { outputTypes.clear();  MPI_File_close(&outMPI);   return 0; }
			
				char buf[5];
				strcpy(buf, "done"); 
				//tell parent you are done.
				MPI_Send(buf, 5, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
			}
			
			MPI_File_close(&outMPI);
			
		}else { //lower triangle format
			if (pid == 0) { //you are the root process 
			
				//do your part
				string outputMyPart;
				unsigned long long mySize;
				
				if (output != "square"){ driverMPI(start, end, outputFile, mySize); }
				else { driverMPI(start, end, outputFile, mySize, output); }
	
				if (m->control_pressed) {  outputTypes.clear();   return 0; }
				
				int amode=MPI_MODE_APPEND|MPI_MODE_WRONLY|MPI_MODE_CREATE; //
				MPI_File outMPI;
				MPI_File inMPI;

				//char* filename = new char[outputFile.length()];
				//memcpy(filename, outputFile.c_str(), outputFile.length());
				
				char filename[1024];
				strcpy(filename, outputFile.c_str());

				MPI_File_open(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &outMPI);
				//delete filename;

				//wait on chidren
				for(int b = 1; b < processors; b++) { 
					unsigned long long fileSize;
					
					if (m->control_pressed) { outputTypes.clear();  MPI_File_close(&outMPI);   return 0; }
					
					MPI_Recv(&fileSize, 1, MPI_LONG, b, tag, MPI_COMM_WORLD, &status); 
					
					string outTemp = outputFile + toString(b) + ".temp";

					char* buf = new char[outTemp.length()];
					memcpy(buf, outTemp.c_str(), outTemp.length());
					
					MPI_File_open(MPI_COMM_SELF, buf, MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);
					delete buf;

					int count = 0;
					while (count < fileSize) { 
						char buf2[1];
						MPI_File_read(inMPI, buf2, 1, MPI_CHAR, &status);
						MPI_File_write(outMPI, buf2, 1, MPI_CHAR, &status);
						count += 1;
					}
					
					MPI_File_close(&inMPI); //deleted on close
				}
				
				MPI_File_close(&outMPI);
			}else { //you are a child process
				//do your part
				unsigned long long size;
				if (output != "square"){ driverMPI(start, end, (outputFile + toString(pid) + ".temp"), size); }
				else { driverMPI(start, end, (outputFile + toString(pid) + ".temp"), size, output); }
				
				if (m->control_pressed) {  return 0; }
			
				//tell parent you are done.
				MPI_Send(&size, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#else		
				
	//#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		//if you don't need to fork anything
		if(processors == 1){
			if (output != "square") {  driver(0, numSeqs, outputFile, cutoff); }
			else { driver(0, numSeqs, outputFile, "square");  }
		}else{ //you have multiple processors
			
			unsigned long long numDists = 0;
			
			if (output == "square") {
				 numDists = numSeqs * numSeqs;
			}else {
				for(int i=0;i<numSeqs;i++){
					for(int j=0;j<i;j++){
						numDists++;
						if (numDists > processors) { break; }
					}
				}
			}
			
			if (numDists < processors) { processors = numDists; }
			
			for (int i = 0; i < processors; i++) {
				distlinePair tempLine;
				lines.push_back(tempLine);
				if (output != "square") {
					lines[i].start = int (sqrt(float(i)/float(processors)) * numSeqs);
					lines[i].end = int (sqrt(float(i+1)/float(processors)) * numSeqs);
				}else{
					lines[i].start = int ((float(i)/float(processors)) * numSeqs);
					lines[i].end = int ((float(i+1)/float(processors)) * numSeqs);
				}
				
			}
			
			createProcesses(outputFile); 
		}
	//#else
		//ifstream inFASTA;
		//if (output != "square") {  driver(0, numSeqs, outputFile, cutoff); }
		//else { driver(0, numSeqs, outputFile, "square");  }
	//#endif
	
#endif
		if (m->control_pressed) { outputTypes.clear();  m->mothurRemove(outputFile); return 0; }
		
		#ifdef USE_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		//if (output == "square") {  convertMatrix(outputFile); }
		
		ifstream fileHandle;
		fileHandle.open(outputFile.c_str());
		if(fileHandle) {
			m->gobble(fileHandle);
			if (fileHandle.eof()) { m->mothurOut(outputFile + " is blank. This can result if there are no distances below your cutoff.");  m->mothurOutEndLine(); }
		}
		
		//append the old column file to the new one
		if ((oldfastafile != "") && (column != ""))  {
			//we had to rename the column file so we didnt overwrite above, but we want to keep old name
			if (outputFile == column) { 
				string tempcolumn = column + ".old";
				m->appendFiles(tempcolumn, outputFile);
				m->mothurRemove(tempcolumn);
			}else{
				m->appendFiles(outputFile, column);
				m->mothurRemove(outputFile);
				outputFile = column;
			}
			
			if (outputDir != "") { 
				string newOutputName = outputDir + m->getSimpleName(outputFile);
				rename(outputFile.c_str(), newOutputName.c_str());
				m->mothurRemove(outputFile);
				outputFile = newOutputName;
			}
		}

		
		#ifdef USE_MPI
			}
		#endif
		
		if (m->control_pressed) { outputTypes.clear();  m->mothurRemove(outputFile); return 0; }
		
		//set phylip file as new current phylipfile
		string current = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setPhylipFile(current); }
		}
		
		//set column file as new current columnfile
		itTypes = outputTypes.find("column");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setColumnFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFile); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - startTime) + " to calculate the distances for " + toString(numSeqs) + " sequences."); m->mothurOutEndLine();


		if (m->isTrue(compress)) {
			m->mothurOut("Compressing..."); m->mothurOutEndLine();
			m->mothurOut("(Replacing " + outputFile + " with " + outputFile + ".gz)"); m->mothurOutEndLine();
			system(("gzip -v " + outputFile).c_str());
			outputNames.push_back(outputFile + ".gz");
		}else { outputNames.push_back(outputFile); }

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
void DistanceCommand::createProcesses(string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				if (output != "square") {  driver(lines[process].start, lines[process].end, filename + toString(getpid()) + ".temp", cutoff); }
				else { driver(lines[process].start, lines[process].end, filename + toString(getpid()) + ".temp", "square"); }
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes. Error code: " + toString(pid)); m->mothurOutEndLine(); 
				perror(" : ");
				for (int i=0;i<processIDS.size();i++) {  int temp = processIDS[i]; kill (temp, SIGINT); }
				exit(0);
			}
		}
		
		//parent does its part
		if (output != "square") {  driver(lines[0].start, lines[0].end, filename, cutoff); }
		else { driver(lines[0].start, lines[0].end, filename, "square"); }
		
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the distanceData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//that's why the distance calculator was moved inside of the driver to make separate copies.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<distanceData*> pDataArray; //[processors-1];
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor-1 worker threads.
		for( int i=0; i<processors-1; i++ ){
			
			// Allocate memory for thread data.
			distanceData* tempDist = new distanceData(lines[i+1].start, lines[i+1].end, (filename + toString(i) + ".temp"), cutoff, alignDB, Estimators, m, output, numNewFasta, countends);
			pDataArray.push_back(tempDist);
			processIDS.push_back(i);
			
			//MyDistThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyDistThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
		//do your part
		if (output != "square") {  driver(lines[0].start, lines[0].end, filename, cutoff); }
		else { driver(lines[0].start, lines[0].end, filename, "square"); }
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif
		
		//append and remove temp files
		for (int i=0;i<processIDS.size();i++) { 
			m->appendFiles((filename + toString(processIDS[i]) + ".temp"), filename);
			m->mothurRemove((filename + toString(processIDS[i]) + ".temp"));
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driver(int startLine, int endLine, string dFileName, float cutoff){
	try {
		ValidCalculators validCalculator;
		Dist* distCalculator;
		if (m->isTrue(countends) == true) {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
					else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
					else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
				}
			}
		}else {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
					else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
					else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
				}
			}
		}
		
		int startTime = time(NULL);
		
		//column file
		ofstream outFile(dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		if((output == "lt") && startLine == 0){	outFile << alignDB.getNumSeqs() << endl;	}
		
		for(int i=startLine;i<endLine;i++){
			if(output == "lt")	{	
				string name = alignDB.get(i).getName();
				if (name.length() < 10) { //pad with spaces to make compatible
					while (name.length() < 10) {  name += " ";  }
				}
				outFile << name << '\t';	
			}
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) { delete distCalculator; outFile.close(); return 0;  }
				
				//if there was a column file given and we are appending, we don't want to calculate the distances that are already in the column file
				//the alignDB contains the new sequences and then the old, so if i an oldsequence and j is an old sequence then break out of this loop
				if ((i >= numNewFasta) && (j >= numNewFasta)) { break; }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					if (output == "column") { outFile << alignDB.get(i).getName() << ' ' << alignDB.get(j).getName() << ' ' << dist << endl; }
				}
				if (output == "lt") {  outFile << dist << '\t'; }
			}
			
			if (output == "lt") { outFile << endl; }
			
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
			}
			
		}
		m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		
		outFile.close();
		delete distCalculator;
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driver(int startLine, int endLine, string dFileName, string square){
	try {
		ValidCalculators validCalculator;
		Dist* distCalculator;
		if (m->isTrue(countends) == true) {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
					else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
					else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
				}
			}
		}else {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
					else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
					else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
				}
			}
		}
		
		int startTime = time(NULL);
		
		//column file
		ofstream outFile(dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		if(startLine == 0){	outFile << alignDB.getNumSeqs() << endl;	}
		
		for(int i=startLine;i<endLine;i++){
				
			string name = alignDB.get(i).getName();
			//pad with spaces to make compatible
			if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } }
				
			outFile << name << '\t';	
			
			for(int j=0;j<alignDB.getNumSeqs();j++){
				
				if (m->control_pressed) { delete distCalculator; outFile.close(); return 0;  }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				outFile << dist << '\t'; 
			}
			
			outFile << endl; 
			
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
			}
			
		}
		m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		
		outFile.close();
		delete distCalculator;
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driver");
		exit(1);
	}
}
#ifdef USE_MPI
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driverMPI(int startLine, int endLine, MPI_File& outMPI, float cutoff){
	try {
		
		ValidCalculators validCalculator;
		Dist* distCalculator;
		if (m->isTrue(countends) == true) {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
					else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
					else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
				}
			}
		}else {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
					else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
					else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
				}
			}
		}
		
		
		MPI_Status status;
		int startTime = time(NULL);
		
		string outputString = "";
		
		for(int i=startLine;i<endLine;i++){
	
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) {  delete distCalculator; return 0;  }
				
				//if there was a column file given and we are appending, we don't want to calculate the distances that are already in the column file
				//the alignDB contains the new sequences and then the old, so if i an oldsequence and j is an old sequence then break out of this loop
				if ((i >= numNewFasta) && (j >= numNewFasta)) { break; }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					 outputString += (alignDB.get(i).getName() + ' ' + alignDB.get(j).getName() + ' ' + toString(dist) + '\n'); 
				}
			}
			
			if(i % 100 == 0){
				//m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
				cout << i << '\t' << (time(NULL) - startTime) << endl;
			}
			
			 
			//send results to parent
			int length = outputString.length();

			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write_shared(outMPI, buf, length, MPI_CHAR, &status);
			outputString = "";
			delete buf;
			
		}
		
		//m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		cout << (endLine-1) << '\t' << (time(NULL) - startTime) << endl;	
		delete distCalculator;
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driverMPI");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driverMPI(int startLine, int endLine, string file, unsigned long long& size){
	try {
		ValidCalculators validCalculator;
		Dist* distCalculator;
		if (m->isTrue(countends) == true) {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
					else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
					else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
				}
			}
		}else {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
					else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
					else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
				}
			}
		}
		
		
		MPI_Status status;
		
		MPI_File outMPI;
		int amode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 

		//char* filename = new char[file.length()];
		//memcpy(filename, file.c_str(), file.length());
		
		char filename[1024];
		strcpy(filename, file.c_str());

		MPI_File_open(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &outMPI);
		//delete filename;

		int startTime = time(NULL);
		
		string outputString = "";
		size = 0;
		
		if(startLine == 0){	outputString += toString(alignDB.getNumSeqs()) + "\n";	}
		
		for(int i=startLine;i<endLine;i++){
				
			string name = alignDB.get(i).getName();
			if (name.length() < 10) { //pad with spaces to make compatible
				while (name.length() < 10) {  name += " ";  }
			}
			outputString += name + "\t";	
			
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) { delete distCalculator; return 0;  }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				outputString += toString(dist) + "\t"; 
			}
			
			outputString += "\n"; 

		
			if(i % 100 == 0){
				//m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
				cout << i << '\t' << (time(NULL) - startTime) << endl;
			}
			
			
			//send results to parent
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write(outMPI, buf, length, MPI_CHAR, &status);
			size += outputString.length();
			outputString = "";
			delete buf;
		}
		
		//m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		cout << (endLine-1) << '\t' << (time(NULL) - startTime) << endl;
		MPI_File_close(&outMPI);
		delete distCalculator;
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driverMPI");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driverMPI(int startLine, int endLine, string file, unsigned long long& size, string square){
	try {
		ValidCalculators validCalculator;
		Dist* distCalculator;
		if (m->isTrue(countends) == true) {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
					else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
					else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
				}
			}
		}else {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
					if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
					else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
					else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
				}
			}
		}
		
		MPI_Status status;
		
		MPI_File outMPI;
		int amode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 

		//char* filename = new char[file.length()];
		//memcpy(filename, file.c_str(), file.length());
		
		char filename[1024];
		strcpy(filename, file.c_str());

		MPI_File_open(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &outMPI);
		//delete filename;

		int startTime = time(NULL);
		
		string outputString = "";
		size = 0;
		
		if(startLine == 0){	outputString += toString(alignDB.getNumSeqs()) + "\n";	}
		
		for(int i=startLine;i<endLine;i++){
				
			string name = alignDB.get(i).getName();
			if (name.length() < 10) { //pad with spaces to make compatible
				while (name.length() < 10) {  name += " ";  }
			}
			outputString += name + "\t";	
			
			for(int j=0;j<alignDB.getNumSeqs();j++){
				
				if (m->control_pressed) { delete distCalculator; return 0;  }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				outputString += toString(dist) + "\t"; 
			}
			
			outputString += "\n"; 

		
			if(i % 100 == 0){
				//m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
				cout << i << '\t' << (time(NULL) - startTime) << endl;
			}
			
			
			//send results to parent
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write(outMPI, buf, length, MPI_CHAR, &status);
			size += outputString.length();
			outputString = "";
			delete buf;
		}
		
		//m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		cout << (endLine-1) << '\t' << (time(NULL) - startTime) << endl;
		MPI_File_close(&outMPI);
		delete distCalculator;
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driverMPI");
		exit(1);
	}
}
#endif
/**************************************************************************************************
int DistanceCommand::convertMatrix(string outputFile) {
	try{

		//sort file by first column so the distances for each row are together
		string outfile = m->getRootName(outputFile) + "sorted.dist.temp";
		
		//use the unix sort 
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			string command = "sort -n " + outputFile + " -o " + outfile;
			system(command.c_str());
		#else //sort using windows sort
			string command = "sort " + outputFile + " /O " + outfile;
			system(command.c_str());
		#endif
		

		//output to new file distance for each row and save positions in file where new row begins
		ifstream in;
		m->openInputFile(outfile, in);
		
		ofstream out;
		m->openOutputFile(outputFile, out);
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);

		out << alignDB.getNumSeqs() << endl;
		
		//get first currentRow
		string first, currentRow, second;
		float dist;
		map<string, float> rowDists; //take advantage of the fact that maps are already sorted by key 
		map<string, float>::iterator it;
		
		in >> first;
		currentRow = first;
		
		rowDists[first] = 0.00; //distance to yourself is 0.0
		
		in.seekg(0);
		//m->openInputFile(outfile, in);
		
		while(!in.eof()) {
			if (m->control_pressed) { in.close(); m->mothurRemove(outfile); out.close(); return 0; }
			
			in >> first >> second >> dist; m->gobble(in);
				
			if (first != currentRow) {
				//print out last row
				out << currentRow << '\t'; //print name

				//print dists
				for (it = rowDists.begin(); it != rowDists.end(); it++) {
					out << it->second << '\t';
				}
				out << endl;
				
				//start new row
				currentRow = first;
				rowDists.clear();
				rowDists[first] = 0.00;
				rowDists[second] = dist;
			}else{
				rowDists[second] = dist;
			}
		}
		//print out last row
		out << currentRow << '\t'; //print name
				
		//print dists
		for (it = rowDists.begin(); it != rowDists.end(); it++) {
			out << it->second << '\t';
		}
		out << endl;
		
		in.close();
		out.close();
		
		m->mothurRemove(outfile);
		
		return 1;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "convertMatrix");
		exit(1);
	}
}
/**************************************************************************************************
int DistanceCommand::convertToLowerTriangle(string outputFile) {
	try{

		//sort file by first column so the distances for each row are together
		string outfile = m->getRootName(outputFile) + "sorted.dist.temp";
		
		//use the unix sort 
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			string command = "sort -n " + outputFile + " -o " + outfile;
			system(command.c_str());
		#else //sort using windows sort
			string command = "sort " + outputFile + " /O " + outfile;
			system(command.c_str());
		#endif
		

		//output to new file distance for each row and save positions in file where new row begins
		ifstream in;
		m->openInputFile(outfile, in);
		
		ofstream out;
		m->openOutputFile(outputFile, out);
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);

		out << alignDB.getNumSeqs() << endl;
		
		//get first currentRow
		string first, currentRow, second;
		float dist;
		int i, j;
		i = 0; j = 0;
		map<string, float> rowDists; //take advantage of the fact that maps are already sorted by key 
		map<string, float>::iterator it;
		
		in >> first;
		currentRow = first;
		
		rowDists[first] = 0.00; //distance to yourself is 0.0
		
		in.seekg(0);
		//m->openInputFile(outfile, in);
		
		while(!in.eof()) {
			if (m->control_pressed) { in.close(); m->mothurRemove(outfile); out.close(); return 0; }
			
			in >> first >> second >> dist; m->gobble(in);
				
			if (first != currentRow) {
				//print out last row
				out << currentRow << '\t'; //print name

				//print dists
				for (it = rowDists.begin(); it != rowDists.end(); it++) {
					if (j >= i) { break; }
					out << it->second << '\t';
					j++;
				}
				out << endl;
				
				//start new row
				currentRow = first;
				rowDists.clear();
				rowDists[first] = 0.00;
				rowDists[second] = dist;
				j = 0;
				i++;
			}else{
				rowDists[second] = dist;
			}
		}
		//print out last row
		out << currentRow << '\t'; //print name
				
		//print dists
		for (it = rowDists.begin(); it != rowDists.end(); it++) {
			out << it->second << '\t';
		}
		out << endl;
		
		in.close();
		out.close();
		
		m->mothurRemove(outfile);
		
		return 1;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "convertToLowerTriangle");
		exit(1);
	}
}
/**************************************************************************************************/
//its okay if the column file does not contain all the names in the fasta file, since some distance may have been above a cutoff,
//but no sequences can be in the column file that are not in oldfasta. also, if a distance is above the cutoff given then remove it.
//also check to make sure the 2 files have the same alignment length.
bool DistanceCommand::sanityCheck() {
	try{
		bool good = true;
		
		//make sure the 2 fasta files have the same alignment length
		ifstream in;
		m->openInputFile(fastafile, in);
		int fastaAlignLength = 0;
		if (in) { 
			Sequence tempIn(in);
			fastaAlignLength = tempIn.getAligned().length();
		}
		in.close();
		
		ifstream in2;
		m->openInputFile(oldfastafile, in2);
		int oldfastaAlignLength = 0;
		if (in2) { 
			Sequence tempIn2(in2);
			oldfastaAlignLength = tempIn2.getAligned().length();
		}
		in2.close();
		
		if (fastaAlignLength != oldfastaAlignLength) { m->mothurOut("fasta files do not have the same alignment length."); m->mothurOutEndLine(); return false;  }
		
		//read fasta file and save names as well as adding them to the alignDB
		set<string> namesOldFasta;
		
		ifstream inFasta;
		m->openInputFile(oldfastafile, inFasta);
		
		while (!inFasta.eof()) {
			if (m->control_pressed) {  inFasta.close(); return good;  }
		
			Sequence temp(inFasta);
			
			if (temp.getName() != "") {
				namesOldFasta.insert(temp.getName());  //save name
				alignDB.push_back(temp);  //add to DB
			}
			
			m->gobble(inFasta);
		}
		
		inFasta.close();
		
		//read through the column file checking names and removing distances above the cutoff
		ifstream inDist;
		m->openInputFile(column, inDist);
		
		ofstream outDist;
		string outputFile = column + ".temp";
		m->openOutputFile(outputFile, outDist);
		
		string name1, name2;
		float dist;
		while (!inDist.eof()) {
			if (m->control_pressed) {  inDist.close(); outDist.close(); m->mothurRemove(outputFile); return good;  }
		
			inDist >> name1 >> name2 >> dist; m->gobble(inDist);
			
			//both names are in fasta file and distance is below cutoff
			if ((namesOldFasta.count(name1) == 0) || (namesOldFasta.count(name2) == 0)) {  good = false; break;  }
			else{
				if (dist <= cutoff) {
					outDist << name1 << '\t' << name2 << '\t' << dist << endl;
				}
			}
		}
		
		inDist.close();
		outDist.close();
		
		if (good) {
			m->mothurRemove(column);
			rename(outputFile.c_str(), column.c_str());
		}else{
			m->mothurRemove(outputFile); //temp file is bad because file mismatch above
		}
		
		return good;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "sanityCheck");
		exit(1);
	}
}
/**************************************************************************************************/




