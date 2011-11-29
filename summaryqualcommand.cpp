/*
 *  summaryqualcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/28/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "summaryqualcommand.h"


//**********************************************************************************************************************
vector<string> SummaryQualCommand::setParameters(){	
	try {
		CommandParameter pqual("qfile", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pqual);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryQualCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The summary.qual command reads a quality file and an optional name file, and summarizes the quality information.\n";
		helpString += "The summary.tax command parameters are qfile, name and processors. qfile is required, unless you have a valid current quality file.\n";
		helpString += "The name parameter allows you to enter a name file associated with your quality file. \n";
		helpString += "The summary.qual command should be in the following format: \n";
		helpString += "summary.qual(qfile=yourQualityFile) \n";
		helpString += "Note: No spaces between parameter labels (i.e. qfile), '=' and parameters (i.e.yourQualityFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
SummaryQualCommand::SummaryQualCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "SummaryQualCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SummaryQualCommand::SummaryQualCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
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
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; }
			else if (qualfile == "not found") { 				
				qualfile = m->getQualFile(); 
				if (qualfile != "") { m->mothurOut("Using " + qualfile + " as input file for the qfile parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current quality file and the qfile parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setQualFile(qualfile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { m->setNameFile(namefile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(qualfile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "SummaryQualCommand");
		exit(1);
	}
}
//***************************************************************************************************************
int SummaryQualCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		int start = time(NULL);
		int numSeqs = 0;
		
		vector<int> position;
		vector<int> averageQ;
		vector< vector<int> > scores;
				
		if (m->control_pressed) { return 0; }
		
		if (namefile != "") { nameMap = m->readNames(namefile); }
		
		vector<unsigned long long> positions; 
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		positions = m->divideFile(qualfile, processors);
		for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else	
		if (processors == 1) {
			lines.push_back(linePair(0, 1000)); 
		}else {
			positions = m->setFilePosFasta(qualfile, numSeqs); 
			
			//figure out how many sequences you have to process
			int numSeqsPerProcessor = numSeqs / processors;
			for (int i = 0; i < processors; i++) {
				int startIndex =  i * numSeqsPerProcessor;
				if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
				lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
			}
		}
#endif
		
		
		if(processors == 1){ numSeqs = driverCreateSummary(position, averageQ, scores, qualfile, lines[0]);  }
		else{  numSeqs = createProcessesCreateSummary(position, averageQ, scores, qualfile);  }
		
		if (m->control_pressed) {  return 0; }
		
		//print summary file
		string summaryFile = outputDir + m->getRootName(m->getSimpleName(qualfile)) + "qual.summary";
		printQual(summaryFile, position, averageQ, scores);
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
		//output results to screen
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		m->mothurOutEndLine();
		m->mothurOut("Position\tNumSeqs\tAverageQ"); m->mothurOutEndLine();
		for (int i = 0; i < position.size(); i+=100) {
			float average = averageQ[i] / (float) position[i];
			cout << i << '\t' << position[i] << '\t' << average;
			m->mothurOutJustToLog(toString(i) + "\t" + toString(position[i]) + "\t" + toString(average)); m->mothurOutEndLine();
		}
		
		m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(summaryFile); m->mothurOutEndLine();	outputNames.push_back(summaryFile); outputTypes["summary"].push_back(summaryFile);
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
int SummaryQualCommand::driverCreateSummary(vector<int>& position, vector<int>& averageQ, vector< vector<int> >& scores, string filename, linePair filePos) {	
	try {
		ifstream in;
		m->openInputFile(filename, in);
		
		in.seekg(filePos.start);
		
		bool done = false;
		int count = 0;
		
		while (!done) {
			
			if (m->control_pressed) { in.close(); return 1; }
			
			QualityScores current(in); m->gobble(in);
			
			if (current.getName() != "") {
				
				int num = 1;
				if (namefile != "") {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(current.getName());
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				vector<int> thisScores = current.getQualityScores();
				
				//resize to num of positions setting number of seqs with that size to 1
				if (position.size() < thisScores.size()) { position.resize(thisScores.size(), 0); }
				if (averageQ.size() < thisScores.size()) { averageQ.resize(thisScores.size(), 0); }
				if (scores.size() < thisScores.size()) { 
					scores.resize(thisScores.size()); 
					for (int i = 0; i < scores.size(); i++) { scores[i].resize(41, 0); }
				}
				
				//increase counts of number of seqs with this position
				//average is really the total, we will average in execute
				for (int i = 0; i < thisScores.size(); i++) { 
					position[i] += num; 
					averageQ[i] += (thisScores[i] * num); //weighting for namesfile
					if (thisScores[i] > 40) { m->mothurOut("[ERROR]: " + current.getName() + " has a quality scores of " + toString(thisScores[i]) + ", expecting values to be less than 40."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { scores[i][thisScores[i]] += num; }  
				}
				
				count += num;
			}
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			unsigned long long pos = in.tellg();
			if ((pos == -1) || (pos >= filePos.end)) { break; }
#else
			if (in.eof()) { break; }
#endif
		}
		
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "driverCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/
int SummaryQualCommand::createProcessesCreateSummary(vector<int>& position, vector<int>& averageQ, vector< vector<int> >& scores, string filename) {
	try {
		int process = 1;
		int numSeqs = 0;
		processIDS.clear();
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				numSeqs = driverCreateSummary(position, averageQ, scores, qualfile, lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = qualfile + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				
				out << numSeqs << endl;
				out << position.size() << endl;
				for (int k = 0; k < position.size(); k++)			{		out << position[k] << '\t'; }  out << endl;
				for (int k = 0; k < averageQ.size(); k++)			{		out << averageQ[k] << '\t'; }  out << endl;
				for (int k = 0; k < scores.size(); k++)	{		
					for (int j = 0; j < 41; j++) {
						out << scores[k][j] << '\t'; 
					}
					out << endl;
				}  
				out << endl;
				
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do your part
		numSeqs = driverCreateSummary(position, averageQ, scores, qualfile, lines[0]);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = qualfile + toString(processIDS[i]) + ".num.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			int temp, tempNum;
			in >> tempNum; m->gobble(in); numSeqs += tempNum;
			in >> tempNum; m->gobble(in);
			
			if (position.size() < tempNum) { position.resize(tempNum, 0); }
			if (averageQ.size() < tempNum) { averageQ.resize(tempNum, 0); }
			if (scores.size() < tempNum) { 
				scores.resize(tempNum); 
				for (int i = 0; i < scores.size(); i++) { scores[i].resize(41, 0); }
			}
			
			for (int k = 0; k < tempNum; k++)			{		in >> temp; position[k]	+= temp;			}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; averageQ[k] += temp; 		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{	
				for (int j = 0; j < 41; j++) {
					in >> temp; scores[k][j] += temp;
					m->gobble(in);
				}	
			}
			
			in.close();
			m->mothurRemove(tempFilename);
		}
		
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumQualData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<seqSumQualData*> pDataArray; 
		DWORD   dwThreadIdArray[processors];
		HANDLE  hThreadArray[processors]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors; i++ ){
			
			// Allocate memory for thread data.
			seqSumQualData* tempSum = new seqSumQualData(&position, &averageQ, &scores, filename, m, lines[i].start, lines[i].end, namefile, nameMap);
			pDataArray.push_back(tempSum);
			processIDS.push_back(i);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MySeqSumQualThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			numSeqs += pDataArray[i]->count;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif		
		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "createProcessesCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/
int SummaryQualCommand::printQual(string sumFile, vector<int>& position, vector<int>& averageQ, vector< vector<int> >& scores) {
	try {
		ofstream out;
		m->openOutputFile(sumFile, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		outputNames.push_back(sumFile); outputTypes["summary"].push_back(sumFile);
		
		//print headings
		out << "Position\tnumSeqs\tAverageQ\t";
		for (int i = 0; i < 41; i++) { out << "q" << i << '\t'; }
		out << endl;
		
		for (int i = 0; i < position.size(); i++) {
			
			if (m->control_pressed) { out.close(); return 0; }
			
			float average = averageQ[i] / (float) position[i];
			out << i << '\t' << position[i] << '\t' << average << '\t';
			
			for (int j = 0; j < 41; j++) {
				out << scores[i][j] << '\t';
			}
			out << endl;
		}
		
		out.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryQualCommand", "printQual");
		exit(1);
	}
}

/**************************************************************************************/


