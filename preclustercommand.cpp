/*
 *  preclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "preclustercommand.h"
#include "deconvolutecommand.h"

//**********************************************************************************************************************
vector<string> PreClusterCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter pdiffs("diffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PreClusterCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The pre.cluster command groups sequences that are within a given number of base mismatches.\n";
		helpString += "The pre.cluster command outputs a new fasta and name file.\n";
		helpString += "The pre.cluster command parameters are fasta, names and diffs. The fasta parameter is required. \n";
		helpString += "The names parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n";
		helpString += "The group parameter allows you to provide a group file so you can cluster by group. \n";
		helpString += "The diffs parameter allows you to specify maximum number of mismatched bases allowed between sequences in a grouping. The default is 1.\n";
		helpString += "The pre.cluster command should be in the following format: \n";
		helpString += "pre.cluster(fasta=yourFastaFile, names=yourNamesFile, diffs=yourMaxDiffs) \n";
		helpString += "Example pre.cluster(fasta=amazon.fasta, diffs=2).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
PreClusterCommand::PreClusterCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["map"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

PreClusterCommand::PreClusterCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["map"] = tempOutNames;
		
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }	
			else { m->setFastaFile(fastafile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { namefile = ""; abort = true; }	
			else {  m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not found") { groupfile =  "";  bygroup = false; }
			else if (groupfile == "not open") { abort = true; groupfile =  ""; }	
			else {   m->setGroupFile(groupfile); bygroup = true;  }
			
			string temp	= validParameter.validFile(parameters, "diffs", false);		if(temp == "not found"){	temp = "1"; }
			m->mothurConvert(temp, diffs); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			if (namefile == "") {
				vector<string> files; files.push_back(fastafile);
				parser.getNameFile(files);
			}
		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int PreClusterCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		int start = time(NULL);
		
		string fileroot = outputDir + m->getRootName(m->getSimpleName(fastafile));
		string newFastaFile = fileroot + "precluster" + m->getExtension(fastafile);
		string newNamesFile = fileroot + "precluster.names";
		string newMapFile = fileroot + "precluster.map"; //add group name if by group
		outputNames.push_back(newFastaFile); outputTypes["fasta"].push_back(newFastaFile);
		outputNames.push_back(newNamesFile); outputTypes["name"].push_back(newNamesFile);
		
		
		if (bygroup) {
			//clear out old files
			ofstream outFasta; m->openOutputFile(newFastaFile, outFasta); outFasta.close();
			ofstream outNames; m->openOutputFile(newNamesFile, outNames);  outNames.close();
			newMapFile = fileroot + "precluster.";
			
			//parse fasta and name file by group
			SequenceParser* parser;
			if (namefile != "") { parser = new SequenceParser(groupfile, fastafile, namefile);	}
			else				{ parser = new SequenceParser(groupfile, fastafile);			}
			
			vector<string> groups = parser->getNamesOfGroups();
			
			if(processors == 1)	{	driverGroups(parser, newFastaFile, newNamesFile, newMapFile, 0, groups.size(), groups);	}
			else				{	createProcessesGroups(parser, newFastaFile, newNamesFile, newMapFile, groups);			}
			
			delete parser;
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}	 return 0; }
			
			//run unique.seqs for deconvolute results
			string inputString = "fasta=" + newFastaFile;
			if (namefile != "") { inputString += ", name=" + newNamesFile; }
			m->mothurOutEndLine(); 
			m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
			m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
			
			Command* uniqueCommand = new DeconvoluteCommand(inputString);
			uniqueCommand->execute();
			
			map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
			
			delete uniqueCommand;
			
			m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
			
			m->renameFile(filenames["fasta"][0], newFastaFile);
			
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run pre.cluster."); m->mothurOutEndLine(); 
				
		}else {
			if (namefile != "") { readNameFile(); }
		
			//reads fasta file and return number of seqs
			int numSeqs = readFASTA(); //fills alignSeqs and makes all seqs active
		
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
	
			if (numSeqs == 0) { m->mothurOut("Error reading fasta file...please correct."); m->mothurOutEndLine(); return 0;  }
			if (diffs > length) { m->mothurOut("Error: diffs is greater than your sequence length."); m->mothurOutEndLine(); return 0;  }
			
			int count = process(newMapFile);
			outputNames.push_back(newMapFile); outputTypes["map"].push_back(newMapFile);
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }	
			
			m->mothurOut("Total number of sequences before precluster was " + toString(alignSeqs.size()) + "."); m->mothurOutEndLine();
			m->mothurOut("pre.cluster removed " + toString(count) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine(); 
			printData(newFastaFile, newNamesFile);
			
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); 
		}
				
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}		
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::createProcessesGroups(SequenceParser* parser, string newFName, string newNName, string newMFile, vector<string> groups) {
	try {
		
		vector<int> processIDS;
		int process = 1;
		int num = 0;
		
		//sanity check
		if (groups.size() < processors) { processors = groups.size(); }
		
		//divide the groups between the processors
		vector<linePair> lines;
		int numGroupsPerProcessor = groups.size() / processors;
		for (int i = 0; i < processors; i++) {
			int startIndex =  i * numGroupsPerProcessor;
			int endIndex = (i+1) * numGroupsPerProcessor;
			if(i == (processors - 1)){	endIndex = groups.size(); 	}
			lines.push_back(linePair(startIndex, endIndex));
		}
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)		
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverGroups(parser, newFName + toString(getpid()) + ".temp", newNName + toString(getpid()) + ".temp", newMFile, lines[process].start, lines[process].end, groups);
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		num = driverGroups(parser, newFName, newNName, newMFile, lines[0].start, lines[0].end, groups);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
#else
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the preClusterData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<preClusterData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
			// Allocate memory for thread data.
			string extension = toString(i) + ".temp";
			
			preClusterData* tempPreCluster = new preClusterData(fastafile, namefile, groupfile, (newFName+extension), (newNName+extension), newMFile, groups, m, lines[i].start, lines[i].end, diffs, i);
			pDataArray.push_back(tempPreCluster);
			processIDS.push_back(i);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyPreclusterThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
				
		//using the main process as a worker saves time and memory
		num = driverGroups(parser, newFName, newNName, newMFile, lines[0].start, lines[0].end, groups);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			for (int j = 0; j < pDataArray[i]->mapFileNames.size(); j++) {
				outputNames.push_back(pDataArray[i]->mapFileNames[j]); outputTypes["map"].push_back(pDataArray[i]->mapFileNames[j]); 
			}
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
#endif		
		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			//newFName = m->getFullPathName(".\\" + newFName);
			//newNName = m->getFullPathName(".\\" + newNName);
			
			m->appendFiles((newFName + toString(processIDS[i]) + ".temp"), newFName);
			m->mothurRemove((newFName + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((newNName + toString(processIDS[i]) + ".temp"), newNName);
			m->mothurRemove((newNName + toString(processIDS[i]) + ".temp"));
		}
		
		return num;	
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "createProcessesGroups");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::driverGroups(SequenceParser* parser, string newFFile, string newNFile, string newMFile, int start, int end, vector<string> groups){
	try {
		
		int numSeqs = 0;
		
		//precluster each group
		for (int i = start; i < end; i++) {
			
			start = time(NULL);
			
			if (m->control_pressed) {  return 0; }
			
			m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[i] + ":"); m->mothurOutEndLine();
			
			map<string, string> thisNameMap;
			if (namefile != "") { thisNameMap = parser->getNameMap(groups[i]); }
			vector<Sequence> thisSeqs = parser->getSeqs(groups[i]);
			
			//fill alignSeqs with this groups info.
			numSeqs = loadSeqs(thisNameMap, thisSeqs);
			
			if (m->control_pressed) {   return 0; }
			
			if (diffs > length) { m->mothurOut("Error: diffs is greater than your sequence length."); m->mothurOutEndLine(); m->control_pressed = true; return 0;  }
			
			int count = process(newMFile+groups[i]+".map");
			outputNames.push_back(newMFile+groups[i]+".map"); outputTypes["map"].push_back(newMFile+groups[i]+".map");
			
			if (m->control_pressed) {  return 0; }
			
			m->mothurOut("Total number of sequences before pre.cluster was " + toString(alignSeqs.size()) + "."); m->mothurOutEndLine();
			m->mothurOut("pre.cluster removed " + toString(count) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine(); 
			printData(newFFile, newNFile);
			
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); 
			
		}
		
		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "driverGroups");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::process(string newMapFile){
	try {
		ofstream out;
		m->openOutputFile(newMapFile, out);
		
		//sort seqs by number of identical seqs
		sort(alignSeqs.begin(), alignSeqs.end(), comparePriority);
		
		int count = 0;
		int numSeqs = alignSeqs.size();
		
		//think about running through twice...
		for (int i = 0; i < numSeqs; i++) {
			
			//are you active
			//			itActive = active.find(alignSeqs[i].seq.getName());
			
			if (alignSeqs[i].active) {  //this sequence has not been merged yet
				
				string chunk = alignSeqs[i].seq.getName() + "\t" + toString(alignSeqs[i].numIdentical) + "\t" + toString(0) + "\t" + alignSeqs[i].seq.getAligned() + "\n";
				
				//try to merge it with all smaller seqs
				for (int j = i+1; j < numSeqs; j++) {
					
					if (m->control_pressed) { out.close(); return 0; }
					
					if (alignSeqs[j].active) {  //this sequence has not been merged yet
						//are you within "diff" bases
						int mismatch = calcMisMatches(alignSeqs[i].seq.getAligned(), alignSeqs[j].seq.getAligned());
						
						if (mismatch <= diffs) {
							//merge
							alignSeqs[i].names += ',' + alignSeqs[j].names;
							alignSeqs[i].numIdentical += alignSeqs[j].numIdentical;
							
							chunk += alignSeqs[j].seq.getName() + "\t" + toString(alignSeqs[j].numIdentical) + "\t" + toString(mismatch) + "\t" + alignSeqs[j].seq.getAligned() + "\n";
							
							alignSeqs[j].active = 0;
							alignSeqs[j].numIdentical = 0;
							count++;
						}
					}//end if j active
				}//end for loop j
				
				//remove from active list 
				alignSeqs[i].active = 0;
				
				out << "ideal_seq_" << (i+1) << '\t' << alignSeqs[i].numIdentical << endl << chunk << endl;;
				
			}//end if active i
			if(i % 100 == 0)	{ m->mothurOut(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); m->mothurOutEndLine();	}
		}
		out.close();
		
		if(numSeqs % 100 != 0)	{ m->mothurOut(toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); m->mothurOutEndLine();	}	
		
		return count;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::readFASTA(){
	try {
		//ifstream inNames;
		ifstream inFasta;
		
		//m->openInputFile(namefile, inNames);
		m->openInputFile(fastafile, inFasta);
		
		//string firstCol, secondCol, nameString;
		set<int> lengths;
		
		while (!inFasta.eof()) {
			
			if (m->control_pressed) { inFasta.close(); return 0; }
			
			//inNames >> firstCol >> secondCol;
			//nameString = secondCol;
			
			//m->gobble(inNames);
			//int size = 1;
			//while (secondCol.find_first_of(',') != -1) { 
			//	size++;
			//	secondCol = secondCol.substr(secondCol.find_first_of(',')+1, secondCol.length());
			//}
			
			Sequence seq(inFasta);  m->gobble(inFasta);
			
			if (seq.getName() != "") {  //can get "" if commented line is at end of fasta file
				if (namefile != "") {
					itSize = sizes.find(seq.getName());
					
					if (itSize == sizes.end()) { m->mothurOut(seq.getName() + " is not in your names file, please correct."); m->mothurOutEndLine(); exit(1); }
					else{
						seqPNode tempNode(itSize->second, seq, names[seq.getName()]);
						alignSeqs.push_back(tempNode);
						lengths.insert(seq.getAligned().length());
					}	
				}else { //no names file, you are identical to yourself 
					seqPNode tempNode(1, seq, seq.getName());
					alignSeqs.push_back(tempNode);
					lengths.insert(seq.getAligned().length());
				}
			}
		}
		inFasta.close();
		//inNames.close();
        
        if (lengths.size() > 1) { m->control_pressed = true; m->mothurOut("[ERROR]: your sequences are not all the same length. pre.cluster requires sequences to be aligned."); m->mothurOutEndLine(); }
        else if (lengths.size() == 1) { length = *(lengths.begin()); }
        
		return alignSeqs.size();
	}
	
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "readFASTA");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::loadSeqs(map<string, string>& thisName, vector<Sequence>& thisSeqs){
	try {
		set<int> lengths;
		alignSeqs.clear();
		map<string, string>::iterator it;
		bool error = false;
		 	
		for (int i = 0; i < thisSeqs.size(); i++) {
			
			if (m->control_pressed) { return 0; }
						
			if (namefile != "") {
				it = thisName.find(thisSeqs[i].getName());
				
				//should never be true since parser checks for this
				if (it == thisName.end()) { m->mothurOut(thisSeqs[i].getName() + " is not in your names file, please correct."); m->mothurOutEndLine(); error = true; }
				else{
					//get number of reps
					int numReps = 1;
					for(int j=0;j<(it->second).length();j++){
						if((it->second)[j] == ','){	numReps++;	}
					}
					
					seqPNode tempNode(numReps, thisSeqs[i], it->second);
					alignSeqs.push_back(tempNode);
                    lengths.insert(thisSeqs[i].getAligned().length());
				}	
			}else { //no names file, you are identical to yourself 
				seqPNode tempNode(1, thisSeqs[i], thisSeqs[i].getName());
				alignSeqs.push_back(tempNode);
				lengths.insert(thisSeqs[i].getAligned().length());
			}
		}
		
        if (lengths.size() > 1) { error = true; m->mothurOut("[ERROR]: your sequences are not all the same length. pre.cluster requires sequences to be aligned."); m->mothurOutEndLine(); }
        else if (lengths.size() == 1) { length = *(lengths.begin()); }
        
		//sanity check
		if (error) { m->control_pressed = true; }
		
		thisSeqs.clear();
		
		return alignSeqs.size();
	}
	
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "loadSeqs");
		exit(1);
	}
}
				
/**************************************************************************************************/

int PreClusterCommand::calcMisMatches(string seq1, string seq2){
	try {
		int numBad = 0;
		
		for (int i = 0; i < seq1.length(); i++) {
			//do they match
			if (seq1[i] != seq2[i]) { numBad++; }
			if (numBad > diffs) { return length;  } //to far to cluster
		}
		
		return numBad;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "calcMisMatches");
		exit(1);
	}
}

/**************************************************************************************************/

void PreClusterCommand::printData(string newfasta, string newname){
	try {
		ofstream outFasta;
		ofstream outNames;
		
		if (bygroup) {
			m->openOutputFileAppend(newfasta, outFasta);
			m->openOutputFileAppend(newname, outNames);
		}else {
			m->openOutputFile(newfasta, outFasta);
			m->openOutputFile(newname, outNames);
		}
		
		for (int i = 0; i < alignSeqs.size(); i++) {
			if (alignSeqs[i].numIdentical != 0) {
				alignSeqs[i].seq.printSequence(outFasta); 
				outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl;
			}
		}
		
		outFasta.close();
		outNames.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "printData");
		exit(1);
	}
}
/**************************************************************************************************/

void PreClusterCommand::readNameFile(){
	try {
		ifstream in;
		m->openInputFile(namefile, in);
		string firstCol, secondCol;
				
		while (!in.eof()) {
			in >> firstCol >> secondCol; m->gobble(in);
			names[firstCol] = secondCol;
			int size = 1;

			for(int i=0;i<secondCol.size();i++){
				if(secondCol[i] == ','){	size++;	}
			}
			sizes[firstCol] = size;
		}
		in.close();
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "readNameFile");
		exit(1);
	}
}

/**************************************************************************************************/


