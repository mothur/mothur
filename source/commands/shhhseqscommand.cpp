/*
 *  shhhseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "shhhseqscommand.h"



//**********************************************************************************************************************
vector<string> ShhhSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-map",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none","name",false,true,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter psigma("sigma", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(psigma);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ShhhSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The shhh.seqs command reads a fasta and name file and ....\n";
		helpString += "The shhh.seqs command parameters are fasta, name, group, sigma and processors.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file associated with your fasta file. It is required. \n";
		helpString += "The group parameter allows you to provide a group file.  When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
		helpString += "The sigma parameter ....  The default is 0.01. \n";
		helpString += "The shhh.seqs command should be in the following format: \n";
		helpString += "shhh.seqs(fasta=yourFastaFile, name=yourNameFile) \n";
		helpString += "Example: shhh.seqs(fasta=AD.align, name=AD.names) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ShhhSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],shhh_seqs.fasta";   }
        else if (type == "name")    {   pattern = "[filename],shhh_seqs.names";   }
        else if (type == "map")        {   pattern = "[filename],shhh_seqs.map";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ShhhSeqsCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************

ShhhSeqsCommand::ShhhSeqsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["map"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "ShhhSeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
ShhhSeqsCommand::ShhhSeqsCommand(string option) {
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
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }	
			else { current->setFastaFile(fastafile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found") { 			
				namefile = current->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current namefile and the name parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (namefile == "not open") { namefile =  ""; abort = true; }	
			else {  current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found") { groupfile =  "";   }
			else if (groupfile == "not open") { abort = true; groupfile =  ""; }	
			else {   current->setGroupFile(groupfile);  }
			
			string temp	= validParameter.valid(parameters, "sigma");		if(temp == "not found"){	temp = "0.01"; }
			util.mothurConvert(temp, sigma); 
			sigma = 1/sigma;
            
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
                    
			if (namefile == "") {
				vector<string> files; files.push_back(fastafile);
				if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "ShhhSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int ShhhSeqsCommand::execute() {
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (outputDir == "") { outputDir = util.hasPath(fastafile);  }//if user entered a file with a path then preserve it		
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
		string outputFileName = getOutputFileName("fasta",variables);
		string nameFileName = getOutputFileName("name",variables);
		string mapFileName = getOutputFileName("map",variables);
		
		if (groupfile != "") {
			//Parse sequences by group
            vector<string> temp;
			SequenceParser parser(groupfile, fastafile, namefile, temp);
			vector<string> groups = parser.getNamesOfGroups();
			
			if (m->getControl_pressed()) {  return 0; }
			
			//clears files
			ofstream out, out1, out2;
			util.openOutputFile(outputFileName, out); out.close(); 
			util.openOutputFile(nameFileName, out1); out1.close();
			mapFileName = outputDir + util.getRootName(util.getSimpleName(fastafile))  + "shhh.";
			
			vector<string> mapFileNames;
			if(processors == 1)	{	mapFileNames = driverGroups(parser, outputFileName, nameFileName, mapFileName, 0, groups.size(), groups);	}
			else				{	mapFileNames = createProcessesGroups(parser, outputFileName, nameFileName, mapFileName, groups);			}
			
			if (m->getControl_pressed()) {    return 0;	}	
			
			for (int j = 0; j < mapFileNames.size(); j++) { outputNames.push_back(mapFileNames[j]); outputTypes["map"].push_back(mapFileNames[j]); }
			
			//deconvolute results by running unique.seqs
			deconvoluteResults(outputFileName, nameFileName);
			
			if (m->getControl_pressed()) {   return 0;	}				
			
		}else{	
			vector<string> sequences;
			vector<string> uniqueNames;
			vector<string> redundantNames;
			vector<int> seqFreq;
			
			seqNoise noise;
			correctDist* correct = new correctDist(processors);
			
			//reads fasta and name file and loads them in order
			readData(correct, noise, sequences, uniqueNames, redundantNames, seqFreq);
			if (m->getControl_pressed()) { return 0; }
			
			//calc distances for cluster
			string distFileName = outputDir + util.getRootName(util.getSimpleName(fastafile)) + "shhh.dist";
			correct->execute(distFileName);
			delete correct;
			
			if (m->getControl_pressed()) { util.mothurRemove(distFileName); return 0; }
			
			driver(noise, sequences, uniqueNames, redundantNames, seqFreq, distFileName, outputFileName, nameFileName, mapFileName); 
			outputNames.push_back(mapFileName); outputTypes["map"].push_back(mapFileName);
		}
		
		if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
		
		outputNames.push_back(outputFileName); outputTypes["fasta"].push_back(outputFileName);
		outputNames.push_back(nameFileName); outputTypes["name"].push_back(nameFileName);
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int ShhhSeqsCommand::readData(correctDist* correct, seqNoise& noise, vector<string>& seqs, vector<string>& uNames, vector<string>& rNames, vector<int>& freq) {
	try {
		map<string, string> nameMap; 
		map<string, string>::iterator it;
		util.readNames(namefile, nameMap);
		bool error = false;
		
		ifstream in;
		util.openInputFile(fastafile, in);
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
			Sequence seq(in); util.gobble(in);
			
			if (seq.getName() != "") {
				correct->addSeq(seq.getName(), seq.getAligned());
				
				it = nameMap.find(seq.getName());
				if (it != nameMap.end()) {
					noise.addSeq(seq.getAligned(), seqs);
					noise.addRedundantName(it->first, it->second, uNames, rNames, freq);
				}else {
					m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your namefile, please correct.");
					error = true;
				}
			}
		}
		in.close();
		
		if (error) { m->setControl_pressed(true); }
		
		return seqs.size();
		
	}catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "readData");
		exit(1);
	}
}
//**********************************************************************************************************************
int ShhhSeqsCommand::loadData(correctDist* correct, seqNoise& noise, vector<string>& seqs, vector<string>& uNames, vector<string>& rNames, vector<int>& freq, map<string, string>& nameMap, vector<Sequence>& sequences) {
	try {
		bool error = false;
		map<string, string>::iterator it;
		
		for (int i = 0; i < sequences.size(); i++) {
			
			if (m->getControl_pressed()) { return 0; }
			
			if (sequences[i].getName() != "") {
				correct->addSeq(sequences[i].getName(), sequences[i].getAligned());
				
				it = nameMap.find(sequences[i].getName());
				if (it != nameMap.end()) {
					noise.addSeq(sequences[i].getAligned(), seqs);
					noise.addRedundantName(it->first, it->second, uNames, rNames, freq);
				}else {
					m->mothurOut("[ERROR]: " + sequences[i].getName() + " is in your fasta file and not in your namefile, please correct.");
					error = true;
				}
			}
		}
				
		if (error) { m->setControl_pressed(true); }
		
		return seqs.size();
		
	}catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "loadData");
		exit(1);
	}
}
/**************************************************************************************************/
vector<string> ShhhSeqsCommand::createProcessesGroups(SequenceParser& parser, string newFName, string newNName, string newMName, vector<string> groups) {
	try {
		
		vector<int> processIDS;
		int process = 1;
		vector<string> mapfileNames;
        bool recalc = false;
		
		//sanity check
		if (groups.size() < processors) { processors = groups.size(); }
		
		//divide the groups between the processors
		vector<linePair> lines;
		int remainingPairs = groups.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }

		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				mapfileNames = driverGroups(parser, newFName + toString(process) + ".temp", newNName + toString(process) + ".temp", newMName, lines[process].start, lines[process].end, groups);
				
				//pass filenames to parent
				ofstream out;
				string tempFile = newMName + toString(process) + ".temp";
				util.openOutputFile(tempFile, out);
				out << mapfileNames.size() << endl;
				for (int i = 0; i < mapfileNames.size(); i++) {
					out << mapfileNames[i] << endl;
				}
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->setControl_pressed(false);
                for (int i=0;i<processIDS.size();i++) {
                    util.mothurRemove(newFName + (toString(processIDS[i]) + ".temp"));
                    util.mothurRemove(newNName + (toString(processIDS[i]) + ".temp"));
                    util.mothurRemove(newMName + (toString(processIDS[i]) + ".temp"));
                }
                recalc = true;
                break;
			}
		}
		
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->setControl_pressed(false);  processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            lines.clear();
            int remainingPairs = groups.size();
            int startIndex = 0;
            for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
                int numPairs = remainingPairs; //case for last processor
                if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
                lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
                startIndex = startIndex + numPairs;
                remainingPairs = remainingPairs - numPairs;
            }
            
            mapfileNames.clear();
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    mapfileNames = driverGroups(parser, newFName + toString(process) + ".temp", newNName + toString(process) + ".temp", newMName, lines[process].start, lines[process].end, groups);
                    
                    //pass filenames to parent
                    ofstream out;
                    string tempFile = newMName + toString(process) + ".temp";
                    util.openOutputFile(tempFile, out);
                    out << mapfileNames.size() << endl;
                    for (int i = 0; i < mapfileNames.size(); i++) {
                        out << mapfileNames[i] << endl;
                    }
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

        
		//do my part
		mapfileNames = driverGroups(parser, newFName, newNName, newMName, lines[0].start, lines[0].end, groups);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			ifstream in;
			string tempFile =  newMName + toString(processIDS[i]) + ".temp";
			util.openInputFile(tempFile, in);
			if (!in.eof()) { 
				int tempNum = 0; in >> tempNum;  util.gobble(in);
				for (int j = 0; j < tempNum; j++) {
					string filename;
					in >> filename; util.gobble(in);
					mapfileNames.push_back(filename);
				}
			}
			in.close(); util.mothurRemove(tempFile);
			
		}
#else
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the shhhseqsData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<shhhseqsData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
			// Allocate memory for thread data.
			string extension = toString(i) + ".temp";

			shhhseqsData* tempShhhseqs = new shhhseqsData(fastafile, namefile, groupfile, (newFName+extension), (newNName+extension), newMName, groups, m, lines[i].start, lines[i].end, sigma, i);
			pDataArray.push_back(tempShhhseqs);
			processIDS.push_back(i);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyShhhSeqsThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
		
		//using the main process as a worker saves time and memory
		mapfileNames = driverGroups(parser, newFName, newNName, newMName, lines[0].start, lines[0].end, groups);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            if (pDataArray[i]->count != (pDataArray[i]->end-pDataArray[i]->start)) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end-pDataArray[i]->start) + " groups assigned to it, quitting. \n"); m->setControl_pressed(true); 
            }
			for (int j = 0; j < pDataArray[i]->mapfileNames.size(); j++) {
				mapfileNames.push_back(pDataArray[i]->mapfileNames[j]);
			}
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
#endif		
		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			util.appendFiles((newFName + toString(processIDS[i]) + ".temp"), newFName);
			util.mothurRemove((newFName + toString(processIDS[i]) + ".temp"));
			
			util.appendFiles((newNName + toString(processIDS[i]) + ".temp"), newNName);
			util.mothurRemove((newNName + toString(processIDS[i]) + ".temp"));
		}
		
		return mapfileNames;	
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "createProcessesGroups");
		exit(1);
	}
}
/**************************************************************************************************/
vector<string> ShhhSeqsCommand::driverGroups(SequenceParser& parser, string newFFile, string newNFile, string newMFile, int start, int end, vector<string> groups){
	try {
		
		vector<string> mapFileNames;
		
		for (int i = start; i < end; i++) {
			
			start = time(NULL);
			
			if (m->getControl_pressed()) {  return mapFileNames; }
            
            string lowerCaseName = groups[i];
            for (int j = 0; j < lowerCaseName.length(); j++) { lowerCaseName[j] = tolower(lowerCaseName[j]);    }
            
            if (lowerCaseName == "ignore") {   }
            else {
                
                m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[i] + ":"); m->mothurOutEndLine();
                
                map<string, string> thisNameMap;
                thisNameMap = parser.getNameMap(groups[i]);
                vector<Sequence> thisSeqs = parser.getSeqs(groups[i]);
                
                vector<string> sequences;
                vector<string> uniqueNames;
                vector<string> redundantNames;
                vector<int> seqFreq;
                
                seqNoise noise;
                correctDist* correct = new correctDist(1); //we use one processor since we already split up the work load.
                
                //load this groups info in order
                loadData(correct, noise, sequences, uniqueNames, redundantNames, seqFreq, thisNameMap, thisSeqs);
                if (m->getControl_pressed()) { return mapFileNames; }
                
                //calc distances for cluster
                string distFileName = outputDir + util.getRootName(util.getSimpleName(fastafile)) + groups[i] + ".shhh.dist";
                correct->execute(distFileName);
                delete correct;
                
                if (m->getControl_pressed()) { util.mothurRemove(distFileName); return mapFileNames; }
                
                driver(noise, sequences, uniqueNames, redundantNames, seqFreq, distFileName, newFFile+groups[i], newNFile+groups[i], newMFile+groups[i]+".map");
                
                if (m->getControl_pressed()) { return mapFileNames; }
                
                util.appendFiles(newFFile+groups[i], newFFile); util.mothurRemove(newFFile+groups[i]);
                util.appendFiles(newNFile+groups[i], newNFile); util.mothurRemove(newNFile+groups[i]);
                mapFileNames.push_back(newMFile+groups[i]+".map");
                
                m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process group " + groups[i] + "."); m->mothurOutEndLine();
            }
		}
		
		return mapFileNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "driverGroups");
		exit(1);
	}
}
//**********************************************************************************************************************
int ShhhSeqsCommand::driver(seqNoise& noise, 
							vector<string>& sequences, 
							vector<string>& uniqueNames, 
							vector<string>& redundantNames, 
							vector<int>& seqFreq, 
							string distFileName, string outputFileName, string nameFileName, string mapFileName) {
	try {
		double cutOff = 0.08;
		int minIter = 10;
		int maxIter = 1000;
		double minDelta = 1e-6;
		int numIters = 0;
		double maxDelta = 1e6;
		int numSeqs = sequences.size();
				
		//run cluster command
		string inputString = "phylip=" + distFileName + ", method=furthest, cutoff=0.08";
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		m->mothurOut("Running command: cluster(" + inputString + ")"); m->mothurOutEndLine(); 
		
		Command* clusterCommand = new ClusterCommand(inputString);
		clusterCommand->execute();
		
		map<string, vector<string> > filenames = clusterCommand->getOutputFiles();
		string listFileName = filenames["list"][0];
		string rabundFileName = filenames["rabund"][0]; util.mothurRemove(rabundFileName);
		string sabundFileName = filenames["sabund"][0]; util.mothurRemove(sabundFileName);
	
		delete clusterCommand;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		
		if (m->getControl_pressed()) { util.mothurRemove(distFileName); util.mothurRemove(listFileName); return 0; }
		
		vector<double> distances(numSeqs * numSeqs);
		noise.getDistanceData(distFileName, distances);
		util.mothurRemove(distFileName); 
		if (m->getControl_pressed()) { util.mothurRemove(listFileName); return 0; }
		
		vector<int> otuData(numSeqs);
		vector<int> otuFreq;
		vector<vector<int> > otuBySeqLookUp;
		noise.getListData(listFileName, cutOff, otuData, otuFreq, otuBySeqLookUp);
		util.mothurRemove(listFileName);
		if (m->getControl_pressed()) { return 0; }
		
		int numOTUs = otuFreq.size();
		
		vector<double> weights(numOTUs, 0);
		vector<int> change(numOTUs, 1);
		vector<int> centroids(numOTUs, -1);
		vector<int> cumCount(numOTUs, 0);
		
		vector<double> tau(numSeqs, 1);
		vector<int> anP(numSeqs, 0);
		vector<int> anI(numSeqs, 0);
		vector<int> anN(numSeqs, 0);
		vector<vector<int> > aanI = otuBySeqLookUp;
		
		while(numIters < minIter || ((maxDelta > minDelta) && (numIters < maxIter))){
			
			if (m->getControl_pressed()) { return 0; }
			
			noise.updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount); if (m->getControl_pressed()) { return 0; }
			maxDelta = noise.calcNewWeights(weights, seqFreq, anI, cumCount, anP, otuFreq, tau);  if (m->getControl_pressed()) { return 0; }
			
			noise.calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau); if (m->getControl_pressed()) { return 0; }
			noise.checkCentroids(weights, centroids); if (m->getControl_pressed()) { return 0; }
			 
			otuFreq.assign(numOTUs, 0);
			
			int total = 0;
			
			for(int i=0;i<numSeqs;i++){
				if (m->getControl_pressed()) { return 0; }
				
				double offset = 1e6;
				double norm = 0.0000;
				double minWeight = 0.1;
				vector<double> currentTau(numOTUs);
				
				for(int j=0;j<numOTUs;j++){
					if (m->getControl_pressed()) { return 0; }
					if(weights[j] > minWeight && distances[i * numSeqs+centroids[j]] < offset){
						offset = distances[i * numSeqs+centroids[j]];
					}
				}
				
				for(int j=0;j<numOTUs;j++){
					if (m->getControl_pressed()) { return 0; }
					if(weights[j] > minWeight){
						currentTau[j] = exp(sigma * (-distances[(i * numSeqs + centroids[j])] + offset)) * weights[j];
						norm += currentTau[j];
					}
					else{
						currentTau[j] = 0.0000;
					}
				}			
				
				for(int j=0;j<numOTUs;j++){
					if (m->getControl_pressed()) { return 0; }
					currentTau[j] /= norm;
				}
				
				for(int j=0;j<numOTUs;j++){
					if (m->getControl_pressed()) { return 0; }
					
					if(currentTau[j] > 1.0e-4){
						int oldTotal = total;
						total++;
						
						tau.resize(oldTotal+1);
						tau[oldTotal] = currentTau[j];
						otuBySeqLookUp[j][otuFreq[j]] = oldTotal;
						aanI[j][otuFreq[j]] = i;
						otuFreq[j]++;
						
					}
				}
				
				anP.resize(total);
				anI.resize(total);
			}
			
			numIters++;
		}
		
		noise.updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount);  if (m->getControl_pressed()) { return 0; }
		
		vector<double> percentage(numSeqs);
		noise.setUpOTUData(otuData, percentage, cumCount, tau, otuFreq, anP, anI);  if (m->getControl_pressed()) { return 0; }
		noise.finishOTUData(otuData, otuFreq, anP, anI, cumCount, otuBySeqLookUp, aanI, tau);  if (m->getControl_pressed()) { return 0; }
		
		change.assign(numOTUs, 1);
		noise.calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau); if (m->getControl_pressed()) { return 0; }
		
		
		vector<int> finalTau(numOTUs, 0);
		for(int i=0;i<numSeqs;i++){
			if (m->getControl_pressed()) { return 0; }
			finalTau[otuData[i]] += int(seqFreq[i]);
		}
		
		noise.writeOutput(outputFileName, nameFileName, mapFileName, finalTau, centroids, otuData, sequences, uniqueNames, redundantNames, seqFreq, distances);
		
		return 0;
		
	}catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "driver");
		exit(1);
	}
}	
//**********************************************************************************************************************
int ShhhSeqsCommand::deconvoluteResults(string fastaFile, string nameFile){
	try {
		m->mothurOutEndLine(); m->mothurOut("Deconvoluting results:"); m->mothurOutEndLine(); m->mothurOutEndLine();
		
		//use unique.seqs to create new name and fastafile
		string inputString = "fasta=" + fastaFile + ", name=" + nameFile;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
		current->setMothurCalling(true);
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		current->setMothurCalling(false);
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		
		string newnameFile = filenames["name"][0];
		string newfastaFile = filenames["fasta"][0];
		
		util.mothurRemove(fastaFile); rename(newfastaFile.c_str(), fastaFile.c_str()); 
		if (nameFile != newnameFile) { util.mothurRemove(nameFile); rename(newnameFile.c_str(), nameFile.c_str()); }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhhSeqsCommand", "deconvoluteResults");
		exit(1);
	}
}
//**********************************************************************************************************************



