/*
 *  countseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "countseqscommand.h"
#include "sharedutilities.h"
#include "counttable.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> CountSeqsCommand::setParameters(){	
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "NameSHared-sharedGroup", "NameSHared", "none","count",false,false,true); parameters.push_back(pshared);
		CommandParameter pname("name", "InputTypes", "", "", "NameSHared", "NameSHared", "none","count",false,false,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "sharedGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The count.seqs aka. make.table command reads a name or shared file and outputs a .count_table file.  You may also provide a group with the names file to get the counts broken down by group.\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include in the counts, by default all groups in your groupfile are used.\n";
		helpString += "When you use the groups parameter and a sequence does not represent any sequences from the groups you specify it is not included in the .count.summary file.\n";
        helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The count.seqs command should be in the following format: count.seqs(name=yourNameFile).\n";
		helpString += "Example count.seqs(name=amazon.names) or make.table(name=amazon.names).\n";
		helpString += "Note: No spaces between parameter labels (i.e. name), '=' and parameters (i.e.yourNameFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        if (type == "count") {  pattern = "[filename],count_table-[filename],[distance],count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "CountSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
CountSeqsCommand::CountSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "CountSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

CountSeqsCommand::CountSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;
        allLines = 1;
		
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
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["count"] = tempOutNames;
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
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
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){	namefile = ""; }
            else { m->setNameFile(namefile); }
            
            sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found"){	sharedfile = ""; }
            else { m->setSharedFile(sharedfile); }
            
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { m->setGroupFile(groupfile); }
            
            if ((namefile == "") && (sharedfile == "")) {
                namefile = m->getNameFile();
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else {
                    sharedfile = m->getSharedFile();
                    if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
                    else {
                        m->mothurOut("You have no current namefile or sharedfile and the name or shared parameter is required."); m->mothurOutEndLine(); abort = true;
                    }
                }
			}

			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "all"; }
			m->splitAtDash(groups, Groups);
            m->setGroups(Groups);
            
            string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}

		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "CountSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
        processors=1;
#endif
		
        map<string, string> variables;

        if (namefile != "") {
            unsigned long long total = 0;
            int start = time(NULL);
            if (outputDir == "") { outputDir = m->hasPath(namefile); }
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(namefile));
            string outputFileName = getOutputFileName("count", variables);
            
            total = process(outputFileName);
            
            if (m->control_pressed) { m->mothurRemove(outputFileName); return 0; }
            
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create a table for " + toString(total) + " sequences.");
            m->mothurOutEndLine(); m->mothurOutEndLine();
            
            m->mothurOutEndLine();
            m->mothurOut("Total number of sequences: " + toString(total)); m->mothurOutEndLine();
 
        }else {
            if (outputDir == "") { outputDir = m->hasPath(sharedfile); }
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
            
            InputData input(sharedfile, "sharedfile");
            SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
            string lastLabel = lookup->getLabel();
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->control_pressed) { delete lookup; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
                
                if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
                    
                    m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                    processShared(data, variables);
                    for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear();
                    
                    processedLabels.insert(lookup->getLabel());
                    userLabels.erase(lookup->getLabel());
                }
                
                if ((m->anyLabelsToProcess(lookup->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = lookup->getLabel();
                    
                    delete lookup;
                    lookup = input.getSharedRAbundVectors(lastLabel);
                    m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    
                    vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                    processShared(data, variables);
                    for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear();
                    
                    processedLabels.insert(lookup->getLabel());
                    userLabels.erase(lookup->getLabel());
                    
                    //restore real lastlabel to save below
                    lookup->setLabels(saveLabel);
                }
                
                lastLabel = lookup->getLabel();
                //prevent memory leak
                delete lookup;
                
                if (m->control_pressed) { return 0; }
                
                //get next line to process
                lookup = input.getSharedRAbundVectors();
            }
            
            if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }
            
            //output error messages about any remaining user labels
            set<string>::iterator it;
            bool needToRun = false;
            for (it = userLabels.begin(); it != userLabels.end(); it++) {
                m->mothurOut("Your file does not include the label " + *it);
                if (processedLabels.count(lastLabel) != 1) {
                    m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
                    needToRun = true;
                }else {
                    m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
                }
            }
            
            //run last label if you need to
            if (needToRun == true)  {
                delete lookup;
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                
                vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                processShared(data, variables);
                for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear();
                
               delete lookup;
            }
            
        }
        
        //set rabund file as new current rabundfile
		itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { string current = (itTypes->second)[0]; m->setCountTableFile(current); }
		}
        
        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
		m->mothurOutEndLine();
        
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

unsigned long long CountSeqsCommand::processShared(vector<SharedRAbundVector*>& lookup, map<string, string> variables){
    try {
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("count", variables);
        outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
        
        ofstream out;
        m->openOutputFile(outputFileName, out);
        
        out << "OTU_Label\ttotal";
        for (int i = 0; i < lookup.size(); i++) { out << '\t' << lookup[i]->getGroup(); } out << endl;
        
        for (int j = 0; j < lookup[0]->getNumBins(); j++) {
            if (m->control_pressed) { break; }
            
            int total = 0;
            string output = "";
            for (int i = 0; i < lookup.size(); i++) {
                total += lookup[i]->get(j);
                output += '\t' + toString(lookup[i]->get(j));
            }
            out << m->currentSharedBinLabels[j] << '\t' << total << output << endl;
        }
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountSeqsCommand", "processShared");
        exit(1);
    }
}
//**********************************************************************************************************************

unsigned long long CountSeqsCommand::process(string outputFileName){
	try {
        ofstream out;
        m->openOutputFile(outputFileName, out); outputTypes["count"].push_back(outputFileName);
        outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
		out << "Representative_Sequence\ttotal";
        
        GroupMap* groupMap;
		if (groupfile != "") { 
			groupMap = new GroupMap(groupfile); groupMap->readMap(); 
			
			//make sure groups are valid. takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			vector<string> nameGroups = groupMap->getNamesOfGroups();
			util->setGroups(Groups, nameGroups);
			delete util;
			
			//sort groupNames so that the group title match the counts below, this is needed because the map object automatically sorts
			sort(Groups.begin(), Groups.end());
			
			//print groupNames
			for (int i = 0; i < Groups.size(); i++) {
				out << '\t' << Groups[i];
			}
		}
		out << endl;
        out.close();
        
        unsigned long long total = createProcesses(groupMap, outputFileName);
        
        if (groupfile != "") { delete groupMap; }
        
        return total;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "processSmall");
		exit(1);
	}
}
/**************************************************************************************************/
unsigned long long CountSeqsCommand::createProcesses(GroupMap*& groupMap, string outputFileName) {
	try {
		
		vector<int> processIDS;
		int process = 0;
        vector<unsigned long long> positions;
        vector<linePair> lines;
        unsigned long long numSeqs = 0;
        bool recalc = false;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		positions = m->divideFilePerLine(namefile, processors);
		for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
#else
		if(processors == 1){ lines.push_back(linePair(0, 1000));  }
        else {
            unsigned long long numSeqs = 0;
            positions = m->setFilePosEachLine(namefile, numSeqs);
            if (positions.size() < processors) { processors = positions.size(); }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif

        		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		
		//loop through and create all the processes you want
		while (process != processors-1) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
                string filename = m->mothurGetpid(process) + ".temp";
				numSeqs = driver(lines[process].start, lines[process].end, filename, groupMap);
                
                string tempFile = m->mothurGetpid(process) + ".num.temp";
                ofstream outTemp;
                m->openOutputFile(tempFile, outTemp);
                
                outTemp << numSeqs << endl;
                outTemp.close();
                
				exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".temp"));
                    m->mothurRemove((toString(processIDS[i]) + ".num.temp"));
                }
                m->control_pressed = false;
                recalc = true;
                break;
            }
		}
		
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".temp"));m->mothurRemove((toString(processIDS[i]) + ".num.temp"));}m->control_pressed = false;  processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            positions.clear();
            lines.clear();
            positions = m->divideFilePerLine(namefile, processors);
            for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
            
            numSeqs = 0;
            processIDS.resize(0);
            process = 0;
            
            while (process != processors-1) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    string filename = m->mothurGetpid(process) + ".temp";
                    numSeqs = driver(lines[process].start, lines[process].end, filename, groupMap);
                    
                    string tempFile = m->mothurGetpid(process) + ".num.temp";
                    ofstream outTemp;
                    m->openOutputFile(tempFile, outTemp);
                    
                    outTemp << numSeqs << endl;
                    outTemp.close();
                    
                    exit(0);
                }else {
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

        
		string filename = m->mothurGetpid(process) + ".temp";
        numSeqs = driver(lines[processors-1].start, lines[processors-1].end, filename, groupMap);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) {
			int temp = processIDS[i];
			wait(&temp);
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            string tempFile = toString(processIDS[i]) +  ".num.temp";
            ifstream intemp;
            m->openInputFile(tempFile, intemp);
            
            int num;
            intemp >> num; intemp.close();
            numSeqs += num;
            m->mothurRemove(tempFile);
        }
#else		
		vector<countData*> pDataArray;
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1];
        vector<GroupMap*> copies;
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
			string filename = toString(i) + ".temp";
            
            GroupMap* copyGroup = new GroupMap();
            copyGroup->getCopy(groupMap);
            copies.push_back(copyGroup);
            vector<string> cGroups = Groups;
           
			countData* temp = new countData(filename, copyGroup, m, lines[i].start, lines[i].end, groupfile, namefile, cGroups);
			pDataArray.push_back(temp);
			processIDS.push_back(i);
			
			hThreadArray[i] = CreateThread(NULL, 0, MyCountThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);
		}
		
		string filename = toString(processors-1) + ".temp";
        numSeqs = driver(lines[processors-1].start, lines[processors-1].end, filename, groupMap);
        		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            numSeqs += pDataArray[i]->total;
            delete copies[i];
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif
		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((toString(processIDS[i]) + ".temp"), outputFileName);
			m->mothurRemove((toString(processIDS[i]) + ".temp"));
		}
        m->appendFiles(filename, outputFileName);
        m->mothurRemove(filename);

        
        //sanity check
        if (groupfile != "") {
            if (numSeqs != groupMap->getNumSeqs()) {
                m->mothurOut("[ERROR]: processes reported processing " + toString(numSeqs) + " sequences, but group file indicates you have " + toString(groupMap->getNumSeqs()) + " sequences.");
                if (processors == 1) { m->mothurOut(" Could you have a file mismatch?\n"); }
                else { m->mothurOut(" Either you have a file mismatch or a process failed to complete the task assigned to it.\n"); m->control_pressed = true; }
            }
		}
		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
unsigned long long CountSeqsCommand::driver(unsigned long long start, unsigned long long end, string outputFileName, GroupMap*& groupMap) {
	try {
        
        ofstream out;
        m->openOutputFile(outputFileName, out);
        
        ifstream in;
		m->openInputFile(namefile, in);
		in.seekg(start);
        
        //adjust start if null strings
        if (start == 0) {  m->zapGremlins(in); m->gobble(in);  }

		bool done = false;
        unsigned long long total = 0;
		while (!done) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol; m->gobble(in); in >> secondCol; m->gobble(in);
            //cout << firstCol << '\t' << secondCol << endl;
            m->checkName(firstCol);
            m->checkName(secondCol);
			//cout << firstCol << '\t' << secondCol << endl;
            
			vector<string> names;
			m->splitAtChar(secondCol, names, ',');
			
			if (groupfile != "") {
				//set to 0
				map<string, int> groupCounts;
				int total = 0;
				for (int i = 0; i < Groups.size(); i++) { groupCounts[Groups[i]] = 0; }
				
				//get counts for each of the users groups
				for (int i = 0; i < names.size(); i++) {
					string group = groupMap->getGroup(names[i]);
					
					if (group == "not found") { m->mothurOut("[ERROR]: " + names[i] + " is not in your groupfile, please correct."); m->mothurOutEndLine(); }
					else {
						map<string, int>::iterator it = groupCounts.find(group);
						
						//if not found, then this sequence is not from a group we care about
						if (it != groupCounts.end()) {
							it->second++;
							total++;
						}
					}
				}
				
				if (total != 0) {
					out << firstCol << '\t' << total;
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << '\t' << it->second;
					}
					out << endl;
				}
			}else {
				out << firstCol << '\t' << names.size() << endl;
			}
			
			total += names.size();
            
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= end)) { break; }
#else
            if (in.eof()) { break; }
#endif

		}
		in.close();
        out.close();
        
        return total;

    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
map<int, string> CountSeqsCommand::processNameFile(string name) {
	try {
        map<int, string> indexToNames;
        
        ofstream out;
        m->openOutputFile(name, out);
        
        //open input file
		ifstream in;
		m->openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        int count = 0;
        
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    m->checkName(firstCol);
                    m->checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    m->splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  out << theseNames[i] << '\t' << count << endl;  }
                    indexToNames[count] = firstCol;
                    pairDone = false; 
                    count++;
                }
            }
		}
		in.close();
       
		
        if (rest != "") {
            vector<string> pieces = m->splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    m->checkName(firstCol);
                    m->checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    m->splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  out << theseNames[i] << '\t' << count << endl;  }
                    indexToNames[count] = firstCol;
                    pairDone = false; 
                    count++;
                }
            }

        }
        out.close();
        
        return indexToNames;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "processNameFile");
		exit(1);
	}
}
/**************************************************************************************************/
map<int, string> CountSeqsCommand::getGroupNames(string filename, set<string>& namesOfGroups) {
	try {
        map<int, string> indexToGroups;
        map<string, int> groupIndex;
        map<string, int>::iterator it;
        
        ofstream out;
        m->openOutputFile(filename, out);
        
        //open input file
		ifstream in;
		m->openInputFile(groupfile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        int count = 0;
        
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    m->checkName(firstCol);
                    it = groupIndex.find(secondCol);
                    if (it == groupIndex.end()) { //add group, assigning the group and number so we can use vectors above
                        groupIndex[secondCol] = count;
                        count++;
                    }
                    out << firstCol << '\t' << groupIndex[secondCol] << endl; 
                    namesOfGroups.insert(secondCol);
                    pairDone = false; 
                }
            }
		}
		in.close();
        
        
        if (rest != "") {
            vector<string> pieces = m->splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    m->checkName(firstCol);
                    it = groupIndex.find(secondCol);
                    if (it == groupIndex.end()) { //add group, assigning the group and number so we can use vectors above
                        groupIndex[secondCol] = count;
                        count++;
                    }
                    out << firstCol << '\t' << groupIndex[secondCol] << endl; 
                    namesOfGroups.insert(secondCol);
                    pairDone = false; 
                }
            }
        }
        out.close();
		
        for (it = groupIndex.begin(); it != groupIndex.end(); it++) {  indexToGroups[it->second] = it->first;  }
        
        return indexToGroups;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getGroupNames");
		exit(1);
	}
}
//**********************************************************************************************************************



