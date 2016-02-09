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
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pdiffs("diffs", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter palign("align", "Multiple", "needleman-gotoh-blast-noalign", "needleman", "", "", "","",false,false); parameters.push_back(palign);
        CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
        CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
        CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
        CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);

        CommandParameter ptopdown("topdown", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(ptopdown);

		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The pre.cluster command parameters are fasta, name, group, count, topdown, processors and diffs. The fasta parameter is required. \n";
		helpString += "The name parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n";
		helpString += "The group parameter allows you to provide a group file so you can cluster by group. \n";
        helpString += "The count parameter allows you to provide a count file so you can cluster by group. \n";
		helpString += "The diffs parameter allows you to specify maximum number of mismatched bases allowed between sequences in a grouping. The default is 1.\n";
        helpString += "The topdown parameter allows you to specify whether to cluster from largest abundance to smallest or smallest to largest.  Default=T, meaning largest to smallest.\n";
        helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n";
        helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
        helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
        helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
        helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
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
string PreClusterCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],precluster,[extension]"; } 
        else if (type == "name") {  pattern = "[filename],precluster.names"; } 
        else if (type == "count") {  pattern = "[filename],precluster.count_table"; }
        else if (type == "map") {  pattern =  "[filename],precluster.map"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "getOutputPattern");
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
        outputTypes["count"] = tempOutNames;
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
            outputTypes["count"] = tempOutNames;
		
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
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
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not found") { countfile =  "";   }
			else if (countfile == "not open") { abort = true; countfile =  ""; }	
			else {   
                m->setCountTableFile(countfile); 
                ct.readTable(countfile, true, false);
                if (ct.hasGroupInfo()) { bygroup = true; }
                else { bygroup = false;  }
            }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			
			string temp	= validParameter.validFile(parameters, "diffs", false);		if(temp == "not found"){	temp = "1"; }
			m->mothurConvert(temp, diffs); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
            temp = validParameter.validFile(parameters, "topdown", false);		if(temp == "not found"){  temp = "T"; }
			topdown = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
            m->mothurConvert(temp, match);
            
            temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
            m->mothurConvert(temp, misMatch);
            if (misMatch > 0) { m->mothurOut("[ERROR]: mismatch must be negative.\n"); abort=true; }
            
            temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
            m->mothurConvert(temp, gapOpen);
            if (gapOpen > 0) { m->mothurOut("[ERROR]: gapopen must be negative.\n"); abort=true; }
            
            temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
            m->mothurConvert(temp, gapExtend);
            if (gapExtend > 0) { m->mothurOut("[ERROR]: gapextend must be negative.\n"); abort=true; }

            align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
            
            method = "unaligned";
            
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    parser.getNameFile(files);
                }
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
        
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, 1000);	}
        else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 1000);			}
        else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
        else if(align == "noalign")		{	alignment = new NoAlign();													}
        else {
            m->mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
            m->mothurOutEndLine();
            alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 1000);
        }
		
		string fileroot = outputDir + m->getRootName(m->getSimpleName(fastafile));
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
		string newNamesFile = getOutputFileName("name",variables);
        string newCountFile = getOutputFileName("count",variables);
		string newMapFile = getOutputFileName("map",variables); //add group name if by group
        variables["[extension]"] = m->getExtension(fastafile);
		string newFastaFile = getOutputFileName("fasta", variables);
		outputNames.push_back(newFastaFile); outputTypes["fasta"].push_back(newFastaFile);
		if (countfile == "") { outputNames.push_back(newNamesFile); outputTypes["name"].push_back(newNamesFile); }
		else { outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile); }
		
		if (bygroup) {
			//clear out old files
			ofstream outFasta; m->openOutputFile(newFastaFile, outFasta); outFasta.close();
			ofstream outNames; m->openOutputFile(newNamesFile, outNames);  outNames.close();
			newMapFile = fileroot + "precluster.";
			
			//parse fasta and name file by group
            vector<string> groups;
			if (countfile != "") {
                cparser = new SequenceCountParser(countfile, fastafile);
                groups = cparser->getNamesOfGroups();
            }else {
                if (namefile != "") { parser = new SequenceParser(groupfile, fastafile, namefile);	}
                else				{ parser = new SequenceParser(groupfile, fastafile);			}
                groups = parser->getNamesOfGroups();
			}
            
			if(processors == 1)	{	driverGroups(newFastaFile, newNamesFile, newMapFile, 0, groups.size(), groups);	}
			else				{	createProcessesGroups(newFastaFile, newNamesFile, newMapFile, groups);			}
			
			if (countfile != "") { 
                mergeGroupCounts(newCountFile, newNamesFile, newFastaFile);
                delete cparser; 
            }else {  
                delete parser; 
                //run unique.seqs for deconvolute results
                string inputString = "fasta=" + newFastaFile;
                if (namefile != "") { inputString += ", name=" + newNamesFile; }
                m->mothurOutEndLine(); 
                m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
                m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
                m->mothurCalling = true;
                
                Command* uniqueCommand = new DeconvoluteCommand(inputString);
                uniqueCommand->execute();
                
                map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
                
                delete uniqueCommand;
                m->mothurCalling = false;
                m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
                
                m->renameFile(filenames["fasta"][0], newFastaFile);
                m->renameFile(filenames["name"][0], newNamesFile); 
			}
            if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}	 delete alignment; return 0; }
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run pre.cluster."); m->mothurOutEndLine(); 
				
		}else {
            if (processors != 1) { m->mothurOut("When using running without group information mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
			if (namefile != "") { readNameFile(); }
		
			//reads fasta file and return number of seqs
			int numSeqs = readFASTA(); //fills alignSeqs and makes all seqs active
		
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} delete alignment; return 0; }
	
			if (numSeqs == 0) { m->mothurOut("Error reading fasta file...please correct."); m->mothurOutEndLine(); delete alignment; return 0;  }
			if (diffs > length) { m->mothurOut("Error: diffs is greater than your sequence length."); m->mothurOutEndLine(); delete alignment; return 0;  }
			
			int count = process(newMapFile);
			outputNames.push_back(newMapFile); outputTypes["map"].push_back(newMapFile);
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} delete alignment; return 0; }
			
			m->mothurOut("Total number of sequences before precluster was " + toString(alignSeqs.size()) + "."); m->mothurOutEndLine();
			m->mothurOut("pre.cluster removed " + toString(count) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine(); 
			if (countfile != "") { newNamesFile = newCountFile; }
            printData(newFastaFile, newNamesFile, "");
            			
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); 
		}
				
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} delete alignment; return 0; }
        
        delete alignment;
		
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
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::createProcessesGroups(string newFName, string newNName, string newMFile, vector<string> groups) {
	try {
		
		vector<int> processIDS;
		int process = 1;
		int num = 0;
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
                outputNames.clear();
				num = driverGroups(newFName + m->mothurGetpid(process) + ".temp", newNName + m->mothurGetpid(process) + ".temp", newMFile, lines[process].start, lines[process].end, groups);
                
                string tempFile = m->mothurGetpid(process) + ".outputNames.temp";
                ofstream outTemp;
                m->openOutputFile(tempFile, outTemp);
                
                outTemp << outputNames.size();
                for (int i = 0; i < outputNames.size(); i++) { outTemp << outputNames[i] << endl; }
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
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".outputNames.temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".outputNames.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            lines.clear();
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            int remainingPairs = groups.size();
            int startIndex = 0;
            for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
                int numPairs = remainingPairs; //case for last processor
                if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
                lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
                startIndex = startIndex + numPairs;
                remainingPairs = remainingPairs - numPairs;
            }
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    outputNames.clear();
                    num = driverGroups(newFName + m->mothurGetpid(process) + ".temp", newNName + m->mothurGetpid(process) + ".temp", newMFile, lines[process].start, lines[process].end, groups);
                    
                    string tempFile = m->mothurGetpid(process) + ".outputNames.temp";
                    ofstream outTemp;
                    m->openOutputFile(tempFile, outTemp);
                    
                    outTemp << outputNames.size();
                    for (int i = 0; i < outputNames.size(); i++) { outTemp << outputNames[i] << endl; }
                    outTemp.close();
                    
                    exit(0);
                }else {
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

		
		//do my part
		num = driverGroups(newFName, newNName, newMFile, lines[0].start, lines[0].end, groups);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            string tempFile = toString(processIDS[i]) +  ".outputNames.temp";
            ifstream intemp;
            m->openInputFile(tempFile, intemp);
            
            int num;
            intemp >> num;
            for (int k = 0; k < num; k++) {
                string name = "";
                intemp >> name; m->gobble(intemp);
                
                outputNames.push_back(name); outputTypes["map"].push_back(name);
            }
            intemp.close();
            m->mothurRemove(tempFile);
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
			
			preClusterData* tempPreCluster = new preClusterData(fastafile, namefile, groupfile, countfile, (newFName+extension), (newNName+extension), newMFile, groups, m, lines[i].start, lines[i].end, diffs, topdown, i, method, align, match, misMatch, gapOpen, gapExtend);
			pDataArray.push_back(tempPreCluster);
			processIDS.push_back(i);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyPreclusterThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
				
		//using the main process as a worker saves time and memory
		num = driverGroups(newFName, newNName, newMFile, lines[0].start, lines[0].end, groups);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            if (pDataArray[i]->count != (pDataArray[i]->end-pDataArray[i]->start)) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end-pDataArray[i]->start) + " groups assigned to it, quitting. \n"); m->control_pressed = true; 
            }
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
int PreClusterCommand::driverGroups(string newFFile, string newNFile, string newMFile, int start, int end, vector<string> groups){
	try {
		
		int numSeqs = 0;
		
		//precluster each group
		for (int i = start; i < end; i++) {
			
			start = time(NULL);
			
			if (m->control_pressed) {  return 0; }
			
			m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[i] + ":"); m->mothurOutEndLine();
			
			map<string, string> thisNameMap;
            vector<Sequence> thisSeqs;
			if (groupfile != "") { 
                thisSeqs = parser->getSeqs(groups[i]);
            }else if (countfile != "") {
                thisSeqs = cparser->getSeqs(groups[i]);
            }
			if (namefile != "") {  thisNameMap = parser->getNameMap(groups[i]); }
            
			//fill alignSeqs with this groups info.
			numSeqs = loadSeqs(thisNameMap, thisSeqs, groups[i]);
			
			if (m->control_pressed) {   return 0; }
			
            if (method == "aligned") { if (diffs > length) { m->mothurOut("Error: diffs is greater than your sequence length."); m->mothurOutEndLine(); m->control_pressed = true; return 0;  } }
			
			int count= process(newMFile+groups[i]+".map");
			outputNames.push_back(newMFile+groups[i]+".map"); outputTypes["map"].push_back(newMFile+groups[i]+".map");
			
			if (m->control_pressed) {  return 0; }
			
			m->mothurOut("Total number of sequences before pre.cluster was " + toString(alignSeqs.size()) + "."); m->mothurOutEndLine();
			m->mothurOut("pre.cluster removed " + toString(count) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine(); 
			printData(newFFile, newNFile, groups[i]);
			
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
        if (topdown) { sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityTopDown);  }
        else {  sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityDownTop);  }
		
		int count = 0;
		int numSeqs = alignSeqs.size();
		
        if (topdown) {
            //think about running through twice...
            for (int i = 0; i < numSeqs; i++) {
                
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
                                alignSeqs[j].diffs = mismatch;
                                count++;
                            }
                        }//end if j active
                    }//end for loop j
                    
                    //remove from active list 
                    alignSeqs[i].active = 0;
                    
                    out << "ideal_seq_" << (i+1) << '\t' << alignSeqs[i].numIdentical << endl << chunk << endl;
                    
                }//end if active i
                if(i % 100 == 0)	{ m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
            }
        }else {
            map<int, string> mapFile;
            map<int, int> originalCount;
            map<int, int>::iterator itCount;
            for (int i = 0; i < numSeqs; i++) { mapFile[i] = ""; originalCount[i] = alignSeqs[i].numIdentical; }
            
            //think about running through twice...
            for (int i = 0; i < numSeqs; i++) {
                
                //try to merge it into larger seqs
                for (int j = i+1; j < numSeqs; j++) {
                    
                    if (m->control_pressed) { out.close(); return 0; }
                    
                    if (originalCount[j] > originalCount[i]) {  //this sequence is more abundant than I am
                        //are you within "diff" bases
                        int mismatch = calcMisMatches(alignSeqs[i].seq.getAligned(), alignSeqs[j].seq.getAligned());
                        
                        if (mismatch <= diffs) {
                            //merge
                            alignSeqs[j].names += ',' + alignSeqs[i].names;
                            alignSeqs[j].numIdentical += alignSeqs[i].numIdentical;
                            
                            mapFile[j] = alignSeqs[i].seq.getName() + "\t" + toString(alignSeqs[i].numIdentical) + "\t" + toString(mismatch) + "\t" + alignSeqs[i].seq.getAligned() + "\n" + mapFile[i];
                            alignSeqs[i].numIdentical = 0;
                            originalCount.erase(i);
                            mapFile[i] = "";
                            count++;
                            j+=numSeqs; //exit search, we merged this one in.
                        }
                    }//end abundance check
                }//end for loop j
                
                if(i % 100 == 0)	{ m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
            }
            
            for (int i = 0; i < numSeqs; i++) {
                if (alignSeqs[i].numIdentical != 0) {
                    out << "ideal_seq_" << (i+1) << '\t' << alignSeqs[i].numIdentical << endl  << alignSeqs[i].seq.getName() + "\t" + toString(alignSeqs[i].numIdentical) + "\t" + toString(0) + "\t" + alignSeqs[i].seq.getAligned() + "\n" << mapFile[i] << endl;
                }
            }
            
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
		
		m->openInputFile(fastafile, inFasta);
		set<int> lengths;
		
		while (!inFasta.eof()) {
			
			if (m->control_pressed) { inFasta.close(); return 0; }
						
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
                    int numRep = 1;
                    if (countfile != "") { numRep = ct.getNumSeqs(seq.getName()); }
					seqPNode tempNode(numRep, seq, seq.getName());
					alignSeqs.push_back(tempNode);
					lengths.insert(seq.getAligned().length());
				}
			}
		}
		inFasta.close();
        
        if (lengths.size() > 1) { method = "unaligned"; }
        else if (lengths.size() == 1) {  method = "aligned"; }
        
        length = *(lengths.begin());
        
		return alignSeqs.size();
	}
	
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "readFASTA");
		exit(1);
	}
}
/**************************************************************************************************/
int PreClusterCommand::loadSeqs(map<string, string>& thisName, vector<Sequence>& thisSeqs, string group){
	try {
		set<int> lengths;
		alignSeqs.clear();
		map<string, string>::iterator it;
		bool error = false;
        map<string, int> thisCount;
        if (countfile != "") { thisCount = cparser->getCountTable(group);  }
        	
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
                int numRep = 1;
                if (countfile != "") { 
                    map<string, int>::iterator it2 = thisCount.find(thisSeqs[i].getName());
                    
                    //should never be true since parser checks for this
                    if (it2 == thisCount.end()) { m->mothurOut(thisSeqs[i].getName() + " is not in your count file, please correct."); m->mothurOutEndLine(); error = true; }
                    else { numRep = it2->second;  }
                }
				seqPNode tempNode(numRep, thisSeqs[i], thisSeqs[i].getName());
				alignSeqs.push_back(tempNode);
				lengths.insert(thisSeqs[i].getAligned().length());
			}
		}
    
        if (lengths.size() > 1) { method = "unaligned"; }
        else if (lengths.size() == 1) {  method = "aligned"; }
        
        length = *(lengths.begin());
        
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
		
        if (method == "unaligned") {
            //align to eachother
            Sequence seqI("seq1", seq1);
            Sequence seqJ("seq2", seq2);
            
            //align seq2 to seq1 - less abundant to more abundant
            alignment->align(seqJ.getUnaligned(), seqI.getUnaligned());
            seq2 = alignment->getSeqAAln();
            seq1 = alignment->getSeqBAln();
            
            //chop gap ends
            int startPos = 0;
            int endPos = seq2.length()-1;
            for (int i = 0; i < seq2.length(); i++) {  if (isalpha(seq2[i])) { startPos = i; break; } }
            for (int i = seq2.length()-1; i >= 0; i--) {  if (isalpha(seq2[i])) { endPos = i; break; } }
            
            //count number of diffs
            for (int i = startPos; i <= endPos; i++) {
                if (seq2[i] != seq1[i]) { numBad++; }
                if (numBad > diffs) { return length;  } //to far to cluster
            }

        }else {
            //count diffs
            for (int i = 0; i < seq1.length(); i++) {
                //do they match
                if (seq1[i] != seq2[i]) { numBad++; }
                if (numBad > diffs) { return length;  } //to far to cluster
            }
        }
		return numBad;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "calcMisMatches");
		exit(1);
	}
}
/**************************************************************************************************/

int PreClusterCommand::mergeGroupCounts(string newcount, string newname, string newfasta){
	try {
		ifstream inNames;
        m->openInputFile(newname, inNames);
        
        string group, first, second;
        set<string> uniqueNames;
        while (!inNames.eof()) {
            if (m->control_pressed) { break; }
            inNames >> group; m->gobble(inNames);
            inNames >> first; m->gobble(inNames);
            inNames >> second; m->gobble(inNames);
            
            vector<string> names;
            m->splitAtComma(second, names);
            
            uniqueNames.insert(first);
            
            int total = ct.getGroupCount(first, group);
            for (int i = 1; i < names.size(); i++) {
                total += ct.getGroupCount(names[i], group);
                ct.setAbund(names[i], group, 0);
            }
            ct.setAbund(first, group, total);
        }
        inNames.close();
        
        vector<string> namesOfSeqs = ct.getNamesOfSeqs();
        for (int i = 0; i < namesOfSeqs.size(); i++) {
            if (ct.getNumSeqs(namesOfSeqs[i]) == 0) {
                ct.remove(namesOfSeqs[i]);
            }
        }
        
        ct.printTable(newcount); 
        m->mothurRemove(newname);
        
        if (bygroup) { //if by group, must remove the duplicate seqs that are named the same
            ifstream in;
            m->openInputFile(newfasta, in);
            
            ofstream out;
            m->openOutputFile(newfasta+"temp", out);
            
            int count = 0;
            set<string> already;
            while(!in.eof()) {
                if (m->control_pressed) { break; }
                
                Sequence seq(in); m->gobble(in);
                
                if (seq.getName() != "") {
                    count++;
                    if (already.count(seq.getName()) == 0) {
                        seq.printSequence(out);
                        already.insert(seq.getName());
                    }
                }
            }
            in.close();
            out.close();
            m->mothurRemove(newfasta);
            m->renameFile(newfasta+"temp", newfasta);
        }
		        return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "mergeGroupCounts");
		exit(1);
	}
}

/**************************************************************************************************/

void PreClusterCommand::printData(string newfasta, string newname, string group){
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
		
        if ((countfile != "") && (group == ""))  { outNames << "Representative_Sequence\ttotal\n";  }
		for (int i = 0; i < alignSeqs.size(); i++) {
			if (alignSeqs[i].numIdentical != 0) {
				alignSeqs[i].seq.printSequence(outFasta); 
				if (countfile != "") {  
                    if (group != "") {  outNames << group << '\t' << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl; }
                    else {  outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].numIdentical << endl;  }
                }else {  outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl;  }
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
            
            m->checkName(firstCol);
            m->checkName(secondCol);
            int size = m->getNumNames(secondCol);
            
			names[firstCol] = secondCol;
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


