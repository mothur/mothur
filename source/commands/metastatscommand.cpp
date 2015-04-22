/*
 *  metastatscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/16/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "metastatscommand.h"
#include "sharedutilities.h"
#include "sharedrabundfloatvector.h"


//**********************************************************************************************************************
vector<string> MetaStatsCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","metastats",false,true,true); parameters.push_back(pshared);
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pdesign);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pthreshold("threshold", "Number", "", "0.05", "", "", "","",false,false); parameters.push_back(pthreshold);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter psets("sets", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MetaStatsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "This command is based on the Metastats program, White, J.R., Nagarajan, N. & Pop, M. Statistical methods for detecting differentially abundant features in clinical metagenomic samples. PLoS Comput Biol 5, e1000352 (2009).\n";
		helpString += "The metastats command outputs a .metastats file. \n";
		helpString += "The metastats command parameters are shared, iters, threshold, groups, label, design, sets and processors.  The shared and design parameters are required, unless you have valid current files.\n";
		helpString += "The design parameter allows you to assign your groups to sets when you are running metastat. mothur will run all pairwise comparisons of the sets. It is required. \n";
		helpString += "The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the group name and the second column is the set the group belongs to.\n";
		helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile.\n";
		helpString += "The iters parameter allows you to set number of bootstrap permutations for estimating null distribution of t statistic.  The default is 1000. \n";
		helpString += "The threshold parameter allows you to set the significance level to reject null hypotheses (default 0.05).\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
		helpString += "The metastats command should be in the following format: metastats(design=yourDesignFile).\n";
		helpString += "Example metastats(design=temp.design, groups=A-B-C).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MetaStatsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "metastats") {  pattern = "[filename],[distance],[group],metastats"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MetaStatsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MetaStatsCommand::MetaStatsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["metastats"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "MetaStatsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MetaStatsCommand::MetaStatsCommand(string option) {
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
			
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["metastats"] = tempOutNames;
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
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
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") {  				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
			
			//check for required parameters
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { abort = true; }
			else if (designfile == "not found") {  				
				//if there is a current design file, use it
				designfile = m->getDesignFile(); 
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current designfile and the design parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setDesignFile(designfile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(sharedfile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
			
			sets = validParameter.validFile(parameters, "sets", false);			
			if (sets == "not found") { sets = ""; }
			else { 
				m->splitAtDash(sets, Sets);
			}

			
			string temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "threshold", false);			if (temp == "not found") { temp = "0.05"; }
			m->mothurConvert(temp, threshold); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
		}

	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "MetaStatsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MetaStatsCommand::execute(){
	try {
	
        if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        //just used to convert files to test metastats online
        /****************************************************/
        bool convertInputToShared = false;
        convertSharedToInput = false;
        if (convertInputToShared) { convertToShared(sharedfile); return 0; }
        /****************************************************/
		
		designMap = new GroupMap(designfile);
		designMap->readDesignMap();
	
		input = new InputData(sharedfile, "sharedfile");
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//setup the pairwise comparions of sets for metastats
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		//make sure sets are all in designMap
		SharedUtil* util = new SharedUtil(); 
		vector<string> dGroups = designMap->getNamesOfGroups();
		util->setGroups(Sets, dGroups);  
		delete util;
		
		int numGroups = Sets.size();
		for (int a=0; a<numGroups; a++) { 
			for (int l = 0; l < a; l++) {	
				vector<string> groups; groups.push_back(Sets[a]); groups.push_back(Sets[l]);
				namesOfGroupCombos.push_back(groups);
			}
		}
	
		
		//only 1 combo
		if (numGroups == 2) { processors = 1; }
		else if (numGroups < 2)	{ m->mothurOut("Not enough sets, I need at least 2 valid sets. Unable to complete command."); m->mothurOutEndLine(); m->control_pressed = true; }

        if(processors != 1){
            int remainingPairs = namesOfGroupCombos.size();
            int startIndex = 0;
            for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
                int numPairs = remainingPairs; //case for last processor
                if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
                lines.push_back(linePair(startIndex, numPairs)); //startIndex, numPairs
                startIndex = startIndex + numPairs;
                remainingPairs = remainingPairs - numPairs;
            }
        }
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  outputTypes.clear(); for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } m->clearGroups(); delete input; delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
	
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
			
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				lookup = input->getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
			
			if (m->control_pressed) {  outputTypes.clear(); m->clearGroups(); delete input;  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

			//get next line to process
			lookup = input->getSharedRAbundVectors();				
		}
		
		if (m->control_pressed) { outputTypes.clear(); m->clearGroups(); delete input; delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }

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
			for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
			lookup = input->getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			
			process(lookup);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
	
		//reset groups parameter
		m->clearGroups();  
		delete input; 
		delete designMap;
		
        if (m->control_pressed) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int MetaStatsCommand::process(vector<SharedRAbundVector*>& thisLookUp){
	try {
		
		
				if(processors == 1){
					driver(0, namesOfGroupCombos.size(), thisLookUp);
				}else{
					int process = 1;
					vector<int> processIDS;
                    bool recalc = false;
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
					//loop through and create all the processes you want
					while (process != processors) {
						pid_t pid = fork();
			
						if (pid > 0) {
							processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
							process++;
						}else if (pid == 0){
							driver(lines[process].start, lines[process].end, thisLookUp);
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
                            recalc = true;
                            break;
						}
					}
					
                    if (recalc) {
                        //test line, also set recalc to true.
                        //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
                        
                        //redo file divide
                        lines.clear();
                        int remainingPairs = namesOfGroupCombos.size();
                        int startIndex = 0;
                        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
                            int numPairs = remainingPairs; //case for last processor
                            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
                            lines.push_back(linePair(startIndex, numPairs)); //startIndex, numPairs
                            startIndex = startIndex + numPairs;
                            remainingPairs = remainingPairs - numPairs;
                        }
                        
                        processIDS.resize(0);
                        process = 1;
                        
                        //loop through and create all the processes you want
                        while (process != processors) {
                            pid_t pid = fork();
                            
                            if (pid > 0) {
                                processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                                process++;
                            }else if (pid == 0){
                                driver(lines[process].start, lines[process].end, thisLookUp);
                                exit(0);
                            }else {
                                m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                                exit(0);
                            }
                        }
                    }

                    
					//do my part
					driver(lines[0].start, lines[0].end, thisLookUp);
		
					//force parent to wait until all the processes are done
					for (int i=0;i<(processors-1);i++) { 
						int temp = processIDS[i];
						wait(&temp);
					}
        #else
                    
                    //////////////////////////////////////////////////////////////////////////////////////////////////////
                    //Windows version shared memory, so be careful when passing variables through the summarySharedData struct. 
                    //Above fork() will clone, so memory is separate, but that's not the case with windows, 
                    //Taking advantage of shared memory to pass results vectors.
                    //////////////////////////////////////////////////////////////////////////////////////////////////////
                    
                    vector<metastatsData*> pDataArray; 
                    DWORD   dwThreadIdArray[processors-1];
                    HANDLE  hThreadArray[processors-1]; 
                    
                    //Create processor worker threads.
                    for( int i=1; i<processors; i++ ){
                        
                        //make copy of lookup so we don't get access violations
                        vector<SharedRAbundVector*> newLookup;
                        vector<string> designMapGroups;
                        for (int k = 0; k < thisLookUp.size(); k++) {
                            SharedRAbundVector* temp = new SharedRAbundVector();
                            temp->setLabel(thisLookUp[k]->getLabel());
                            temp->setGroup(thisLookUp[k]->getGroup());
                            newLookup.push_back(temp);
                            designMapGroups.push_back(designMap->getGroup(thisLookUp[k]->getGroup()));
                        }
                        
                        //for each bin
                        for (int k = 0; k < thisLookUp[0]->getNumBins(); k++) {
                            if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
                            for (int j = 0; j < thisLookUp.size(); j++) { newLookup[j]->push_back(thisLookUp[j]->getAbundance(k), thisLookUp[j]->getGroup()); }
                        }
                        
                        // Allocate memory for thread data.
                        metastatsData* tempSum = new metastatsData(sharedfile, outputDir, m, lines[i].start, lines[i].end, namesOfGroupCombos, newLookup, designMapGroups, iters, threshold);
                        pDataArray.push_back(tempSum);
                        processIDS.push_back(i);
                        
                        hThreadArray[i-1] = CreateThread(NULL, 0, MyMetastatsThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
                    }
                    
                    //do my part
					driver(lines[0].start, lines[0].end, thisLookUp);
                    
                    //Wait until all threads have terminated.
                    WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
                    
                    //Close all thread handles and free memory allocations.
                    for(int i=0; i < pDataArray.size(); i++){
                        if (pDataArray[i]->count != (pDataArray[i]->num)) {
                            m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->num) + " groups assigned to it, quitting. \n"); m->control_pressed = true; 
                        }
                        for (int j = 0; j < pDataArray[i]->thisLookUp.size(); j++) {  delete pDataArray[i]->thisLookUp[j];  } 
                        for (int j = 0; j < pDataArray[i]->outputNames.size(); j++) {  
                            outputNames.push_back(pDataArray[i]->outputNames[j]);
                            outputTypes["metastats"].push_back(pDataArray[i]->outputNames[j]);
                        }
                                                
                        CloseHandle(hThreadArray[i]);
                        delete pDataArray[i];
                    }
        #endif

				}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************
int MetaStatsCommand::driver(unsigned long long start, unsigned long long num, vector<SharedRAbundVector*>& thisLookUp) {
	try {
	
		//for each combo
		for (int c = start; c < (start+num); c++) {
			
			//get set names
			string setA = namesOfGroupCombos[c][0];
			string setB = namesOfGroupCombos[c][1];
		
			//get filename
            map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
            variables["[distance]"] = thisLookUp[0]->getLabel();
            variables["[group]"] = setA + "-" + setB;
			string outputFileName = getOutputFileName("metastats",variables);
			outputNames.push_back(outputFileName); outputTypes["metastats"].push_back(outputFileName);
			//int nameLength = outputFileName.length();
			//char * output = new char[nameLength];
			//strcpy(output, outputFileName.c_str());
	
			//build matrix from shared rabunds
			//double** data;
			//data = new double*[thisLookUp[0]->getNumBins()];
			
			vector< vector<double> > data2; data2.resize(thisLookUp[0]->getNumBins());
			
			vector<SharedRAbundVector*> subset;
			int setACount = 0;
			int setBCount = 0;
			for (int i = 0; i < thisLookUp.size(); i++) {
				string thisGroup = thisLookUp[i]->getGroup();
				
				//is this group for a set we want to compare??
				//sorting the sets by putting setB at the back and setA in the front
				if ((designMap->getGroup(thisGroup) == setB)) {  
					subset.push_back(thisLookUp[i]);
					setBCount++;
				}else if ((designMap->getGroup(thisGroup) == setA)) {
					subset.insert(subset.begin()+setACount, thisLookUp[i]);
					setACount++;
				}
			}
						
			if ((setACount == 0) || (setBCount == 0))  { 
				m->mothurOut("Missing shared info for " + setA + " or " + setB + ". Skipping comparison."); m->mothurOutEndLine(); 
				outputNames.pop_back();
			}else {
                
				//fill data
				for (int j = 0; j < thisLookUp[0]->getNumBins(); j++) {
					//data[j] = new double[subset.size()];
					data2[j].resize(subset.size(), 0.0);
                   
					for (int i = 0; i < subset.size(); i++) {
						data2[j][i] = (subset[i]->getAbundance(j));
					}
				}
				
				m->mothurOut("Comparing " + setA + " and " + setB + "..."); m->mothurOutEndLine(); 
				//metastat_main(output, thisLookUp[0]->getNumBins(), subset.size(), threshold, iters, data, setACount);
                if (convertSharedToInput) { convertToInput(subset, outputFileName);  }
				
				m->mothurOutEndLine();
				MothurMetastats mothurMeta(threshold, iters);
				mothurMeta.runMetastats(outputFileName , data2, setACount);
				m->mothurOutEndLine();
				m->mothurOutEndLine(); 
			}
			
			//free memory
			//delete output;
			//for(int i = 0; i < thisLookUp[0]->getNumBins(); i++)  {  delete[] data[i];  }
			//delete[] data; 
		}
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
/*Metastats files look like:
 13_0	14_0	13_52	14_52	70S	71S	72S	M1	M2	M3	C11	C12	C21	C15	C16	C19	C3	C4	C9
 Alphaproteobacteria	0	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0
 Mollicutes	0	0	2	0	0	59	5	11	4	1	0	2	8	1	0	1	0	3	0
 Verrucomicrobiae	0	0	0	0	0	1	6	0	0	0	0	0	0	0	0	0	0	0	0
 Deltaproteobacteria	0	0	0	0	0	6	1	0	1	0	1	1	7	0	0	0	0	0	0
 Cyanobacteria	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
 Epsilonproteobacteria	0	0	0	0	0	0	0	0	6	0	0	3	1	0	0	0	0	0	0
 Clostridia	75	65	207	226	801	280	267	210	162	197	81	120	106	148	120	94	84	98	121
 Bacilli	3	2	16	8	21	52	31	70	46	65	4	28	5	23	62	26	20	30	25
 Bacteroidetes (class)	21	25	22	64	226	193	296	172	98	55	19	149	201	85	50	76	113	92	82
 Gammaproteobacteria	0	0	0	0	0	1	0	0	0	0	1	1	0	0	0	1	0	0	0
 TM7_genera_incertae_sedis	0	0	0	0	0	0	0	0	1	0	1	2	0	2	0	0	0	0	0
 Actinobacteria (class)	1	1	1	2	0	0	0	9	3	7	1	1	1	3	1	2	1	2	3
 Betaproteobacteria	0	0	3	3	0	0	9	1	1	0	1	2	3	1	1	0	0	0	0
*/
//this function is just used to convert files to test the differences between the metastats version and mothurs version
int MetaStatsCommand::convertToShared(string filename) {
	try {
        ifstream in;
        m->openInputFile(filename, in);
        
        string header = m->getline(in); m->gobble(in);
        
        vector<string> groups = m->splitWhiteSpace(header);
        vector<SharedRAbundFloatVector*> newLookup;
        cout << groups.size() << endl;
        for (int i = 0; i < groups.size(); i++) {
            cout << "creating group " << groups[i] << endl;
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
            temp->setLabel("0.03");
            temp->setGroup(groups[i]);
            newLookup.push_back(temp);
        }
        
        int otuCount = 0;
        while (!in.eof()) {
            if (m->control_pressed) { break; }
            
            string otuname;
            in >> otuname; m->gobble(in);
            otuCount++;
            cout << otuname << endl;
            for (int i = 0; i < groups.size(); i++) {
                double temp;
                in >> temp; m->gobble(in);
                newLookup[i]->push_back(temp, groups[i]);
            }
            m->gobble(in);
        }
        in.close();
    
        ofstream out;
        m->openOutputFile(filename+".shared", out);
        
        out << "label\tgroup\tnumOTUs\t";
        
        string snumBins = toString(otuCount);
        for (int i = 0; i < otuCount; i++) {
            string binLabel = "Otu";
            string sbinNumber = toString(i+1);
            if (sbinNumber.length() < snumBins.length()) {
                int diff = snumBins.length() - sbinNumber.length();
                for (int h = 0; h < diff; h++) { binLabel += "0"; }
            }
            binLabel += sbinNumber;
            out << binLabel << '\t';
        }
        out << endl;
        
        for (int i = 0; i < groups.size(); i++) {
            out << "0.03" << '\t' << groups[i] << '\t';
            newLookup[i]->print(out);
        }
        out.close();
        
        cout << filename+".shared" << endl;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "convertToShared");
		exit(1);
	}
}
//**********************************************************************************************************************

int MetaStatsCommand::convertToInput(vector<SharedRAbundVector*>& subset, string thisfilename) {
	try {
        ofstream out;
        m->openOutputFile(thisfilename+".matrix", out);
        
        out << "\t";
        for (int i = 0; i < subset.size()-1; i++) {
            out << subset[i]->getGroup() << '\t';
        }
        out << subset[subset.size()-1]->getGroup() << endl;
        
        for (int i = 0; i < subset[0]->getNumBins(); i++) {
            out << m->currentSharedBinLabels[i] << '\t';
            for (int j = 0; j < subset.size()-1; j++) {
                out << subset[j]->getAbundance(i) << '\t';
            }
            out << subset[subset.size()-1]->getAbundance(i) << endl;
        }
        out.close();
        
        cout << thisfilename+".matrix" << endl;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MetaStatsCommand", "convertToInput");
		exit(1);
	}
}

//**********************************************************************************************************************
