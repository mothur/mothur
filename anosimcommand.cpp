/*
 *  anosimcommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "anosimcommand.h"
#include "sharedutilities.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "sharedjackknife.h"
#include "whittaker.h"
#include "odum.h"
#include "canberra.h"
#include "structeuclidean.h"
#include "structchord.h"
#include "hellinger.h"
#include "manhattan.h"
#include "structpearson.h"
#include "soergel.h"
#include "spearman.h"
#include "structkulczynski.h"
#include "structchi2.h"
#include "speciesprofile.h"
#include "hamming.h"
#include "gower.h"
#include "memchi2.h"
#include "memchord.h"
#include "memeuclidean.h"
#include "mempearson.h"

//**********************************************************************************************************************
vector<string> AnosimCommand::getValidParameters(){	
	try {
		string Array[] =  {"groups","label","outputdir","iters","phylip","design","sets","processors","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
AnosimCommand::AnosimCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["anosim"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "AnosimCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AnosimCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"design"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AnosimCommand::getRequiredFiles(){	
	try {
		string Array[] =  {};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

AnosimCommand::AnosimCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false; calledHelp = false;   
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label","outputdir","iters","phylip","design","sets","processors","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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
			outputTypes["anosim"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
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
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
			}
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  globaldata->newRead(); format = "phylip";  globaldata->setPhylipFile(phylipfile);	}
			
			//check for required parameters
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { abort = true; }
			else if (designfile == "not found") {  designfile = "";  m->mothurOut("You must provide an design file."); m->mothurOutEndLine(); abort = true; }	
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getSharedFile() == "")) {
				if ((phylipfile == "")) {
					m->mothurOut("You must read a list and a group, a shared file or provide a distance matrix before you can use the anosim command."); m->mothurOutEndLine(); abort = true; 
				}
			}else { sharedfile = globaldata->getSharedFile(); } 
			
			//use distance matrix if one is provided
			if ((sharedfile != "") && (phylipfile != "")) { sharedfile = ""; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if (label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "";  }
			else { 
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			sets = validParameter.validFile(parameters, "sets", false);			
			if (sets == "not found") { sets = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				m->splitAtDash(sets, Sets);
			}
			
			
			string temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "processors", false);			if (temp == "not found"){	temp = "1";	}
			convert(temp, processors); 
			
			vector<string> Estimators;
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "morisitahorn";  }
			m->splitAtDash(calc, Estimators);
			
			if (abort == false) {
				
				ValidCalculators* validCalculator = new ValidCalculators();
				
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("treegroup", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedsobs") { 
							calculators.push_back(new SharedSobsCS());
						}else if (Estimators[i] == "sharedchao") { 
							calculators.push_back(new SharedChao1());
						}else if (Estimators[i] == "sharedace") { 
							calculators.push_back(new SharedAce());
						}else if (Estimators[i] == "jabund") { 	
							calculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							calculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							calculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							calculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							calculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							calculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							calculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							calculators.push_back(new ThetaN());
						}else if (Estimators[i] == "kstest") { 
							calculators.push_back(new KSTest());
						}else if (Estimators[i] == "sharednseqs") { 
							calculators.push_back(new SharedNSeqs());
						}else if (Estimators[i] == "ochiai") { 
							calculators.push_back(new Ochiai());
						}else if (Estimators[i] == "anderberg") { 
							calculators.push_back(new Anderberg());
						}else if (Estimators[i] == "kulczynski") { 
							calculators.push_back(new Kulczynski());
						}else if (Estimators[i] == "kulczynskicody") { 
							calculators.push_back(new KulczynskiCody());
						}else if (Estimators[i] == "lennon") { 
							calculators.push_back(new Lennon());
						}else if (Estimators[i] == "morisitahorn") { 
							calculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							calculators.push_back(new BrayCurtis());
						}else if (Estimators[i] == "whittaker") { 
							calculators.push_back(new Whittaker());
						}else if (Estimators[i] == "odum") { 
							calculators.push_back(new Odum());
						}else if (Estimators[i] == "canberra") { 
							calculators.push_back(new Canberra());
						}else if (Estimators[i] == "structeuclidean") { 
							calculators.push_back(new StructEuclidean());
						}else if (Estimators[i] == "structchord") { 
							calculators.push_back(new StructChord());
						}else if (Estimators[i] == "hellinger") { 
							calculators.push_back(new Hellinger());
						}else if (Estimators[i] == "manhattan") { 
							calculators.push_back(new Manhattan());
						}else if (Estimators[i] == "structpearson") { 
							calculators.push_back(new StructPearson());
						}else if (Estimators[i] == "soergel") { 
							calculators.push_back(new Soergel());
						}else if (Estimators[i] == "spearman") { 
							calculators.push_back(new Spearman());
						}else if (Estimators[i] == "structkulczynski") { 
							calculators.push_back(new StructKulczynski());
						}else if (Estimators[i] == "speciesprofile") { 
							calculators.push_back(new SpeciesProfile());
						}else if (Estimators[i] == "hamming") { 
							calculators.push_back(new Hamming());
						}else if (Estimators[i] == "structchi2") { 
							calculators.push_back(new StructChi2());
						}else if (Estimators[i] == "gower") { 
							calculators.push_back(new Gower());
						}else if (Estimators[i] == "memchi2") { 
							calculators.push_back(new MemChi2());
						}else if (Estimators[i] == "memchord") { 
							calculators.push_back(new MemChord());
						}else if (Estimators[i] == "memeuclidean") { 
							calculators.push_back(new MemEuclidean());
						}else if (Estimators[i] == "mempearson") { 
							calculators.push_back(new MemPearson());
						}
					}
				}
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "AnosimCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void AnosimCommand::help(){
	try {
		m->mothurOut("The anosim command can only be executed after a successful read.otu command of a list and group or shared file, or by providing a phylip formatted distance matrix.\n");
		m->mothurOut("The anosim command outputs a .anosim file. \n");
		m->mothurOut("The anosim command parameters are phylip, iters, groups, label, design, sets and processors.  The design parameter is required.\n");
		m->mothurOut("The design parameter allows you to assign your groups to sets when you are running anosim. It is required. \n");
		m->mothurOut("The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the group name and the second column is the set the group belongs to.\n");
		m->mothurOut("The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile. To run the pairwise comparisons of all sets use sets=all.\n");
		m->mothurOut("The iters parameter allows you to set number of randomization for the P value.  The default is 1000. \n");
		m->mothurOut("The groups parameter allows you to specify which of the groups you would like included. The group names are separated by dashes. groups=all will find all pairwise comparisons. \n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		m->mothurOut("The anosim command should be in the following format: anosim(design=yourDesignFile).\n");
		m->mothurOut("Example anosim(design=temp.design, groups=A-B-C).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

AnosimCommand::~AnosimCommand(){}

//**********************************************************************************************************************

int AnosimCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read design file
		designMap = new GroupMap(designfile);
		designMap->readDesignMap();
		
		//make sure sets are all in designMap
		SharedUtil* util = new SharedUtil();  
		util->setGroups(Sets, designMap->namesOfGroups);  
		delete util;
		
		//if the user pickedGroups or set sets=all, then we want to figure out how to split the pairwise comparison between processors
		int numGroups = Sets.size();
		if (pickedGroups) {
			for (int a=0; a<numGroups; a++) { 
				for (int l = 0; l < a; l++) {	
					vector<string> groups; groups.push_back(Sets[a]); groups.push_back(Sets[l]);
					namesOfGroupCombos.push_back(groups);
				}
			}
		}else { //one giant compare
			vector<string> groups;
			for (int a=0; a<numGroups; a++) { groups.push_back(Sets[a]); }
			namesOfGroupCombos.push_back(groups);
		}
		
		//only 1 combo
		if (numGroups == 2) { processors = 1; }
		else if (numGroups < 2)	{ m->mothurOut("Not enough sets, I need at least 2 valid sets. Unable to complete command."); m->mothurOutEndLine(); m->control_pressed = true; }
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors != 1){
			int numPairs = namesOfGroupCombos.size();
			int numPairsPerProcessor = numPairs / processors;
			
			for (int i = 0; i < processors; i++) {
				int startPos = i * numPairsPerProcessor;
				if(i == processors - 1){
					numPairsPerProcessor = numPairs - i * numPairsPerProcessor;
				}
				lines.push_back(linePair(startPos, numPairsPerProcessor));
			}
		}
#endif
		
		if (sharedfile != "") { //create distance matrix for each label
			
			if (outputDir == "") { outputDir = m->hasPath(sharedfile); }
			
			//for each calculator
			for(int i = 0 ; i < calculators.size(); i++) {
				
				//create a new filename
				ofstream out;
				string outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + calculators[i]->getName() + ".anosim";				
				m->openOutputFile(outputFile, out);
				outputNames.push_back(outputFile); outputTypes["anosim"].push_back(outputFile);
				
				//print headers
				out << "label\tgroupsCompared\tRValue\tpValue" << endl;  
				m->mothurOut("\nlabel\tgroupsCompared\tRValue\tpValue"); m->mothurOutEndLine();  
				out.close();
			}
			
			InputData input(sharedfile, "sharedfile");
			vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
			string lastLabel = lookup[0]->getLabel();
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			//sanity check between shared file groups and design file groups
			for (int i = 0; i < lookup.size(); i++) { 
				string group = designMap->getGroup(lookup[i]->getGroup());
				
				if (group == "not found") { m->control_pressed = true;  m->mothurOut("[ERROR]: " + lookup[i]->getGroup() + " is not in your design file, please correct."); m->mothurOutEndLine();  }
			}
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } globaldata->Groups.clear();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
				
				if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
					
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					process(lookup);
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup[0]->getLabel();
					
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
					lookup = input.getSharedRAbundVectors(lastLabel);
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
				
				if (m->control_pressed) {  globaldata->Groups.clear();   delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
				
				//get next line to process
				lookup = input.getSharedRAbundVectors();				
			}
			
			if (m->control_pressed) { globaldata->Groups.clear();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }  return 0; }
			
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
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				process(lookup);
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			}
			
			//reset groups parameter
			globaldata->Groups.clear();  
			
			
		}else { //user provided distance matrix
			
			if (outputDir == "") { outputDir = m->hasPath(phylipfile); }
			
			//create a new filename
			ofstream out;
			string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipfile))  + "anosim";				
			m->openOutputFile(outputFile, out);
			outputNames.push_back(outputFile); outputTypes["anosim"].push_back(outputFile);
			
			//print headers
			out << "groupsCompared\tRValue\tpValue" << endl; 
			m->mothurOut("\ngroupsCompared\tRValue\tpValue"); m->mothurOutEndLine();  
			out.close();
			
			ReadPhylipVector readMatrix(phylipfile);
			vector< vector<double> > completeMatrix;
			vector<string> names = readMatrix.read(completeMatrix);
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				driver(0, namesOfGroupCombos.size(), names, "", completeMatrix);
			}else{
				int process = 1;
				vector<int> processIDS;
				
				//loop through and create all the processes you want
				while (process != processors) {
					int pid = fork();
					
					if (pid > 0) {
						processIDS.push_back(pid);  
						process++;
					}else if (pid == 0){
						driver(lines[process].start, lines[process].num, names, toString(getpid()), completeMatrix);
						exit(0);
					}else { 
						m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
						for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
						exit(0);
					}
				}
				
				//do my part
				driver(lines[0].start, lines[0].num, names, "", completeMatrix);
				
				//force parent to wait until all the processes are done
				for (int i=0;i<(processors-1);i++) { 
					int temp = processIDS[i];
					wait(&temp);
				}
				
				//append files
				string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipfile)) + "anosim";				
				for (int j = 0; j < processIDS.size(); j++) {
					m->appendFiles((outputFile + toString(processIDS[j])), outputFile);
					remove((outputFile + toString(processIDS[j])).c_str());
				}
				
			}
#else
			driver(0, namesOfGroupCombos.size(), names, "", completeMatrix);
#endif
			
		}
		
		delete designMap;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int AnosimCommand::process(vector<SharedRAbundVector*> thisLookup) {
	try{
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1){
			driver(0, namesOfGroupCombos.size(), thisLookup, "");
		}else{
			int process = 1;
			vector<int> processIDS;
			
			//loop through and create all the processes you want
			while (process != processors) {
				int pid = fork();
				
				if (pid > 0) {
					processIDS.push_back(pid);  
					process++;
				}else if (pid == 0){
					driver(lines[process].start, lines[process].num, thisLookup, toString(getpid()));
					exit(0);
				}else { 
					m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
					for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
					exit(0);
				}
			}
			
			//do my part
			driver(lines[0].start, lines[0].num, thisLookup, "");
			
			//force parent to wait until all the processes are done
			for (int i=0;i<(processors-1);i++) { 
				int temp = processIDS[i];
				wait(&temp);
			}
			
			//append files
			for(int i = 0 ; i < calculators.size(); i++) {
				string outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + calculators[i]->getName() + ".anosim";				
				
				for (int j = 0; j < processIDS.size(); j++) {
					m->appendFiles((outputFile + toString(processIDS[j])), outputFile);
					remove((outputFile + toString(processIDS[j])).c_str());
				}
			}
		}
#else
		driver(0, namesOfGroupCombos.size(), thisLookUp, "");
#endif
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************

int AnosimCommand::driver(int start, int num, vector<SharedRAbundVector*> thisLookup, string pidValue) {
	try {
		vector<SharedRAbundVector*> subset;
		EstOutput data;
		
		//for each calculator
		for(int i = 0 ; i < calculators.size(); i++) {
			
			//create a new filename
			ofstream out;
			string outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + calculators[i]->getName() + ".anosim" + pidValue;				
			m->openOutputFileAppend(outputFile, out);
			out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
			cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
			
			//for each combo
			for (int c = start; c < (start+num); c++) {
				
				if (m->control_pressed) { out.close(); return 0; }
				
				//get set names
				vector<string> setNames;
				for (int j = 0; j < namesOfGroupCombos[c].size(); j++) { setNames.push_back(namesOfGroupCombos[c][j]); }
				
				vector<SharedRAbundVector*> thisCombosLookup;
				vector<string> thisCombosLookupSets; //what set each sharedRabund is from to be used when calculating SSWithin
				for (int k = 0; k < thisLookup.size(); k++) {
					string thisGroup = thisLookup[k]->getGroup();
					
					//is this group for a set we want to compare??
					if (m->inUsersGroups(designMap->getGroup(thisGroup), setNames)) {  
						thisCombosLookup.push_back(thisLookup[k]);
						thisCombosLookupSets.push_back(designMap->getGroup(thisGroup));
					}
					
				}
				
				int numGroups = thisCombosLookup.size();
				
				//calc the distance matrix
				matrix.clear();
				matrix.resize(numGroups);
				for (int k = 0; k < matrix.size(); k++)	{	for (int j = 0; j < matrix.size(); j++)	{	matrix[k].push_back(1.0);	}	}
				
				if (thisCombosLookup.size() == 0)  { 
					m->mothurOut("[ERROR]: Missing shared info for sets. Skipping comparison."); m->mothurOutEndLine(); 
				}else{
					
					out << thisLookup[0]->getLabel() << '\t';
					m->mothurOut(thisLookup[0]->getLabel() + "\t");
					if (setNames.size() == 2) { out << setNames[0] << "-" << setNames[1] << '\t'; m->mothurOut(setNames[0] + "-" + setNames[1] + "\t");  }
					else { out << "all" << '\t'; m->mothurOut("all\t"); }
					
					for (int k = 0; k < thisCombosLookup.size(); k++) { 
						for (int l = k; l < thisCombosLookup.size(); l++) {
							
							if (m->control_pressed) { out.close(); return 0; }
							
							if (k != l) { //we dont need to similiarity of a groups to itself
								//get estimated similarity between 2 groups
								subset.clear(); //clear out old pair of sharedrabunds
								//add new pair of sharedrabunds
								subset.push_back(thisCombosLookup[k]); subset.push_back(thisCombosLookup[l]); 
								
								//if this calc needs all groups to calculate the pair load all groups
								if (calculators[i]->getNeedsAll()) { 
									//load subset with rest of lookup for those calcs that need everyone to calc for a pair
									for (int w = 0; w < thisCombosLookup.size(); w++) {
										if ((w != k) && (w != l)) { subset.push_back(thisCombosLookup[w]); }
									}
								}
								
								data = calculators[i]->getValues(subset); //saves the calculator outputs
								
								//save values in similarity matrix
								matrix[k][l] = 1.0 - data[0];
								matrix[l][k] = 1.0 - data[0];
							}
						}
					}
					
					//calc anosim
					calcAnosim(out, setNames.size(), thisCombosLookupSets);
				}
			}
			
			out.close();
		}		
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************

int AnosimCommand::driver(int start, int num, vector<string> names, string pidValue, vector< vector<double> >& completeMatrix) {
	try {
		
		//create a new filename
		ofstream out;
		string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipfile)) + "anosim" + pidValue;				
		m->openOutputFileAppend(outputFile, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		//for each combo
		for (int c = start; c < (start+num); c++) {
			
			if (m->control_pressed) { out.close(); return 0; }
			
			//get set names
			vector<string> setNames;
			for (int j = 0; j < namesOfGroupCombos[c].size(); j++) { setNames.push_back(namesOfGroupCombos[c][j]); }
			
			vector<string> thisCombosSets; //what set each line in the distance matrix is from to be used when calculating SSWithin
			set<int> indexes;
			for (int k = 0; k < names.size(); k++) {
				//is this group for a set we want to compare??
				if (m->inUsersGroups(designMap->getGroup(names[k]), setNames)) {  
					thisCombosSets.push_back(designMap->getGroup(names[k]));	
					indexes.insert(k); //save indexes of valid rows in matrix for submatrix
				}
			}
			
			int numGroups = thisCombosSets.size();
			
			//calc the distance matrix
			matrix.clear();
			matrix.resize(numGroups);
			
			for (int k = 0; k < matrix.size(); k++)	{	for (int j = 0; j < matrix.size(); j++)	{	matrix[k].push_back(1.0);	}	}
			
			if (thisCombosSets.size() == 0)  { 
				m->mothurOut("[ERROR]: Missing distance info for sets. Skipping comparison."); m->mothurOutEndLine(); 
			}else{
				
				if (setNames.size() == 2) { out << setNames[0] << "-" << setNames[1] << '\t'; m->mothurOut(setNames[0] + "-" + setNames[1] + "\t");  }
				else { out << "all" << '\t'; m->mothurOut("all\t"); }
				
				//fill submatrix
				int rowCount = 0;
				for (int j = 0; j < completeMatrix.size(); j++) {
					
					if (indexes.count(j) != 0) { //we want at least part of this row
						int columnCount = 0;
						for (int k = 0; k < completeMatrix[j].size(); k++) {
							
							if (indexes.count(k) != 0) { //we want this distance
								matrix[rowCount][columnCount] = completeMatrix[j][k];
								matrix[columnCount][rowCount] = completeMatrix[j][k];
								columnCount++;
							}
						}
						rowCount++;
					}
				}
				
				
				//calc anosim
				calcAnosim(out, setNames.size(), thisCombosSets);
			}
		}
		
		out.close();
		
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
int AnosimCommand::calcAnosim(ofstream& out, int numTreatments, vector<string> thisCombosLookupSets) {
	try {
		//rank distances
		vector<seqDist> rankMatrix = convertToRanks();
		
		double meanBetween, meanWithin, RValue, pValue;
		
		meanWithin = calcWithinBetween(rankMatrix, thisCombosLookupSets, meanBetween);
		
		int N = matrix.size();
		int div = (N * (N-1) / 4);
	
		//calc RValue
		RValue = (meanBetween - meanWithin) / (double) div;
		
		//calc Pvalue
		int count = 0;
		for (int i = 0; i < iters; i++) {
			if (m->control_pressed) { break; }
			
			//randomly shuffle names to randomize the matrix
			vector<string> copyNames = thisCombosLookupSets;
			random_shuffle(copyNames.begin(), copyNames.end());
			
			meanWithin = calcWithinBetween(rankMatrix, copyNames, meanBetween);
			
			//calc RValue
			double randomRValue = (meanBetween - meanWithin) / (double) div;	
			
			if (randomRValue >= RValue) { count++; }
		}
		
		pValue = count / (float) iters;
		
		out << RValue << '\t' << pValue << endl;
		cout << RValue << '\t' << pValue << endl;
		m->mothurOutJustToLog(toString(RValue) + "\t" + toString(pValue) + "\n");
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "calcAnisom");
		exit(1);
	}
}
//**********************************************************************************************************************
double AnosimCommand::calcWithinBetween(vector<seqDist>& thisMatrix, vector<string> thisCombosLookupSets, double& between) {
	try {
		double within = 0.0;
		int count = 0;
		int count2 = 0;
		between = 0.0;
		
		for (int l = 0; l < thisMatrix.size(); l++) {
			//if you are from the same treatment 
			if (thisCombosLookupSets[thisMatrix[l].seq1] == thisCombosLookupSets[thisMatrix[l].seq2]) { 
				within += thisMatrix[l].dist; //rank of this distance
				count++;
			}else { //different treatments
				between += thisMatrix[l].dist; //rank of this distance
				count2++;
			}
		}
		
		within /= (float) count;
		between /= (float) count2;
		
		return within;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "calcWithinBetween");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<seqDist> AnosimCommand::convertToRanks() {
	try {
		vector<seqDist> ranks;
		
		for (int i = 0; i < matrix.size(); i++) {
			for (int j = 0; j < i; j++) {
				seqDist member(i, j, matrix[i][j]);
				ranks.push_back(member);
			}
		}
		
		//sort distances
		sort(ranks.begin(), ranks.end(), compareSequenceDistance); 	
		
		//find ranks of distances
		vector<seqDist*> ties;
		int rankTotal = 0;
		for (int j = 0; j < ranks.size(); j++) {
			rankTotal += (j+1);
			ties.push_back(&ranks[j]);
			
			if (j != (ranks.size()-1)) { // you are not the last so you can look ahead
				if (ranks[j].dist != ranks[j+1].dist) { // you are done with ties, rank them and continue
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						(*ties[k]).dist = thisrank;
					}
					ties.clear();
					rankTotal = 0;
				}
			}else { // you are the last one
				
				for (int k = 0; k < ties.size(); k++) {
					float thisrank = rankTotal / (float) ties.size();
					(*ties[k]).dist = thisrank;
				}
			}
		}
		
		return ranks;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "convertToRanks");
		exit(1);
	}
}
//**********************************************************************************************************************



