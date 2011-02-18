/*
 *  amovacommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "amovacommand.h"
#include "sharedutilities.h"
//#include "sharedsobscollectsummary.h"
//#include "sharedchao1.h"
//#include "sharedace.h"
//#include "sharednseqs.h"
//#include "sharedjabund.h"
//#include "sharedsorabund.h"
//#include "sharedjclass.h"
//#include "sharedsorclass.h"
//#include "sharedjest.h"
//#include "sharedsorest.h"
//#include "sharedthetayc.h"
//#include "sharedthetan.h"
//#include "sharedkstest.h"
//#include "whittaker.h"
//#include "sharedochiai.h"
//#include "sharedanderbergs.h"
//#include "sharedkulczynski.h"
//#include "sharedkulczynskicody.h"
//#include "sharedlennon.h"
//#include "sharedmorisitahorn.h"
//#include "sharedbraycurtis.h"
//#include "sharedjackknife.h"
//#include "whittaker.h"
//#include "odum.h"
//#include "canberra.h"
//#include "structeuclidean.h"
//#include "structchord.h"
//#include "hellinger.h"
//#include "manhattan.h"
//#include "structpearson.h"
//#include "soergel.h"
//#include "spearman.h"
//#include "structkulczynski.h"
//#include "structchi2.h"
//#include "speciesprofile.h"
//#include "hamming.h"
//#include "gower.h"
//#include "memchi2.h"
//#include "memchord.h"
//#include "memeuclidean.h"
//#include "mempearson.h"

//**********************************************************************************************************************
vector<string> AmovaCommand::getValidParameters(){	
	try {
		string Array[] =  {"groups","label","outputdir","iters","phylip","design","alpha", 
			//"sets","processors"
			"inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
AmovaCommand::AmovaCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["amova"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "AmovaCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AmovaCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"design"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AmovaCommand::getRequiredFiles(){	
	try {
		string Array[] =  {};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

AmovaCommand::AmovaCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false; calledHelp = false;   
//		allLines = 1;
//		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label","outputdir","iters","phylip","alpha",
//				"design","sets","processors",
				"inputdir"};
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
			outputTypes["amova"] = tempOutNames;
			
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
			
			phylipFileName = validParameter.validFile(parameters, "phylip", true);
			if (phylipFileName == "not open") { phylipFileName = ""; abort = true; }
			else if (phylipFileName == "not found") { phylipFileName = ""; }	
			else {
				globaldata->newRead();
				globaldata->setPhylipFile(phylipFileName);
			}
			
			//check for required parameters
			designFileName = validParameter.validFile(parameters, "design", true);
			if (designFileName == "not open") { abort = true; }
			else if (designFileName == "not found") {
				designFileName = "";
				m->mothurOut("You must provide an design file.");
				m->mothurOutEndLine();
				abort = true;
			}	

			string temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "alpha", false);			if (temp == "not found") { temp = "0.05"; }
			convert(temp, experimentwiseAlpha); 
			
//			make sure the user has already run the read.otu command
//			if ((globaldata->getSharedFile() == "")) {
//				if ((phylipfile == "")) {
//					m->mothurOut("You must read a list and a group, a shared file or provide a distance matrix before you can use the amova command.");
//					m->mothurOutEndLine();
//					abort = true; 
//				}
//			}else { sharedfile = globaldata->getSharedFile(); } 
//				
//			//use distance matrix if one is provided
//			if ((sharedfile != "") && (phylipfile != "")) { sharedfile = ""; }
//			
//			//check for optional parameter and set defaults
//			// ...at some point should added some additional type checking...
//			label = validParameter.validFile(parameters, "label", false);			
//			if (label == "not found") { label = ""; }
//			else { 
//				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
//				else { allLines = 1;  }
//			}
//			
//			//if the user has not specified any labels use the ones from read.otu
//			if (label == "") {  
//				allLines = globaldata->allLines; 
//				labels = globaldata->labels; 
//			}
//			
//			groups = validParameter.validFile(parameters, "groups", false);			
//			if (groups == "not found") { groups = "";  }
//			else { 
//				m->splitAtDash(groups, Groups);
//				globaldata->Groups = Groups;
//			}
//			
//			sets = validParameter.validFile(parameters, "sets", false);			
//			if (sets == "not found") { sets = ""; pickedGroups = false; }
//			else { 
//				pickedGroups = true;
//				m->splitAtDash(sets, Sets);
//			}
//		
//			temp = validParameter.validFile(parameters, "processors", false);			if (temp == "not found"){	temp = "1";	}
//			convert(temp, processors); 
//			
//			vector<string> Estimators;
//			calc = validParameter.validFile(parameters, "calc", false);			
//			if (calc == "not found") { calc = "morisitahorn";  }
//			m->splitAtDash(calc, Estimators);
//				
//			if (abort == false) {
//				
//				ValidCalculators* validCalculator = new ValidCalculators();
//				
//				for (int i=0; i<Estimators.size(); i++) {
//					if (validCalculator->isValidCalculator("treegroup", Estimators[i]) == true) { 
//						if (Estimators[i] == "sharedsobs") { 
//							calculators.push_back(new SharedSobsCS());
//						}else if (Estimators[i] == "sharedchao") { 
//							calculators.push_back(new SharedChao1());
//						}else if (Estimators[i] == "sharedace") { 
//							calculators.push_back(new SharedAce());
//						}else if (Estimators[i] == "jabund") { 	
//							calculators.push_back(new JAbund());
//						}else if (Estimators[i] == "sorabund") { 
//							calculators.push_back(new SorAbund());
//						}else if (Estimators[i] == "jclass") { 
//							calculators.push_back(new Jclass());
//						}else if (Estimators[i] == "sorclass") { 
//							calculators.push_back(new SorClass());
//						}else if (Estimators[i] == "jest") { 
//							calculators.push_back(new Jest());
//						}else if (Estimators[i] == "sorest") { 
//							calculators.push_back(new SorEst());
//						}else if (Estimators[i] == "thetayc") { 
//							calculators.push_back(new ThetaYC());
//						}else if (Estimators[i] == "thetan") { 
//							calculators.push_back(new ThetaN());
//						}else if (Estimators[i] == "kstest") { 
//							calculators.push_back(new KSTest());
//						}else if (Estimators[i] == "sharednseqs") { 
//							calculators.push_back(new SharedNSeqs());
//						}else if (Estimators[i] == "ochiai") { 
//							calculators.push_back(new Ochiai());
//						}else if (Estimators[i] == "anderberg") { 
//							calculators.push_back(new Anderberg());
//						}else if (Estimators[i] == "kulczynski") { 
//							calculators.push_back(new Kulczynski());
//						}else if (Estimators[i] == "kulczynskicody") { 
//							calculators.push_back(new KulczynskiCody());
//						}else if (Estimators[i] == "lennon") { 
//							calculators.push_back(new Lennon());
//						}else if (Estimators[i] == "morisitahorn") { 
//							calculators.push_back(new MorHorn());
//						}else if (Estimators[i] == "braycurtis") { 
//							calculators.push_back(new BrayCurtis());
//						}else if (Estimators[i] == "whittaker") { 
//							calculators.push_back(new Whittaker());
//						}else if (Estimators[i] == "odum") { 
//							calculators.push_back(new Odum());
//						}else if (Estimators[i] == "canberra") { 
//							calculators.push_back(new Canberra());
//						}else if (Estimators[i] == "structeuclidean") { 
//							calculators.push_back(new StructEuclidean());
//						}else if (Estimators[i] == "structchord") { 
//							calculators.push_back(new StructChord());
//						}else if (Estimators[i] == "hellinger") { 
//							calculators.push_back(new Hellinger());
//						}else if (Estimators[i] == "manhattan") { 
//							calculators.push_back(new Manhattan());
//						}else if (Estimators[i] == "structpearson") { 
//							calculators.push_back(new StructPearson());
//						}else if (Estimators[i] == "soergel") { 
//							calculators.push_back(new Soergel());
//						}else if (Estimators[i] == "spearman") { 
//							calculators.push_back(new Spearman());
//						}else if (Estimators[i] == "structkulczynski") { 
//							calculators.push_back(new StructKulczynski());
//						}else if (Estimators[i] == "speciesprofile") { 
//							calculators.push_back(new SpeciesProfile());
//						}else if (Estimators[i] == "hamming") { 
//							calculators.push_back(new Hamming());
//						}else if (Estimators[i] == "structchi2") { 
//							calculators.push_back(new StructChi2());
//						}else if (Estimators[i] == "gower") { 
//							calculators.push_back(new Gower());
//						}else if (Estimators[i] == "memchi2") { 
//							calculators.push_back(new MemChi2());
//						}else if (Estimators[i] == "memchord") { 
//							calculators.push_back(new MemChord());
//						}else if (Estimators[i] == "memeuclidean") { 
//							calculators.push_back(new MemEuclidean());
//						}else if (Estimators[i] == "mempearson") { 
//							calculators.push_back(new MemPearson());
//						}
//					}
//				}
//			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "AmovaCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void AmovaCommand::help(){
	try {
		m->mothurOut("Referenced: Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecol 26: 32-46.\n");
//		m->mothurOut("The amova command can only be executed after a successful read.otu command of a list and group or shared file, or by providing a phylip formatted distance matrix.\n");
		m->mothurOut("The amova command outputs a .amova file. \n");
		m->mothurOut("The amova command parameters are phylip, iters, and alpha.  The phylip and design parameters are required.\n");
		m->mothurOut("The design parameter allows you to assign your samples to groups when you are running amova. It is required. \n");
		m->mothurOut("The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the sample name and the second column is the group the sample belongs to.\n");
//		m->mothurOut("The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile. To run the pairwise comparisons of all sets use sets=all.\n");
		m->mothurOut("The iters parameter allows you to set number of randomization for the P value.  The default is 1000. \n");
//		m->mothurOut("The groups parameter allows you to specify which of the groups you would like included. The group names are separated by dashes. groups=all will find all pairwise comparisons. \n");
//		m->mothurOut("The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n");
//		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		m->mothurOut("The amova command should be in the following format: amova(phylip=file.dist, design=file.design).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e. 1000).\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

AmovaCommand::~AmovaCommand(){}

//**********************************************************************************************************************

int AmovaCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read design file
		designMap = new GroupMap(designFileName);
		designMap->readDesignMap();

		if (outputDir == "") { outputDir = m->hasPath(phylipFileName); }
						
		//read in distance matrix and square it
		ReadPhylipVector readMatrix(phylipFileName);
		vector<string> sampleNames = readMatrix.read(distanceMatrix);
		int numSamples = sampleNames.size();
		
		for(int i=0;i<distanceMatrix.size();i++){
			for(int j=0;j<i;j++){
				distanceMatrix[i][j] *= distanceMatrix[i][j];	
			}
		}
		
		//link designMap to rows/columns in distance matrix
		map<string, vector<int> > origGroupSampleMap;
		for(int i=0;i<numSamples;i++){
			origGroupSampleMap[designMap->getGroup(sampleNames[i])].push_back(i);
		}
		int numGroups = origGroupSampleMap.size();
		
		//create a new filename
		ofstream AMOVAFile;
		string AMOVAFileName = outputDir + m->getRootName(m->getSimpleName(phylipFileName))  + "amova";				
		m->openOutputFile(AMOVAFileName, AMOVAFile);
		outputNames.push_back(AMOVAFileName); outputTypes["amova"].push_back(AMOVAFileName);

		experimentwiseAlpha = 0.05;
		
		double fullANOVAPValue = runAMOVA(AMOVAFile, origGroupSampleMap, experimentwiseAlpha);
		if(fullANOVAPValue <= experimentwiseAlpha && numGroups > 2){
			
			int numCombos = numGroups * (numGroups-1) / 2;
			double pairwiseAlpha = experimentwiseAlpha / (double) numCombos;
			
			map<string, vector<int> >::iterator itA;
			map<string, vector<int> >::iterator itB;
			
			for(itA=origGroupSampleMap.begin();itA!=origGroupSampleMap.end();itA++){
				itB = itA;itB++;
				for(itB;itB!=origGroupSampleMap.end();itB++){
					
					map<string, vector<int> > pairwiseGroupSampleMap;
					pairwiseGroupSampleMap[itA->first] = itA->second;
					pairwiseGroupSampleMap[itB->first] = itB->second;
					
					runAMOVA(AMOVAFile, pairwiseGroupSampleMap, pairwiseAlpha);
				}			
			}
			cout << "Experiment-wise error rate: " << experimentwiseAlpha << endl;
			cout << "Pair-wise error rate (Bonferroni): " << pairwiseAlpha << endl;
		}
		else{
			cout << "Experiment-wise error rate: " << experimentwiseAlpha << endl;
		}

		cout << "If you have borderline P-values, you should try increasing the number of iterations" << endl;

		
//		//make sure sets are all in designMap
//		SharedUtil* util = new SharedUtil();  
//		util->setGroups(Sets, designMap->namesOfGroups);  
//		delete util;
//		
//		//if the user pickedGroups or set sets=all, then we want to figure out how to split the pairwise comparison between processors
//		if (pickedGroups) {
//			for (int a=0; a<numGroups; a++) { 
//				for (int l = 0; l < a; l++) {	
//					vector<string> groups; groups.push_back(Sets[a]); groups.push_back(Sets[l]);
//					namesOfGroupCombos.push_back(groups);
//				}
//			}
//		}else { //one giant compare
//			vector<string> groups;
//			for (int a=0; a<numGroups; a++) { groups.push_back(Sets[a]); }
//			namesOfGroupCombos.push_back(groups);
//		}
//		
//		//only 1 combo
//		if (numGroups == 2) { processors = 1; }
//		else if (numGroups < 2)	{
//			m->mothurOut("Not enough sets, I need at least 2 valid sets. Unable to complete command.");
//			m->mothurOutEndLine(); m->control_pressed = true;
//		}
//		
//		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
//			if(processors != 1){
//				int numPairs = namesOfGroupCombos.size();
//				int numPairsPerProcessor = numPairs / processors;
//			
//				for (int i = 0; i < processors; i++) {
//					int startPos = i * numPairsPerProcessor;
//					if(i == processors - 1){
//						numPairsPerProcessor = numPairs - i * numPairsPerProcessor;
//					}
//					lines.push_back(linePair(startPos, numPairsPerProcessor));
//				}
//			}
//		#endif
//		
//		if (sharedfile != "") { //create distance matrix for each label
//			
//			if (outputDir == "") { outputDir = m->hasPath(sharedfile); }
//			
//			//for each calculator
//			for(int i = 0 ; i < calculators.size(); i++) {
//				
//				//create a new filename
//				ofstream out;
//				string outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + calculators[i]->getName() + ".amova";
//				m->openOutputFile(outputFile, out);
//				outputNames.push_back(outputFile); outputTypes["amova"].push_back(outputFile);
//				
//				//print headers
//				out << "label\tgroupsCompared\tMeanSquaredWithin\tMeanSquaredAmong\tFstatistic\tpValue" << endl;  
//				out.close();
//				m->mothurOut("label\tgroupsCompared\tMeanSquaredWithin\tMeanSquaredAmong\tFstatistic\tpValue\n");  
//			}
//			
//			InputData input(sharedfile, "sharedfile");
//			vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
//			string lastLabel = lookup[0]->getLabel();
//			
//			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
//			set<string> processedLabels;
//			set<string> userLabels = labels;
//			
//			//sanity check between shared file groups and design file groups
//			for (int i = 0; i < lookup.size(); i++) { 
//				string group = designMap->getGroup(lookup[i]->getGroup());
//				
//				if (group == "not found") { m->control_pressed = true;  m->mothurOut("[ERROR]: " + lookup[i]->getGroup() + " is not in your design file, please correct."); m->mothurOutEndLine();  }
//			}
//			
//			//as long as you are not at the end of the file or done wih the lines you want
//			while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
//				
//				if (m->control_pressed) {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } globaldata->Groups.clear();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
//				
//				if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
//					
//					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
//					process(lookup);
//					
//					processedLabels.insert(lookup[0]->getLabel());
//					userLabels.erase(lookup[0]->getLabel());
//				}
//				
//				if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
//					string saveLabel = lookup[0]->getLabel();
//					
//					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
//					lookup = input.getSharedRAbundVectors(lastLabel);
//					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
//					
//					process(lookup);
//					
//					processedLabels.insert(lookup[0]->getLabel());
//					userLabels.erase(lookup[0]->getLabel());
//					
//					//restore real lastlabel to save below
//					lookup[0]->setLabel(saveLabel);
//				}
//				
//				lastLabel = lookup[0]->getLabel();
//				//prevent memory leak
//				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
//				
//				if (m->control_pressed) {  globaldata->Groups.clear();   delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
//				
//				//get next line to process
//				lookup = input.getSharedRAbundVectors();				
//			}
//			
//			if (m->control_pressed) { globaldata->Groups.clear();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }  return 0; }
//			
//			//output error messages about any remaining user labels
//			set<string>::iterator it;
//			bool needToRun = false;
//			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
//				m->mothurOut("Your file does not include the label " + *it); 
//				if (processedLabels.count(lastLabel) != 1) {
//					m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
//					needToRun = true;
//				}else {
//					m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
//				}
//			}
//			
//			//run last label if you need to
//			if (needToRun == true)  {
//				for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
//				lookup = input.getSharedRAbundVectors(lastLabel);
//				
//				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
//				
//				process(lookup);
//				
//				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
//			}
//			
//			//reset groups parameter
//			globaldata->Groups.clear();  
//			
//			
//		}else { //user provided distance matrix
//			
//			if (outputDir == "") { outputDir = m->hasPath(phylipFile); }
//					
//			//create a new filename
//			ofstream out;
//			string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipFile))  + "amova";				
//			m->openOutputFile(outputFile, out);
//			outputNames.push_back(outputFile); outputTypes["amova"].push_back(outputFile);
//			
//			//print headers
//			out << "groupsCompared\tMeanSquaredWithin\tMeanSquaredAmong\tFstatistic\tpValue" << endl;  
//			out.close();
//			m->mothurOut("groupsCompared\tMeanSquaredWithin\tMeanSquaredAmong\tFstatistic\tpValue\n");  
//			
//			ReadPhylipVector readMatrix(phylipFile);
//			vector< vector<double> > completeMatrix;
//			vector<string> names = readMatrix.read(completeMatrix);
//			
//			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
//				if(processors == 1){
//					driver(0, namesOfGroupCombos.size(), names, "", completeMatrix);
//				}else{
//					int process = 1;
//					vector<int> processIDS;
//					
//					//loop through and create all the processes you want
//					while (process != processors) {
//						int pid = fork();
//						
//						if (pid > 0) {
//							processIDS.push_back(pid);  
//							process++;
//						}else if (pid == 0){
//							driver(lines[process].start, lines[process].num, names, toString(getpid()), completeMatrix);
//							exit(0);
//						}else { 
//							m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
//							for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
//							exit(0);
//						}
//					}
//					
//					//do my part
//					driver(lines[0].start, lines[0].num, names, "", completeMatrix);
//					
//					//force parent to wait until all the processes are done
//					for (int i=0;i<(processors-1);i++) { 
//						int temp = processIDS[i];
//						wait(&temp);
//					}
//					
//					//append files
//					string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipFile)) + "amova";				
//					for (int j = 0; j < processIDS.size(); j++) {
//						m->appendFiles((outputFile + toString(processIDS[j])), outputFile);
//						remove((outputFile + toString(processIDS[j])).c_str());
//					}
//					
//				}
//			#else
//				driver(0, namesOfGroupCombos.size(), names, "", completeMatrix);
//			#endif
//			
//		}
		
		delete designMap;
	 
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

double AmovaCommand::runAMOVA(ofstream& AMOVAFile, map<string, vector<int> > groupSampleMap, double alpha) {
	try {
		map<string, vector<int> >::iterator it;

		int numGroups = groupSampleMap.size();
		int numSamples = 0;

		for(it = groupSampleMap.begin();it!=groupSampleMap.end();it++){
			numSamples += it->second.size();			
		}

		double ssTotalOrig = calcSSTotal(groupSampleMap);
		double ssWithinOrig = calcSSWithin(groupSampleMap);
		double ssAmongOrig = ssTotalOrig - ssWithinOrig;
		
		double counter = 0;
		for(int i=0;i<iters;i++){
			map<string, vector<int> > randomizedGroup = getRandomizedGroups(groupSampleMap);
			double ssWithinRand = calcSSWithin(randomizedGroup);
			if(ssWithinRand < ssWithinOrig){	counter++;	}
		}
		
		double pValue = (double)counter / (double) iters;
		string pString = "";
		if(pValue < 1/(double)iters){	pString = '<' + toString(1/(double)iters);	}
		else						{	pString = toString(pValue);					}
		
		
		//print anova table
		it = groupSampleMap.begin();
		cout << it->first;
		it++;
		for(it;it!=groupSampleMap.end();it++){
			cout << " vs. " << it->first;
		}cout << endl;
		cout << "ANOVA\tAmong\tWithin\tTotal" << endl;
		cout << "SS\t" << ssAmongOrig << '\t' << ssWithinOrig << '\t' << ssTotalOrig << endl;
		
		int dfAmong = numGroups - 1;			double MSAmong = ssAmongOrig / (double) dfAmong;
		int dfWithin = numSamples - numGroups;	double MSWithin = ssWithinOrig / (double) dfWithin;
		int dfTotal = numSamples - 1;			double Fs = MSAmong / MSWithin;
		
		cout << "df\t" << dfAmong << '\t' << dfWithin << '\t' << dfTotal << endl;
		cout << "MS\t" << MSAmong << '\t' << MSWithin << endl << endl;
		cout << "Fs:\t" << Fs << endl;
		cout << "p-value: " << pString;
		if(pValue < alpha){	cout << " ***";	}
		cout << endl << endl;

		return pValue;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "calcAmova");
		exit(1);
	}
}

//**********************************************************************************************************************

map<string, vector<int> > AmovaCommand::getRandomizedGroups(map<string, vector<int> > origMapping){
	try{
		vector<int> sampleIndices;
		vector<int> samplesPerGroup;
		
		map<string, vector<int> >::iterator it;
		for(it=origMapping.begin();it!=origMapping.end();it++){
			vector<int> indices = it->second;
			samplesPerGroup.push_back(indices.size());
			sampleIndices.insert(sampleIndices.end(), indices.begin(), indices.end());
		}
		
		random_shuffle(sampleIndices.begin(), sampleIndices.end());
		
		int index = 0;
		map<string, vector<int> > randomizedGroups = origMapping;
		for(it=randomizedGroups.begin();it!=randomizedGroups.end();it++){
			for(int i=0;i<it->second.size();i++){
				it->second[i] = sampleIndices[index++];				
			}
		}

		return randomizedGroups;		
	}
	catch (exception& e) {
		m->errorOut(e, "AmovaCommand", "randomizeGroups");
		exit(1);
	}
}

//**********************************************************************************************************************

double AmovaCommand::calcSSTotal(map<string, vector<int> >& groupSampleMap) {
	try {
		
		vector<int> indices;
		map<string, vector<int> >::iterator it;
		for(it=groupSampleMap.begin();it!=groupSampleMap.end();it++){
			indices.insert(indices.end(), it->second.begin(), it->second.end());			
		}
		sort(indices.begin(), indices.end());
			
 		int numIndices =indices.size();
		double ssTotal = 0.0;
		
		for(int i=1;i<numIndices;i++){
			int row = indices[i];
			
			for(int j=0;j<i;j++){
				ssTotal += distanceMatrix[row][indices[j]];
			}
		}
		ssTotal /= numIndices;
			
		return ssTotal;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "calcSSTotal");
		exit(1);
	}
}

//**********************************************************************************************************************

double AmovaCommand::calcSSWithin(map<string, vector<int> >& groupSampleMap) {
	try {

		double ssWithin = 0.0;
		
		map<string, vector<int> >::iterator it;
		for(it=groupSampleMap.begin();it!=groupSampleMap.end();it++){
			
			double withinGroup = 0;
			
			vector<int> samples = it->second;
			
			for(int i=0;i<samples.size();i++){
				int row = samples[i];

				for(int j=0;j<samples.size();j++){
					int col = samples[j];

					if(col < row){
						withinGroup += distanceMatrix[row][col];
					}
					
				}
			}

			ssWithin += withinGroup / samples.size();
		}

		return ssWithin;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "calcWithin");
		exit(1);
	}
}

//**********************************************************************************************************************

//int AmovaCommand::process(vector<SharedRAbundVector*> thisLookup) {
//	try{
//		
//#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
//		if(processors == 1){
//			driver(0, namesOfGroupCombos.size(), thisLookup, "");
//		}else{
//			int process = 1;
//			vector<int> processIDS;
//			
//			//loop through and create all the processes you want
//			while (process != processors) {
//				int pid = fork();
//				
//				if (pid > 0) {
//					processIDS.push_back(pid);  
//					process++;
//				}else if (pid == 0){
//					driver(lines[process].start, lines[process].num, thisLookup, toString(getpid()));
//					exit(0);
//				}else { 
//					m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
//					for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
//					exit(0);
//				}
//			}
//			
//			//do my part
//			driver(lines[0].start, lines[0].num, thisLookup, "");
//			
//			//force parent to wait until all the processes are done
//			for (int i=0;i<(processors-1);i++) { 
//				int temp = processIDS[i];
//				wait(&temp);
//			}
//			
//			//append files
//			for(int i = 0 ; i < calculators.size(); i++) {
//				string outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + calculators[i]->getName() + ".amova";				
//
//				for (int j = 0; j < processIDS.size(); j++) {
//					m->appendFiles((outputFile + toString(processIDS[j])), outputFile);
//					remove((outputFile + toString(processIDS[j])).c_str());
//				}
//			}
//		}
//#else
//		driver(0, namesOfGroupCombos.size(), thisLookUp, "");
//#endif
//			
//		return 0;
//	}
//	catch(exception& e) {
//		m->errorOut(e, "AmovaCommand", "process");
//		exit(1);
//	}
//}

//**********************************************************************************************************************

//int AmovaCommand::driver(int start, int num, vector<SharedRAbundVector*> thisLookup, string pidValue) {
//	try {
//		vector<SharedRAbundVector*> subset;
//		EstOutput data;
//		
//		//for each calculator
//		for(int i = 0 ; i < calculators.size(); i++) {
//			
//			//create a new filename
//			ofstream out;
//			string outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + calculators[i]->getName() + ".amova" + pidValue;				
//			m->openOutputFileAppend(outputFile, out);
//			out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
//			
//			//for each combo
//			for (int c = start; c < (start+num); c++) {
//				
//				if (m->control_pressed) { out.close(); return 0; }
//				
//				//get set names
//				vector<string> setNames;
//				for (int j = 0; j < namesOfGroupCombos[c].size(); j++) { setNames.push_back(namesOfGroupCombos[c][j]); }
//				
//				vector<SharedRAbundVector*> thisCombosLookup;
//				vector<string> thisCombosLookupSets; //what set each sharedRabund is from to be used when calculating SSWithin
//				for (int k = 0; k < thisLookup.size(); k++) {
//					string thisGroup = thisLookup[k]->getGroup();
//					
//					//is this group for a set we want to compare??
//					if (m->inUsersGroups(designMap->getGroup(thisGroup), setNames)) {  
//						thisCombosLookup.push_back(thisLookup[k]);
//						thisCombosLookupSets.push_back(designMap->getGroup(thisGroup));
//					}
//					
//				}
//				
//				int numGroups = thisCombosLookup.size();
//			
//				//calc the distance matrix
//				matrix.clear();
//				matrix.resize(numGroups);
//				for (int k = 0; k < matrix.size(); k++)	{	for (int j = 0; j < matrix.size(); j++)	{	matrix[k].push_back(1.0);	}	}
//				
//				if (thisCombosLookup.size() == 0)  { 
//					m->mothurOut("[ERROR]: Missing shared info for sets. Skipping comparison."); m->mothurOutEndLine(); 
//				}else{
//					
//					m->mothurOut(thisLookup[0]->getLabel() + '\t');
//					out << thisLookup[0]->getLabel() << '\t';
//					if (setNames.size() == 2) {
//						m->mothurOut(setNames[0] + '-' + setNames[1] + '\t');
//						out << setNames[0] << "-" << setNames[1] << '\t'; 
//					}
//					else {
//						out << "all" << '\t'; 
//						m->mothurOut("all\t");
//					}
//					
//					for (int k = 0; k < thisCombosLookup.size(); k++) { 
//						for (int l = k; l < thisCombosLookup.size(); l++) {
//							
//							if (m->control_pressed) { out.close(); return 0; }
//							
//							if (k != l) { //we dont need to similiarity of a groups to itself
//								//get estimated similarity between 2 groups
//								subset.clear(); //clear out old pair of sharedrabunds
//								//add new pair of sharedrabunds
//								subset.push_back(thisCombosLookup[k]); subset.push_back(thisCombosLookup[l]); 
//								
//								//if this calc needs all groups to calculate the pair load all groups
//								if (calculators[i]->getNeedsAll()) { 
//									//load subset with rest of lookup for those calcs that need everyone to calc for a pair
//									for (int w = 0; w < thisCombosLookup.size(); w++) {
//										if ((w != k) && (w != l)) { subset.push_back(thisCombosLookup[w]); }
//									}
//								}
//								
//								data = calculators[i]->getValues(subset); //saves the calculator outputs
//								
//								//save values in similarity matrix
//								matrix[k][l] = 1.0 - data[0];
//								matrix[l][k] = 1.0 - data[0];
//							}
//						}
//					}
//					
//					//calc amova
//					calcAmova(out, setNames.size(), thisCombosLookupSets);
//				}
//			}
//			
//			out.close();
//		}		
//		
//		return 0;
//		
//	}
//	catch(exception& e) {
//		m->errorOut(e, "AmovaCommand", "driver");
//		exit(1);
//	}
//}
//
//**********************************************************************************************************************
//
//int AmovaCommand::driver(int start, int num, vector<string> names, string pidValue, vector< vector<double> >& completeMatrix) {
//	try {
//		
//		//create a new filename
//		ofstream out;
//		string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipFile)) + "amova" + pidValue;				
//		m->openOutputFileAppend(outputFile, out);
//		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
//		
//		//for each combo
//		for (int c = start; c < (start+num); c++) {
//			
//			if (m->control_pressed) { out.close(); return 0; }
//			
//			//get set names
//			vector<string> setNames;
//			for (int j = 0; j < namesOfGroupCombos[c].size(); j++) { setNames.push_back(namesOfGroupCombos[c][j]); }
//			
//			vector<string> thisCombosSets; //what set each line in the distance matrix is from to be used when calculating SSWithin
//			set<int> indexes;
//			for (int k = 0; k < names.size(); k++) {
//				//is this group for a set we want to compare??
//				if (m->inUsersGroups(designMap->getGroup(names[k]), setNames)) {  
//					thisCombosSets.push_back(designMap->getGroup(names[k]));	
//					indexes.insert(k); //save indexes of valid rows in matrix for submatrix
//				}
//			}
//			
//			int numGroups = thisCombosSets.size();
//			
//			//calc the distance matrix
//			matrix.clear();
//			matrix.resize(numGroups);
//			
//			for (int k = 0; k < matrix.size(); k++)	{	for (int j = 0; j < matrix.size(); j++)	{	matrix[k].push_back(1.0);	}	}
//			
//			if (thisCombosSets.size() == 0)  { 
//				m->mothurOut("[ERROR]: Missing distance info for sets. Skipping comparison."); m->mothurOutEndLine(); 
//			}else{
//				
//				if (setNames.size() == 2) { 
//					out << setNames[0] << "-" << setNames[1] << '\t'; 
//					m->mothurOut(setNames[0] + '-' + setNames[1] + '\t');
//				}
//				else {
//					out << "all" << '\t'; 
//					m->mothurOut("all\t");
//				}
//				
//				//fill submatrix
//				int rowCount = 0;
//				for (int j = 0; j < completeMatrix.size(); j++) {
//					
//					if (indexes.count(j) != 0) { //we want at least part of this row
//						int columnCount = 0;
//						for (int k = 0; k < completeMatrix[j].size(); k++) {
//							
//							if (indexes.count(k) != 0) { //we want this distance
//								matrix[rowCount][columnCount] = completeMatrix[j][k];
//								matrix[columnCount][rowCount] = completeMatrix[j][k];
//								columnCount++;
//							}
//						}
//						rowCount++;
//					}
//				}
//				
//				//calc amova
//				calcAmova(out, setNames.size(), thisCombosSets);
//			}
//		}
//		
//		out.close();
//		
//	
//		return 0;
//		
//	}
//	catch(exception& e) {
//		m->errorOut(e, "AmovaCommand", "driver");
//		exit(1);
//	}
//}
//
//**********************************************************************************************************************

