/*
 *  summarycommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarycommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "simpsoneven.h"
#include "invsimpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "heip.h"
#include "smithwilson.h"
#include "shannoneven.h"
#include "jackknife.h"
#include "geom.h"
#include "logsd.h"
#include "qstat.h"
#include "bergerparker.h"
#include "bstick.h"
#include "goodscoverage.h"
#include "coverage.h"
#include "efron.h"
#include "boneh.h"
#include "solow.h"
#include "shen.h"

//**********************************************************************************************************************
vector<string> SummaryCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pcalc("calc", "Multiple", "sobs-chao-nseqs-coverage-ace-jack-shannon-shannoneven-np_shannon-heip-smithwilson-simpson-simpsoneven-invsimpson-bootstrap-geometric-qstat-logseries-bergerparker-bstick-goodscoverage-efron-boneh-solow-shen", "sobs-chao-ace-jack-shannon-npshannon-simpson", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pabund("abund", "Number", "", "10", "", "", "",false,false); parameters.push_back(pabund);
		CommandParameter psize("size", "Number", "", "0", "", "", "",false,false); parameters.push_back(psize);
		CommandParameter pgroupmode("groupmode", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pgroupmode);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The summary.single command parameters are list, sabund, rabund, shared, label, calc, abund and groupmode.  list, sabund, rabund or shared is required unless you have a valid current file.\n";
		helpString += "The summary.single command should be in the following format: \n";
		helpString += "summary.single(label=yourLabel, calc=yourEstimators).\n";
		helpString += "Example summary.single(label=unique-.01-.03, calc=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson).\n";
		helpString += validCalculator.printCalc("summary");
		helpString += "The default value calc is sobs-chao-ace-jack-shannon-npshannon-simpson\n";
		helpString += "If you are running summary.single with a shared file and would like your summary results collated in one file, set groupmode=t. (Default=true).\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabels).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
SummaryCommand::SummaryCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "SummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

SummaryCommand::SummaryCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
				
		//allow user to run help
		if(option == "help") {  help();  abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; m->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; m->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; m->setRabundFile(rabundfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; m->setSharedFile(sharedfile); }
			
			if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					listfile = m->getListFile(); 
					if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else { 
						rabundfile = m->getRabundFile(); 
						if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
						else { 
							sabundfile = m->getSabundFile(); 
							if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
							else { 
								m->mothurOut("No valid current files. You must provide a list, sabund, rabund or shared file before you can use the collect.single command."); m->mothurOutEndLine(); 
								abort = true;
							}
						}
					}
				}
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile);		}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
				 if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}

			string temp;
			temp = validParameter.validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "size", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, size); 
			
			temp = validParameter.validFile(parameters, "groupmode", false);		if (temp == "not found") { temp = "T"; }
			groupMode = m->isTrue(temp);
			
	
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "SummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SummaryCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if ((format != "sharedfile")) { inputFileNames.push_back(inputfile);  }
		else {  inputFileNames = parseSharedFile(sharedfile);  format = "rabund"; }
		
		if (m->control_pressed) { return 0; }
		
		int numLines = 0;
		int numCols = 0;
		
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			numLines = 0;
			numCols = 0;
			
			string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(inputFileNames[p])) + "summary";
			outputNames.push_back(fileNameRoot); outputTypes["summary"].push_back(fileNameRoot);
			
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
			
			sumCalculators.clear();
			
			ValidCalculators validCalculator;
			
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("summary", Estimators[i]) == true) { 
					if(Estimators[i] == "sobs"){
						sumCalculators.push_back(new Sobs());
					}else if(Estimators[i] == "chao"){
						sumCalculators.push_back(new Chao1());
					}else if(Estimators[i] == "coverage"){
						sumCalculators.push_back(new Coverage());
					}else if(Estimators[i] == "geometric"){
						sumCalculators.push_back(new Geom());
					}else if(Estimators[i] == "logseries"){
						sumCalculators.push_back(new LogSD());
					}else if(Estimators[i] == "qstat"){
						sumCalculators.push_back(new QStat());
					}else if(Estimators[i] == "bergerparker"){
						sumCalculators.push_back(new BergerParker());
					}else if(Estimators[i] == "bstick"){
						sumCalculators.push_back(new BStick());
					}else if(Estimators[i] == "ace"){
						if(abund < 5)
							abund = 10;
						sumCalculators.push_back(new Ace(abund));
					}else if(Estimators[i] == "jack"){
						sumCalculators.push_back(new Jackknife());
					}else if(Estimators[i] == "shannon"){
						sumCalculators.push_back(new Shannon());
					}else if(Estimators[i] == "shannoneven"){
						sumCalculators.push_back(new ShannonEven());
					}else if(Estimators[i] == "npshannon"){
						sumCalculators.push_back(new NPShannon());
					}else if(Estimators[i] == "heip"){
						sumCalculators.push_back(new Heip());
					}else if(Estimators[i] == "smithwilson"){
						sumCalculators.push_back(new SmithWilson());
					}else if(Estimators[i] == "simpson"){
						sumCalculators.push_back(new Simpson());
					}else if(Estimators[i] == "simpsoneven"){
						sumCalculators.push_back(new SimpsonEven());
					}else if(Estimators[i] == "invsimpson"){
						sumCalculators.push_back(new InvSimpson());
					}else if(Estimators[i] == "bootstrap"){
						sumCalculators.push_back(new Bootstrap());
					}else if (Estimators[i] == "nseqs") { 
						sumCalculators.push_back(new NSeqs());
					}else if (Estimators[i] == "goodscoverage") { 
						sumCalculators.push_back(new GoodsCoverage());
					}else if (Estimators[i] == "efron") { 
						sumCalculators.push_back(new Efron(size));
					}else if (Estimators[i] == "boneh") { 
						sumCalculators.push_back(new Boneh(size));
					}else if (Estimators[i] == "solow") { 
						sumCalculators.push_back(new Solow(size));
					}else if (Estimators[i] == "shen") { 
						sumCalculators.push_back(new Shen(size, abund));
					}
				}
			}
			
			//if the users entered no valid calculators don't execute command
			if (sumCalculators.size() == 0) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
			
			ofstream outputFileHandle;
			m->openOutputFile(fileNameRoot, outputFileHandle);
			outputFileHandle << "label";
		
			input = new InputData(inputFileNames[p], format);
			sabund = input->getSAbundVector();
			string lastLabel = sabund->getLabel();
		
			for(int i=0;i<sumCalculators.size();i++){
				if(sumCalculators[i]->getCols() == 1){
					outputFileHandle << '\t' << sumCalculators[i]->getName();
					numCols++;
				}
				else{
					outputFileHandle << '\t' << sumCalculators[i]->getName() << "\t" << sumCalculators[i]->getName() << "_lci\t" << sumCalculators[i]->getName() << "_hci";
					numCols += 3;
				}
			}
			outputFileHandle << endl;
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (m->control_pressed) {  outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }  delete sabund;  delete input;  return 0;  }
			
			while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) { outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }  delete sabund;  delete input;  return 0;  }
				
				if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
					
					m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
					
					outputFileHandle << sabund->getLabel();
					for(int i=0;i<sumCalculators.size();i++){
						vector<double> data = sumCalculators[i]->getValues(sabund);
						
						if (m->control_pressed) { outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }  delete sabund;  delete input;  return 0;  }

						outputFileHandle << '\t';
						sumCalculators[i]->print(outputFileHandle);
					}
					outputFileHandle << endl;
					numLines++;
				}
				
				if ((m->anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = sabund->getLabel();
					
					delete sabund;
					sabund = input->getSAbundVector(lastLabel);
					
					m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
					
					outputFileHandle << sabund->getLabel();
					for(int i=0;i<sumCalculators.size();i++){
						vector<double> data = sumCalculators[i]->getValues(sabund);
						
						if (m->control_pressed) {  outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete sabund;  delete input;  return 0;  }
						
						outputFileHandle << '\t';
						sumCalculators[i]->print(outputFileHandle);
					}
					outputFileHandle << endl;
					numLines++;
					
					//restore real lastlabel to save below
					sabund->setLabel(saveLabel);
				}		
				
				lastLabel = sabund->getLabel();			
				
				delete sabund;
				sabund = input->getSAbundVector();
			}
			
			if (m->control_pressed) {  outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }   delete input;  return 0;  }

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
				if (sabund != NULL) {	delete sabund;	}
				sabund = input->getSAbundVector(lastLabel);
				
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				outputFileHandle << sabund->getLabel();
				for(int i=0;i<sumCalculators.size();i++){
					vector<double> data = sumCalculators[i]->getValues(sabund);
					
					if (m->control_pressed) {  outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }  delete sabund;  delete input; return 0;  }

					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				outputFileHandle << endl;
				numLines++;
				delete sabund;
			}
			
			outputFileHandle.close();
			
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }   delete input;  return 0;  }

			
			delete input;  
			for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
		}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }  return 0;  }
		
		//create summary file containing all the groups data for each label - this function just combines the info from the files already created.
		if ((sharedfile != "") && (groupMode)) {   outputNames.push_back(createGroupSummaryFile(numLines, numCols, outputNames));  }
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }  return 0;  }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SummaryCommand::parseSharedFile(string filename) {
	try {
		vector<string> filenames;
		
		map<string, ofstream*> filehandles;
		map<string, ofstream*>::iterator it3;
		
		input = new InputData(filename, "sharedfile");
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		
		string sharedFileRoot = m->getRootName(filename);
		
		//clears file before we start to write to it below
		for (int i=0; i<lookup.size(); i++) {
			remove((sharedFileRoot + lookup[i]->getGroup() + ".rabund").c_str());
			filenames.push_back((sharedFileRoot + lookup[i]->getGroup() + ".rabund"));
		}
		
		ofstream* temp;
		for (int i=0; i<lookup.size(); i++) {
			temp = new ofstream;
			filehandles[lookup[i]->getGroup()] = temp;
			groups.push_back(lookup[i]->getGroup());
		}

		while(lookup[0] != NULL) {
		
			for (int i = 0; i < lookup.size(); i++) {
				RAbundVector rav = lookup[i]->getRAbundVector();
				m->openOutputFileAppend(sharedFileRoot + lookup[i]->getGroup() + ".rabund", *(filehandles[lookup[i]->getGroup()]));
				rav.print(*(filehandles[lookup[i]->getGroup()]));
				(*(filehandles[lookup[i]->getGroup()])).close();
			}
		
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		//free memory
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}
		
		delete input;

		return filenames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "parseSharedFile");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryCommand::createGroupSummaryFile(int numLines, int numCols, vector<string>& outputNames) {
	try {
		
		ofstream out;
		string combineFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + "groups.summary";
		
		//open combined file
		m->openOutputFile(combineFileName, out);
		
		//open each groups summary file
		string newLabel = "";
		map<string, vector<string> > files;
		for (int i=0; i<outputNames.size(); i++) {
			vector<string> thisFilesLines;
			
			ifstream temp;
			m->openInputFile(outputNames[i], temp);
			
			//read through first line - labels
			string tempLabel;
			if (i == 0) { //we want to save the labels to output below
				for (int j = 0; j < numCols+1; j++) {  
					temp >> tempLabel; 
					
					if (j == 1) {  newLabel += "group\t" + tempLabel + '\t';
					}else{  newLabel += tempLabel + '\t';	}
				}
			}else{  for (int j = 0; j < numCols+1; j++) {  temp >> tempLabel;  }  }
			
			m->gobble(temp);
			
			//for each label
			for (int k = 0; k < numLines; k++) {
				
				string thisLine = "";
				string tempLabel;
					
				for (int j = 0; j < numCols+1; j++) {  
					temp >> tempLabel; 
						
					//save for later
					if (j == 1) { thisLine += groups[i] + "\t" + tempLabel + "\t";	}
					else{  thisLine += tempLabel + "\t";	}
				}
					
				thisLine += "\n";
				
				thisFilesLines.push_back(thisLine);
					
				m->gobble(temp);
			}
				
			files[outputNames[i]] = thisFilesLines;
			
			temp.close();
			remove(outputNames[i].c_str());
		}
		
		//output label line to new file
		out << newLabel << endl;
		
		//for each label
		for (int k = 0; k < numLines; k++) {
		
			//grab summary data for each group
			for (int i=0; i<outputNames.size(); i++) {
				out << files[outputNames[i]][k];
			}
		}	
		
		outputNames.clear();
		
		out.close();
		
		//return combine file name
		return combineFileName;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "createGroupSummaryFile");
		exit(1);
	}
}
//**********************************************************************************************************************
