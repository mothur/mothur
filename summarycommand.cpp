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
vector<string> SummaryCommand::getValidParameters(){	
	try {
		string Array[] =  {"label","calc","abund","size","outputdir","groupmode","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SummaryCommand::SummaryCommand(){	
	try {
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "SummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SummaryCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SummaryCommand::getRequiredFiles(){	
	try {
		string AlignArray[] =  {"shared","list","rabund","sabund","or"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

SummaryCommand::SummaryCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); delete validCalculator; abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"label","calc","abund","size","outputdir","groupmode","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getSharedFile() == "") && (globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { m->mothurOut("You must read a list, sabund, rabund or shared file before you can use the summary.single command."); m->mothurOutEndLine(); abort = true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if(label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
				 if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			m->splitAtDash(calc, Estimators);

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

void SummaryCommand::help(){
	try {
		m->mothurOut("The summary.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION.\n");
		m->mothurOut("The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.\n");
		m->mothurOut("The summary.single command parameters are label, calc, abund and groupmode.  No parameters are required.\n");
		m->mothurOut("The summary.single command should be in the following format: \n");
		m->mothurOut("summary.single(label=yourLabel, calc=yourEstimators).\n");
		m->mothurOut("Example summary.single(label=unique-.01-.03, calc=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson).\n");
		validCalculator->printCalc("summary", cout);
		m->mothurOut("The default value calc is sobs-chao-ace-jack-shannon-npshannon-simpson\n");
		m->mothurOut("If you are running summary.single with a shared file and would like your summary results collated in one file, set groupmode=t. (Default=true).\n");
		m->mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabels).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

SummaryCommand::~SummaryCommand(){}

//**********************************************************************************************************************

int SummaryCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		string hadShared = "";
		if ((globaldata->getFormat() != "sharedfile")) { inputFileNames.push_back(globaldata->inputFileName);  }
		else { hadShared = globaldata->getSharedFile(); inputFileNames = parseSharedFile(globaldata->getSharedFile());  globaldata->setFormat("rabund");  }
		
		if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } return 0; }
		
		int numLines = 0;
		int numCols = 0;
		
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			numLines = 0;
			numCols = 0;
			
			string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(inputFileNames[p])) + "summary";
			globaldata->inputFileName = inputFileNames[p];
			outputNames.push_back(fileNameRoot); outputTypes["summary"].push_back(fileNameRoot);
			
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
			
			sumCalculators.clear();
			
			validCalculator = new ValidCalculators();
			
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("summary", Estimators[i]) == true) { 
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
			if (sumCalculators.size() == 0) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } return 0; }
			
			ofstream outputFileHandle;
			m->openOutputFile(fileNameRoot, outputFileHandle);
			outputFileHandle << "label";
			
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			sabund = globaldata->sabund;
			string lastLabel = sabund->getLabel();
			input = globaldata->ginput;
			
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
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  }  outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read; delete sabund; globaldata->sabund = NULL; delete input; globaldata->ginput = NULL; return 0;  }
			
			while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read; delete sabund; globaldata->sabund = NULL; delete input; globaldata->ginput = NULL; return 0;  }
				
				if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
					
					m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
					
					outputFileHandle << sabund->getLabel();
					for(int i=0;i<sumCalculators.size();i++){
						vector<double> data = sumCalculators[i]->getValues(sabund);
						
						if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read; delete sabund; globaldata->sabund = NULL; delete input; globaldata->ginput = NULL; return 0;  }

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
						
						if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read; delete sabund; globaldata->sabund = NULL; delete input; globaldata->ginput = NULL; return 0;  }
						
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
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read;  delete input; globaldata->ginput = NULL; return 0;  }

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
					
					if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } outputFileHandle.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read; delete sabund; globaldata->sabund = NULL; delete input; globaldata->ginput = NULL; return 0;  }

					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				outputFileHandle << endl;
				numLines++;
				delete sabund;
			}
			
			outputFileHandle.close();
			
			if (m->control_pressed) { if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  } for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } delete validCalculator; delete read;  delete input; globaldata->ginput = NULL; return 0;  }

			
			delete input;  globaldata->ginput = NULL;
			delete read;
			delete validCalculator;
			globaldata->sabund = NULL;
			for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
		}
		
		if (hadShared != "") {  globaldata->setSharedFile(hadShared); globaldata->setFormat("sharedfile");  }
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }  return 0;  }
		
		//create summary file containing all the groups data for each label - this function just combines the info from the files already created.
		if ((hadShared != "") && (groupMode)) {   outputNames.push_back(createGroupSummaryFile(numLines, numCols, outputNames));  }
		
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
		
				
		//read first line
		read = new ReadOTUFile(filename);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
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
		delete read;
		delete input;
		globaldata->ginput = NULL;

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
		string combineFileName = outputDir + m->getRootName(m->getSimpleName(globaldata->inputFileName)) + "groups.summary";
		
		//open combined file
		m->openOutputFile(combineFileName, out);
		
		//open each groups summary file
		string newLabel = "";
		ifstream* temp;
		map<string, ifstream*> filehandles;
		for (int i=0; i<outputNames.size(); i++) {
			temp = new ifstream;
			filehandles[outputNames[i]] = temp;
			m->openInputFile(outputNames[i], *(temp));
			
			//read through first line - labels
			string tempLabel;
			if (i == 0) { //we want to save the labels to output below
				for (int j = 0; j < numCols+1; j++) {  
					*(temp) >> tempLabel; 
					
					if (j == 1) {  newLabel += "group\t" + tempLabel + '\t';
					}else{  newLabel += tempLabel + '\t';	}
				}
			}else{  for (int j = 0; j < numCols+1; j++) {  *(temp) >> tempLabel;  }  }
			
			m->gobble(*(temp));
		}
		
		//output label line to new file
		out << newLabel << endl;
		
		//for each label
		for (int i = 0; i < numLines; i++) {
		
			//grab summary data for each group
			for (int i=0; i<outputNames.size(); i++) {
				string tempLabel;
				
				for (int j = 0; j < numCols+1; j++) {  
					*(filehandles[outputNames[i]]) >> tempLabel; 
					
					//print to combined file
					if (j == 1) { out << groups[i] << '\t' << tempLabel << '\t';	}
					else{  out << tempLabel << '\t';	}
				}
				
				out << endl;
				m->gobble(*(filehandles[outputNames[i]]));
			}
		}	
		
		//close each groups summary file
		for (int i=0; i<outputNames.size(); i++) {  (*(filehandles[outputNames[i]])).close();  remove(outputNames[i].c_str());  }
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
