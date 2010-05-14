/*
 *  collectcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectcommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "jackknife.h"
#include "geom.h"
#include "qstat.h"
#include "logsd.h"
#include "bergerparker.h"
#include "bstick.h"
#include "goodscoverage.h"
#include "efron.h"
#include "boneh.h"
#include "solow.h"
#include "shen.h"
#include "coverage.h"


//**********************************************************************************************************************
CollectCommand::CollectCommand(string option)  {
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
			string Array[] =  {"freq","label","calc","abund","size","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getSharedFile() == "") && (globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { m->mothurOut("You must read a list, sabund, rabund or shared file before you can use the collect.single command."); m->mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if (label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
				 if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "size", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, size); 
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "CollectCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

void CollectCommand::help(){
	try {
		m->mothurOut("The collect.single command can only be executed after a successful read.otu command. WITH ONE EXECEPTION. \n");
		m->mothurOut("The collect.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.\n");
		m->mothurOut("The collect.single command parameters are label, freq, calc and abund.  No parameters are required. \n");
		m->mothurOut("The collect.single command should be in the following format: \n");
		m->mothurOut("The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n");
		m->mothurOut("collect.single(label=yourLabel, iters=yourIters, freq=yourFreq, calc=yourEstimators).\n");
		m->mothurOut("Example collect(label=unique-.01-.03, iters=10000, freq=10, calc=sobs-chao-ace-jack).\n");
		m->mothurOut("The default values for freq is 100, and calc are sobs-chao-ace-jack-shannon-npshannon-simpson.\n");
		validCalculator->printCalc("single", cout);
		m->mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

CollectCommand::~CollectCommand(){}

//**********************************************************************************************************************

int CollectCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		vector<string> outputNames;
		
		if ((globaldata->getFormat() != "sharedfile")) { inputFileNames.push_back(globaldata->inputFileName);  }
		else {  inputFileNames = parseSharedFile(globaldata->getSharedFile());  globaldata->setFormat("rabund");  }
	
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}  globaldata->Groups.clear(); return 0; }
			
			if (outputDir == "") { outputDir += hasPath(inputFileNames[p]); }
			string fileNameRoot = outputDir + getRootName(getSimpleName(inputFileNames[p]));
			globaldata->inputFileName = inputFileNames[p];
		
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
		
			validCalculator = new ValidCalculators();
			
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("single", Estimators[i]) == true) { 
					if (Estimators[i] == "sobs") { 
						cDisplays.push_back(new CollectDisplay(new Sobs(), new OneColumnFile(fileNameRoot+"sobs")));
						outputNames.push_back(fileNameRoot+"sobs");
					}else if (Estimators[i] == "chao") { 
						cDisplays.push_back(new CollectDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"chao")));
						outputNames.push_back(fileNameRoot+"chao");
					}else if (Estimators[i] == "nseqs") { 
						cDisplays.push_back(new CollectDisplay(new NSeqs(), new OneColumnFile(fileNameRoot+"nseqs")));
						outputNames.push_back(fileNameRoot+"nseqs");
					}else if (Estimators[i] == "coverage") { 
						cDisplays.push_back(new CollectDisplay(new Coverage(), new OneColumnFile(fileNameRoot+"coverage")));
						outputNames.push_back(fileNameRoot+"coverage");
					}else if (Estimators[i] == "ace") { 
						cDisplays.push_back(new CollectDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"ace")));
						outputNames.push_back(fileNameRoot+"ace");
					}else if (Estimators[i] == "jack") { 
						cDisplays.push_back(new CollectDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"jack")));
						outputNames.push_back(fileNameRoot+"jack");
					}else if (Estimators[i] == "shannon") { 
						cDisplays.push_back(new CollectDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"shannon")));
						outputNames.push_back(fileNameRoot+"shannon");
					}else if (Estimators[i] == "npshannon") { 
						cDisplays.push_back(new CollectDisplay(new NPShannon(), new OneColumnFile(fileNameRoot+"np_shannon")));
						outputNames.push_back(fileNameRoot+"np_shannon");
					}else if (Estimators[i] == "simpson") { 
						cDisplays.push_back(new CollectDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"simpson")));
						outputNames.push_back(fileNameRoot+"simpson");
					}else if (Estimators[i] == "bootstrap") { 
						cDisplays.push_back(new CollectDisplay(new Bootstrap(), new OneColumnFile(fileNameRoot+"bootstrap")));
						outputNames.push_back(fileNameRoot+"bootstrap");
					}else if (Estimators[i] == "geometric") { 
						cDisplays.push_back(new CollectDisplay(new Geom(), new OneColumnFile(fileNameRoot+"geometric")));
						outputNames.push_back(fileNameRoot+"geometric");
					}else if (Estimators[i] == "qstat") { 
						cDisplays.push_back(new CollectDisplay(new QStat(), new OneColumnFile(fileNameRoot+"qstat")));
						outputNames.push_back(fileNameRoot+"qstat");
					}else if (Estimators[i] == "logseries") { 
						cDisplays.push_back(new CollectDisplay(new LogSD(), new OneColumnFile(fileNameRoot+"logseries")));
						outputNames.push_back(fileNameRoot+"logseries");
					}else if (Estimators[i] == "bergerparker") { 
						cDisplays.push_back(new CollectDisplay(new BergerParker(), new OneColumnFile(fileNameRoot+"bergerparker")));
						outputNames.push_back(fileNameRoot+"bergerparker");
					}else if (Estimators[i] == "bstick") { 
						cDisplays.push_back(new CollectDisplay(new BStick(), new ThreeColumnFile(fileNameRoot+"bstick")));
						outputNames.push_back(fileNameRoot+"bstick");
					}else if (Estimators[i] == "goodscoverage") { 
						cDisplays.push_back(new CollectDisplay(new GoodsCoverage(), new OneColumnFile(fileNameRoot+"goodscoverage")));
						outputNames.push_back(fileNameRoot+"goodscoverage");
					}else if (Estimators[i] == "efron") {
						cDisplays.push_back(new CollectDisplay(new Efron(size), new OneColumnFile(fileNameRoot+"efron")));
						outputNames.push_back(fileNameRoot+"efron");
					}else if (Estimators[i] == "boneh") {
						cDisplays.push_back(new CollectDisplay(new Boneh(size), new OneColumnFile(fileNameRoot+"boneh")));
						outputNames.push_back(fileNameRoot+"boneh");
					}else if (Estimators[i] == "solow") {
						cDisplays.push_back(new CollectDisplay(new Solow(size), new OneColumnFile(fileNameRoot+"solow")));
						outputNames.push_back(fileNameRoot+"solow");
					}else if (Estimators[i] == "shen") {
						cDisplays.push_back(new CollectDisplay(new Shen(size, abund), new OneColumnFile(fileNameRoot+"shen")));
						outputNames.push_back(fileNameRoot+"shen");
					}
				}
			}
		
			//if the users entered no valid calculators don't execute command
			if (cDisplays.size() == 0) { return 0; }
			
			read = new ReadOTUFile(inputFileNames[p]);	
			read->read(&*globaldata); 
		
			order = globaldata->gorder;
			string lastLabel = order->getLabel();
			input = globaldata->ginput;
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (m->control_pressed) {  
				for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}
				delete input;  globaldata->ginput = NULL;
				delete read;
				delete order; globaldata->gorder = NULL;
				delete validCalculator;
				globaldata->Groups.clear();
				return 0;
			}


			while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->control_pressed) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}
					delete input;  globaldata->ginput = NULL;
					delete read;
					delete order; globaldata->gorder = NULL;
					delete validCalculator;
					globaldata->Groups.clear();
					return 0;
				}

				
				if(allLines == 1 || labels.count(order->getLabel()) == 1){
				
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					cCurve = new Collect(order, cDisplays);
					cCurve->getCurve(freq);
					delete cCurve;
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					
				}
				//you have a label the user want that is smaller than this label and the last label has not already been processed 
				if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = order->getLabel();
					
					delete order;
					order = (input->getOrderVector(lastLabel));
					
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					cCurve = new Collect(order, cDisplays);
					cCurve->getCurve(freq);
					delete cCurve;
					
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					//restore real lastlabel to save below
					order->setLabel(saveLabel);
				}
				
				lastLabel = order->getLabel();	
				
				delete order;		
				order = (input->getOrderVector());
			}
			
			
			if (m->control_pressed) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}
					delete input;  globaldata->ginput = NULL;
					delete read;
					delete validCalculator;
					globaldata->Groups.clear();
					return 0;
			}
				
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
				if (order != NULL) {	delete order;	}
				order = (input->getOrderVector(lastLabel));
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				
				cCurve = new Collect(order, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
				
				if (m->control_pressed) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}
					delete input;  globaldata->ginput = NULL;
					delete read;
					delete order; globaldata->gorder = NULL;
					delete validCalculator;
					globaldata->Groups.clear();
					return 0;
				}
				delete order;
			}
			
			for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
			cDisplays.clear();
			delete input;  globaldata->ginput = NULL;
			delete read;
			globaldata->gorder = NULL;
			delete validCalculator;
		}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} return 0; }

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
vector<string> CollectCommand::parseSharedFile(string filename) {
	try {
		vector<string> filenames;
		
		map<string, ofstream*> filehandles;
		map<string, ofstream*>::iterator it3;
		
				
		//read first line
		read = new ReadOTUFile(filename);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		
		string sharedFileRoot = getRootName(filename);
		
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
				openOutputFileAppend(sharedFileRoot + lookup[i]->getGroup() + ".rabund", *(filehandles[lookup[i]->getGroup()]));
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
		m->errorOut(e, "CollectCommand", "parseSharedFile");
		exit(1);
	}
}
//**********************************************************************************************************************

