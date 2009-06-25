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
CollectCommand::CollectCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		lines.clear();
		labels.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"freq","line","label","calc","abund","size"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { mothurOut("You must read a list, sabund or rabund before you can use the collect.single command."); mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter.validFile(parameters, "line", false);				
			if (line == "not found") { line = "";  }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//make sure user did not use both the line and label parameters
			if ((line != "") && (label != "")) { mothurOut("You cannot use both the line and label parameters at the same time. "); mothurOutEndLine(); abort = true; }
			//if the user has not specified any line or labels use the ones from read.otu
			else if((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
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
				
			if (abort == false) {
				string fileNameRoot = getRootName(globaldata->inputFileName);
				int i;
				validCalculator = new ValidCalculators();
				
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("single", Estimators[i]) == true) { 
						if (Estimators[i] == "sobs") { 
							cDisplays.push_back(new CollectDisplay(new Sobs(), new OneColumnFile(fileNameRoot+"sobs")));
						}else if (Estimators[i] == "chao") { 
							cDisplays.push_back(new CollectDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"chao")));
						}else if (Estimators[i] == "nseqs") { 
							cDisplays.push_back(new CollectDisplay(new NSeqs(), new OneColumnFile(fileNameRoot+"nseqs")));
						}else if (Estimators[i] == "coverage") { 
							cDisplays.push_back(new CollectDisplay(new Coverage(), new OneColumnFile(fileNameRoot+"coverage")));
						}else if (Estimators[i] == "ace") { 
							cDisplays.push_back(new CollectDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"ace")));
						}else if (Estimators[i] == "jack") { 
							cDisplays.push_back(new CollectDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"jack")));
						}else if (Estimators[i] == "shannon") { 
							cDisplays.push_back(new CollectDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"shannon")));
						}else if (Estimators[i] == "npshannon") { 
							cDisplays.push_back(new CollectDisplay(new NPShannon(), new OneColumnFile(fileNameRoot+"np_shannon")));
						}else if (Estimators[i] == "simpson") { 
							cDisplays.push_back(new CollectDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"simpson")));
						}else if (Estimators[i] == "bootstrap") { 
							cDisplays.push_back(new CollectDisplay(new Bootstrap(), new OneColumnFile(fileNameRoot+"bootstrap")));
						}else if (Estimators[i] == "geometric") { 
							cDisplays.push_back(new CollectDisplay(new Geom(), new OneColumnFile(fileNameRoot+"geometric")));
						}else if (Estimators[i] == "qstat") { 
							cDisplays.push_back(new CollectDisplay(new QStat(), new OneColumnFile(fileNameRoot+"qstat")));
						}else if (Estimators[i] == "logseries") { 
							cDisplays.push_back(new CollectDisplay(new LogSD(), new OneColumnFile(fileNameRoot+"logseries")));
						}else if (Estimators[i] == "bergerparker") { 
							cDisplays.push_back(new CollectDisplay(new BergerParker(), new OneColumnFile(fileNameRoot+"bergerparker")));
						}else if (Estimators[i] == "bstick") { 
							cDisplays.push_back(new CollectDisplay(new BStick(), new ThreeColumnFile(fileNameRoot+"bstick")));
						}else if (Estimators[i] == "goodscoverage") { 
							cDisplays.push_back(new CollectDisplay(new GoodsCoverage(), new OneColumnFile(fileNameRoot+"goodscoverage")));
						}else if (Estimators[i] == "efron") {
							cDisplays.push_back(new CollectDisplay(new Efron(size), new OneColumnFile(fileNameRoot+"efron")));
						}else if (Estimators[i] == "boneh") {
							cDisplays.push_back(new CollectDisplay(new Boneh(size), new OneColumnFile(fileNameRoot+"boneh")));
						}else if (Estimators[i] == "solow") {
							cDisplays.push_back(new CollectDisplay(new Solow(size), new OneColumnFile(fileNameRoot+"solow")));
						}else if (Estimators[i] == "shen") {
							cDisplays.push_back(new CollectDisplay(new Shen(size, abund), new OneColumnFile(fileNameRoot+"shen")));
						}
					}
				}
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "CollectCommand", "CollectCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

void CollectCommand::help(){
	try {
		mothurOut("The collect.single command can only be executed after a successful read.otu command. WITH ONE EXECEPTION. \n");
		mothurOut("The collect.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.\n");
		mothurOut("The collect.single command parameters are label, line, freq, calc and abund.  No parameters are required, but you may not use \n");
		mothurOut("both the line and label parameters at the same time. The collect.single command should be in the following format: \n");
		mothurOut("collect.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, calc=yourEstimators).\n");
		mothurOut("Example collect(label=unique-.01-.03, line=0-5-10, iters=10000, freq=10, calc=sobs-chao-ace-jack).\n");
		mothurOut("The default values for freq is 100, and calc are sobs-chao-ace-jack-shannon-npshannon-simpson.\n");
		validCalculator->printCalc("single", cout);
		mothurOut("The label and line parameters are used to analyze specific lines in your input.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "CollectCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

CollectCommand::~CollectCommand(){
	if (abort == false) {
		//delete order;
		delete input;  globaldata->ginput = NULL;
		delete read;
		delete validCalculator;
		globaldata->gorder = NULL;
	}
}

//**********************************************************************************************************************

int CollectCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
	
		int count = 1;
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }

		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		order = globaldata->gorder;
		string lastLabel = order->getLabel();
		input = globaldata->ginput;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;
		
		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
			if(allLines == 1 || lines.count(count) == 1 || labels.count(order->getLabel()) == 1){
				
				cCurve = new Collect(order, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
			
				mothurOut(order->getLabel() + "\t" + toString(count)); mothurOutEndLine();
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
				userLines.erase(count);
			
			//you have a label the user want that is smaller than this line and the last line has not already been processed 
			}
			
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				
				delete order;
				order = (input->getOrderVector(lastLabel));
				
				cCurve = new Collect(order, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
			
				mothurOut(order->getLabel() + "\t" + toString(count)); mothurOutEndLine();
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			}
			
			lastLabel = order->getLabel();	
			
			delete order;		
			order = (input->getOrderVector());
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
				needToRun = true;
			}else {
				mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			delete order;
			order = (input->getOrderVector(lastLabel));
			
			mothurOut(order->getLabel() + "\t" + toString(count)); mothurOutEndLine();
			
			cCurve = new Collect(order, cDisplays);
			cCurve->getCurve(freq);
			delete cCurve;
			delete order;
		}
		
		
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "CollectCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
