/*
 *  rarefactcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactcommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "jackknife.h"
#include "coverage.h"

//**********************************************************************************************************************


RareFactCommand::RareFactCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"iters","freq","label","calc","abund"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { mothurOut("You must read a list, sabund or rabund before you can use the rarefaction.single command."); mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if(label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs";  }
			else { 
				 if (calc == "default")  {  calc = "sobs";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "abund", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, nIters); 
			
			if (abort == false) {
			
				string fileNameRoot = getRootName(globaldata->inputFileName);
				int i;
				validCalculator = new ValidCalculators();
				
				
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("rarefaction", Estimators[i]) == true) { 
						if (Estimators[i] == "sobs") { 
							rDisplays.push_back(new RareDisplay(new Sobs(), new ThreeColumnFile(fileNameRoot+"rarefaction")));
						}else if (Estimators[i] == "chao") { 
							rDisplays.push_back(new RareDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"r_chao")));
						}else if (Estimators[i] == "ace") { 
							if(abund < 5)
								abund = 10;
							rDisplays.push_back(new RareDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"r_ace")));
						}else if (Estimators[i] == "jack") { 
							rDisplays.push_back(new RareDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"r_jack")));
						}else if (Estimators[i] == "shannon") { 
							rDisplays.push_back(new RareDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"r_shannon")));
						}else if (Estimators[i] == "npshannon") { 
							rDisplays.push_back(new RareDisplay(new NPShannon(), new ThreeColumnFile(fileNameRoot+"r_npshannon")));
						}else if (Estimators[i] == "simpson") { 
							rDisplays.push_back(new RareDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"r_simpson")));
						}else if (Estimators[i] == "bootstrap") { 
							rDisplays.push_back(new RareDisplay(new Bootstrap(), new ThreeColumnFile(fileNameRoot+"r_bootstrap")));
						}else if (Estimators[i] == "coverage") { 
							rDisplays.push_back(new RareDisplay(new Coverage(), new ThreeColumnFile(fileNameRoot+"r_coverage")));
						}else if (Estimators[i] == "nseqs") { 
							rDisplays.push_back(new RareDisplay(new NSeqs(), new ThreeColumnFile(fileNameRoot+"r_nseqs")));
						}
					}
				}
			}
				
		}
		
	}
	catch(exception& e) {
		errorOut(e, "RareFactCommand", "RareFactCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void RareFactCommand::help(){
	try {
		mothurOut("The rarefaction.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION.\n");
		mothurOut("The rarefaction.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.\n");
		mothurOut("The rarefaction.single command parameters are label, iters, freq, calc and abund.  No parameters are required. \n");
		mothurOut("The rarefaction.single command should be in the following format: \n");
		mothurOut("rarefaction.single(label=yourLabel, iters=yourIters, freq=yourFreq, calc=yourEstimators).\n");
		mothurOut("Example rarefaction.single(label=unique-.01-.03, iters=10000, freq=10, calc=sobs-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson).\n");
		mothurOut("The default values for iters is 1000, freq is 100, and calc is rarefaction which calculates the rarefaction curve for the observed richness.\n");
		validCalculator->printCalc("rarefaction", cout);
		mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "RareFactCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

RareFactCommand::~RareFactCommand(){
	if (abort == false) {
		globaldata->gorder = NULL;
		delete input;  globaldata->ginput = NULL;
		delete read;
		delete validCalculator;
	}
}

//**********************************************************************************************************************

int RareFactCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		//if the users entered no valid calculators don't execute command
		if (rDisplays.size() == 0) { return 0; }

		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 

		order = globaldata->gorder;
		string lastLabel = order->getLabel();
		input = globaldata->ginput;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
	
		//as long as you are not at the end of the file or done wih the lines you want
		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(order->getLabel()) == 1){
			
				rCurve = new Rarefact(order, rDisplays);
				rCurve->getCurve(freq, nIters);
				delete rCurve;
			
				mothurOut(order->getLabel()); mothurOutEndLine();
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			}
			
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = order->getLabel();
				
				delete order;
				order = (input->getOrderVector(lastLabel));
				
				rCurve = new Rarefact(order, rDisplays);
				rCurve->getCurve(freq, nIters);
				delete rCurve;
			
				mothurOut(order->getLabel()); mothurOutEndLine();
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
				
				//restore real lastlabel to save below
				order->setLabel(saveLabel);
			}
			
			lastLabel = order->getLabel();		
			
			delete order;
			order = (input->getOrderVector());
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
		
		//run last label if you need to
		if (needToRun == true)  {
			if (order != NULL) {	delete order;	}
			order = (input->getOrderVector(lastLabel));
				
			rCurve = new Rarefact(order, rDisplays);
			rCurve->getCurve(freq, nIters);
			delete rCurve;
			
			mothurOut(order->getLabel()); mothurOutEndLine();
			delete order;
		}
		

		for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "RareFactCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
