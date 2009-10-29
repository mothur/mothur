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
#include "npshannon.h"
#include "shannon.h"
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

SummaryCommand::SummaryCommand(string option){
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
			string Array[] =  {"label","calc","abund","size"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { mothurOut("You must read a list, sabund or rabund before you can use the summary.single command."); mothurOutEndLine(); abort = true; }
			
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
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
				 if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "size", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, size); 
	
			if (abort == false) {
			
				validCalculator = new ValidCalculators();
				int i;
				
				for (i=0; i<Estimators.size(); i++) {
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
						}else if(Estimators[i] == "npshannon"){
							sumCalculators.push_back(new NPShannon());
						}else if(Estimators[i] == "simpson"){
							sumCalculators.push_back(new Simpson());
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
			}
		}

				
	}
	catch(exception& e) {
		errorOut(e, "SummaryCommand", "SummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SummaryCommand::help(){
	try {
		mothurOut("The summary.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION.\n");
		mothurOut("The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.\n");
		mothurOut("The summary.single command parameters are label, calc, abund.  No parameters are required.\n");
		mothurOut("The summary.single command should be in the following format: \n");
		mothurOut("summary.single(label=yourLabel, calc=yourEstimators).\n");
		mothurOut("Example summary.single(label=unique-.01-.03, calc=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson).\n");
		validCalculator->printCalc("summary", cout);
		mothurOut("The default value calc is sobs-chao-ace-jack-shannon-npshannon-simpson\n");
		mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabels).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "SummaryCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

SummaryCommand::~SummaryCommand(){
	if (abort == false) {
		delete input;  globaldata->ginput = NULL;
		delete read;
		delete validCalculator;
		globaldata->sabund = NULL;
	}
}

//**********************************************************************************************************************

int SummaryCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		//if the users entered no valid calculators don't execute command
		if (sumCalculators.size() == 0) { return 0; }

		outputFileName = ((getRootName(globaldata->inputFileName)) + "summary");
		openOutputFile(outputFileName, outputFileHandle);
		outputFileHandle << "label";
	
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		sabund = globaldata->sabund;
		string lastLabel = sabund->getLabel();
		input = globaldata->ginput;
		
		for(int i=0;i<sumCalculators.size();i++){
			if(sumCalculators[i]->getCols() == 1){
				outputFileHandle << '\t' << sumCalculators[i]->getName();
			}
			else{
				outputFileHandle << '\t' << sumCalculators[i]->getName() << "\t" << sumCalculators[i]->getName() << "_lci\t" << sumCalculators[i]->getName() << "_hci";
			}
		}
		outputFileHandle << endl;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
			
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
	
				mothurOut(sabund->getLabel()); mothurOutEndLine();
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
								
				outputFileHandle << sabund->getLabel();
				for(int i=0;i<sumCalculators.size();i++){
					vector<double> data = sumCalculators[i]->getValues(sabund);
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				outputFileHandle << endl;
			}
			
			if ((anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = sabund->getLabel();
				
				delete sabund;
				sabund = input->getSAbundVector(lastLabel);
				
				mothurOut(sabund->getLabel()); mothurOutEndLine();
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				outputFileHandle << sabund->getLabel();
				for(int i=0;i<sumCalculators.size();i++){
					vector<double> data = sumCalculators[i]->getValues(sabund);
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				outputFileHandle << endl;
				
				//restore real lastlabel to save below
				sabund->setLabel(saveLabel);
			}		

			lastLabel = sabund->getLabel();			
			
			delete sabund;
			sabund = input->getSAbundVector();
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
			if (sabund != NULL) {	delete sabund;	}
			sabund = input->getSAbundVector(lastLabel);
			
			mothurOut(sabund->getLabel()); mothurOutEndLine();
			outputFileHandle << sabund->getLabel();
			for(int i=0;i<sumCalculators.size();i++){
				vector<double> data = sumCalculators[i]->getValues(sabund);
				outputFileHandle << '\t';
				sumCalculators[i]->print(outputFileHandle);
			}
			outputFileHandle << endl;
			delete sabund;
		}
		
		outputFileHandle.close();
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "SummaryCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
