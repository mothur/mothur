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
		lines.clear();
		labels.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"line","label","calc","abund","size"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			parser = new OptionParser();
			parser->parse(option, parameters);  delete parser;
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter->isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { cout << "You must read a list, sabund or rabund before you can use the summary.single command." << endl; abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter->validFile(parameters, "line", false);				
			if (line == "not found") { line = "";  }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			label = validParameter->validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//make sure user did not use both the line and label parameters
			if ((line != "") && (label != "")) { cout << "You cannot use both the line and label parameters at the same time. " << endl; abort = true; }
			//if the user has not specified any line or labels use the ones from read.otu
			else if((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
				
			calc = validParameter->validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
				 if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter->validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter->validFile(parameters, "size", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, size); 
	
			delete validParameter;
			
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
		cout << "Standard Error: " << e.what() << " has occurred in the SummaryCommand class Function SummaryCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummaryCommand class function SummaryCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

void SummaryCommand::help(){
	try {
		cout << "The summary.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The summary.single command parameters are label, line, calc, abund.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.single command should be in the following format: " << "\n";
		cout << "summary.single(label=yourLabel, line=yourLines, calc=yourEstimators)." << "\n";
		cout << "Example summary.single(label=unique-.01-.03, line=0,5,10, calc=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson)." << "\n";
		validCalculator->printCalc("summary", cout);
		cout << "The default value calc is sobs-chao-ace-jack-shannon-npshannon-simpson" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummaryCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummaryCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

SummaryCommand::~SummaryCommand(){
	delete sabund;
	delete input;
	delete read;
	delete validCalculator;
}

//**********************************************************************************************************************

int SummaryCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		int count = 1;
		
		//if the users entered no valid calculators don't execute command
		if (sumCalculators.size() == 0) { return 0; }

		outputFileName = ((getRootName(globaldata->inputFileName)) + "summary");
		openOutputFile(outputFileName, outputFileHandle);
		outputFileHandle << "label";
	
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		sabund = globaldata->sabund;
		SAbundVector* lastSAbund = sabund;
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
		set<int> userLines = lines;
		
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(allLines == 1 || lines.count(count) == 1 || labels.count(sabund->getLabel()) == 1){			
	
				cout << sabund->getLabel() << '\t' << count << endl;
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				userLines.erase(count);

				
				outputFileHandle << sabund->getLabel();
				for(int i=0;i<sumCalculators.size();i++){
					vector<double> data = sumCalculators[i]->getValues(sabund);
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				outputFileHandle << endl;
			}
			
			if ((anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastSAbund->getLabel()) != 1)) {

				cout << lastSAbund->getLabel() << '\t' << count << endl;
				processedLabels.insert(lastSAbund->getLabel());
				userLabels.erase(lastSAbund->getLabel());
				
				outputFileHandle << lastSAbund->getLabel();
				for(int i=0;i<sumCalculators.size();i++){
					vector<double> data = sumCalculators[i]->getValues(lastSAbund);
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				outputFileHandle << endl;
			}		

			if (count != 1) { delete lastSAbund; }
			lastSAbund = sabund;			

			sabund = input->getSAbundVector();
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastSAbund->getLabel()) != 1) {
				cout << ". I will use " << lastSAbund->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastSAbund->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastSAbund->getLabel() << '\t' << count << endl;
			outputFileHandle << lastSAbund->getLabel();
			for(int i=0;i<sumCalculators.size();i++){
				vector<double> data = sumCalculators[i]->getValues(lastSAbund);
				outputFileHandle << '\t';
				sumCalculators[i]->print(outputFileHandle);
			}
			outputFileHandle << endl;
		}
		
		delete lastSAbund;
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummaryCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummaryCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
