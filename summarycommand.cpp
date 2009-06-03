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

SummaryCommand::SummaryCommand(){
	try {
		globaldata = GlobalData::getInstance();
		validCalculator = new ValidCalculators();
		int i;
		
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("summary", globaldata->Estimators[i]) == true) { 
				if(globaldata->Estimators[i] == "sobs"){
					sumCalculators.push_back(new Sobs());
				}else if(globaldata->Estimators[i] == "chao"){
					sumCalculators.push_back(new Chao1());
				}else if(globaldata->Estimators[i] == "coverage"){
					sumCalculators.push_back(new Coverage());
				}else if(globaldata->Estimators[i] == "geometric"){
					sumCalculators.push_back(new Geom());
				}else if(globaldata->Estimators[i] == "logseries"){
					sumCalculators.push_back(new LogSD());
				}else if(globaldata->Estimators[i] == "qstat"){
					sumCalculators.push_back(new QStat());
				}else if(globaldata->Estimators[i] == "bergerparker"){
					sumCalculators.push_back(new BergerParker());
				}else if(globaldata->Estimators[i] == "bstick"){
					sumCalculators.push_back(new BStick());
				}else if(globaldata->Estimators[i] == "ace"){
					convert(globaldata->getAbund(), abund);
					if(abund < 5)
						abund = 10;
					sumCalculators.push_back(new Ace(abund));
				}else if(globaldata->Estimators[i] == "jack"){
					sumCalculators.push_back(new Jackknife());
				}else if(globaldata->Estimators[i] == "shannon"){
					sumCalculators.push_back(new Shannon());
				}else if(globaldata->Estimators[i] == "npshannon"){
					sumCalculators.push_back(new NPShannon());
				}else if(globaldata->Estimators[i] == "simpson"){
					sumCalculators.push_back(new Simpson());
				}else if(globaldata->Estimators[i] == "bootstrap"){
					sumCalculators.push_back(new Bootstrap());
				}else if (globaldata->Estimators[i] == "nseqs") { 
					sumCalculators.push_back(new NSeqs());
				}else if (globaldata->Estimators[i] == "goodscoverage") { 
					sumCalculators.push_back(new GoodsCoverage());
				}else if (globaldata->Estimators[i] == "efron") { 
					convert(globaldata->getSize(), size);
					sumCalculators.push_back(new Efron(size));
				}else if (globaldata->Estimators[i] == "boneh") { 
					convert(globaldata->getSize(), size);
					sumCalculators.push_back(new Boneh(size));
				}else if (globaldata->Estimators[i] == "solow") { 
					convert(globaldata->getSize(), size);
					sumCalculators.push_back(new Solow(size));
				}else if (globaldata->Estimators[i] == "shen") { 
					convert(globaldata->getSize(), size);
					sumCalculators.push_back(new Shen(size));
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");

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

SummaryCommand::~SummaryCommand(){
	delete sabund;
	delete input;
	delete read;
}

//**********************************************************************************************************************

int SummaryCommand::execute(){
	try {
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
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;
		
		while((sabund != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(sabund->getLabel()) == 1){			
	
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
