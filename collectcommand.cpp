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
CollectCommand::CollectCommand(){
	try {
		globaldata = GlobalData::getInstance();
		string fileNameRoot;
		fileNameRoot = getRootName(globaldata->inputFileName);
		convert(globaldata->getFreq(), freq);
		int i;
		validCalculator = new ValidCalculators();
		
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("single", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sobs") { 
					cDisplays.push_back(new CollectDisplay(new Sobs(), new OneColumnFile(fileNameRoot+"sobs")));
				}else if (globaldata->Estimators[i] == "chao") { 
					cDisplays.push_back(new CollectDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"chao")));
				}else if (globaldata->Estimators[i] == "nseqs") { 
					cDisplays.push_back(new CollectDisplay(new NSeqs(), new OneColumnFile(fileNameRoot+"nseqs")));
				}else if (globaldata->Estimators[i] == "coverage") { 
					cDisplays.push_back(new CollectDisplay(new Coverage(), new OneColumnFile(fileNameRoot+"coverage")));
				}else if (globaldata->Estimators[i] == "ace") { 
					convert(globaldata->getAbund(), abund);
					cDisplays.push_back(new CollectDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"ace")));
				}else if (globaldata->Estimators[i] == "jack") { 
					cDisplays.push_back(new CollectDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"jack")));
				}else if (globaldata->Estimators[i] == "shannon") { 
					cDisplays.push_back(new CollectDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"shannon")));
				}else if (globaldata->Estimators[i] == "npshannon") { 
					cDisplays.push_back(new CollectDisplay(new NPShannon(), new OneColumnFile(fileNameRoot+"np_shannon")));
				}else if (globaldata->Estimators[i] == "simpson") { 
					cDisplays.push_back(new CollectDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"simpson")));
				}else if (globaldata->Estimators[i] == "bootstrap") { 
					cDisplays.push_back(new CollectDisplay(new Bootstrap(), new OneColumnFile(fileNameRoot+"bootstrap")));
				}else if (globaldata->Estimators[i] == "geometric") { 
					cDisplays.push_back(new CollectDisplay(new Geom(), new OneColumnFile(fileNameRoot+"geometric")));
				}else if (globaldata->Estimators[i] == "qstat") { 
					cDisplays.push_back(new CollectDisplay(new QStat(), new OneColumnFile(fileNameRoot+"qstat")));
				}else if (globaldata->Estimators[i] == "logseries") { 
					cDisplays.push_back(new CollectDisplay(new LogSD(), new OneColumnFile(fileNameRoot+"logseries")));
				}else if (globaldata->Estimators[i] == "bergerparker") { 
					cDisplays.push_back(new CollectDisplay(new BergerParker(), new OneColumnFile(fileNameRoot+"bergerparker")));
				}else if (globaldata->Estimators[i] == "bstick") { 
					cDisplays.push_back(new CollectDisplay(new BStick(), new ThreeColumnFile(fileNameRoot+"bstick")));
				}else if (globaldata->Estimators[i] == "goodscoverage") { 
					cDisplays.push_back(new CollectDisplay(new GoodsCoverage(), new OneColumnFile(fileNameRoot+"goodscoverage")));
				}else if (globaldata->Estimators[i] == "efron") {
					convert(globaldata->getSize(), size); 
					cDisplays.push_back(new CollectDisplay(new Efron(size), new OneColumnFile(fileNameRoot+"efron")));
				}else if (globaldata->Estimators[i] == "boneh") {
					convert(globaldata->getSize(), size); 
					cDisplays.push_back(new CollectDisplay(new Boneh(size), new OneColumnFile(fileNameRoot+"boneh")));
				}else if (globaldata->Estimators[i] == "solow") {
					convert(globaldata->getSize(), size); 
					cDisplays.push_back(new CollectDisplay(new Solow(size), new OneColumnFile(fileNameRoot+"solow")));
				}else if (globaldata->Estimators[i] == "shen") {
					convert(globaldata->getSize(), size); 
					cDisplays.push_back(new CollectDisplay(new Shen(size), new OneColumnFile(fileNameRoot+"shen")));
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectCommand class Function CollectCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectCommand class function CollectCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

CollectCommand::~CollectCommand(){
	delete order;
	delete input;
	delete cCurve;
	delete read;
}

//**********************************************************************************************************************

int CollectCommand::execute(){
	try {
		int count = 1;
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }

		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		order = globaldata->gorder;
		lastOrder = order;
		input = globaldata->ginput;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		
		while((order != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0))) {
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
				
				cCurve = new Collect(order, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
			
				cout << order->getLabel() << '\t' << count << endl;
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			
			//you have a label the user want that is smaller than this line and the last line has not already been processed 
			}
			
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastOrder->getLabel()) != 1)) {
				cCurve = new Collect(lastOrder, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
			
				cout << lastOrder->getLabel() << '\t' << count << endl;
				processedLabels.insert(lastOrder->getLabel());
				userLabels.erase(lastOrder->getLabel());
			}
			
			if (count != 1) { delete lastOrder; }
			lastOrder = order;			
			order = (input->getOrderVector());
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastOrder->getLabel()) != 1) {
				cout << ". I will use " << lastOrder->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastOrder->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cCurve = new Collect(lastOrder, cDisplays);
			cCurve->getCurve(freq);
			delete cCurve;
			
			cout << lastOrder->getLabel() << '\t' << count << endl;
		}
		
		delete lastOrder;
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
