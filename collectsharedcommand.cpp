/*
 *  collectsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectsharedcommand.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharednseqs.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "sharedjackknife.h"
#include "sharedwhittaker.h"



//**********************************************************************************************************************

CollectSharedCommand::CollectSharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		string fileNameRoot;
		fileNameRoot = getRootName(globaldata->inputFileName);
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		util = new SharedUtil();
		
		int i;
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("shared", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sharedchao") { 
					cDisplays.push_back(new CollectDisplay(new SharedChao1(), new SharedOneColumnFile(fileNameRoot+"shared.chao")));
				}else if (globaldata->Estimators[i] == "sharedsobs") { 
					cDisplays.push_back(new CollectDisplay(new SharedSobsCS(), new SharedOneColumnFile(fileNameRoot+"shared.sobs")));
				}else if (globaldata->Estimators[i] == "sharedace") { 
					cDisplays.push_back(new CollectDisplay(new SharedAce(), new SharedOneColumnFile(fileNameRoot+"shared.ace")));
				}else if (globaldata->Estimators[i] == "jabund") { 	
					cDisplays.push_back(new CollectDisplay(new JAbund(), new SharedOneColumnFile(fileNameRoot+"jabund")));
				}else if (globaldata->Estimators[i] == "sorabund") { 
					cDisplays.push_back(new CollectDisplay(new SorAbund(), new SharedOneColumnFile(fileNameRoot+"sorabund")));
				}else if (globaldata->Estimators[i] == "jclass") { 
					cDisplays.push_back(new CollectDisplay(new Jclass(), new SharedOneColumnFile(fileNameRoot+"jclass")));
				}else if (globaldata->Estimators[i] == "sorclass") { 
					cDisplays.push_back(new CollectDisplay(new SorClass(), new SharedOneColumnFile(fileNameRoot+"sorclass")));
				}else if (globaldata->Estimators[i] == "jest") { 
					cDisplays.push_back(new CollectDisplay(new Jest(), new SharedOneColumnFile(fileNameRoot+"jest")));
				}else if (globaldata->Estimators[i] == "sorest") { 
					cDisplays.push_back(new CollectDisplay(new SorEst(), new SharedOneColumnFile(fileNameRoot+"sorest")));
				}else if (globaldata->Estimators[i] == "thetayc") { 
					cDisplays.push_back(new CollectDisplay(new ThetaYC(), new SharedOneColumnFile(fileNameRoot+"thetayc")));
				}else if (globaldata->Estimators[i] == "thetan") { 
					cDisplays.push_back(new CollectDisplay(new ThetaN(), new SharedOneColumnFile(fileNameRoot+"thetan")));
				}else if (globaldata->Estimators[i] == "kstest") { 
					cDisplays.push_back(new CollectDisplay(new KSTest(), new SharedOneColumnFile(fileNameRoot+"kstest")));
				}else if (globaldata->Estimators[i] == "whittaker") { 
					cDisplays.push_back(new CollectDisplay(new Whittaker(), new SharedOneColumnFile(fileNameRoot+"whittaker")));
				}else if (globaldata->Estimators[i] == "sharednseqs") { 
					cDisplays.push_back(new CollectDisplay(new SharedNSeqs(), new SharedOneColumnFile(fileNameRoot+"shared.nseqs")));

				}else if (globaldata->Estimators[i] == "ochiai") { 
					cDisplays.push_back(new CollectDisplay(new Ochiai(), new SharedOneColumnFile(fileNameRoot+"ochiai")));
				}else if (globaldata->Estimators[i] == "anderberg") { 
					cDisplays.push_back(new CollectDisplay(new Anderberg(), new SharedOneColumnFile(fileNameRoot+"anderberg")));
				}else if (globaldata->Estimators[i] == "skulczynski") { 
					cDisplays.push_back(new CollectDisplay(new Kulczynski(), new SharedOneColumnFile(fileNameRoot+"kulczynski")));
				}else if (globaldata->Estimators[i] == "kulczynskicody") { 
					cDisplays.push_back(new CollectDisplay(new KulczynskiCody(), new SharedOneColumnFile(fileNameRoot+"kulczynskicody")));
				}else if (globaldata->Estimators[i] == "lennon") { 
					cDisplays.push_back(new CollectDisplay(new Lennon(), new SharedOneColumnFile(fileNameRoot+"lennon")));
				}else if (globaldata->Estimators[i] == "morisitahorn") { 
					cDisplays.push_back(new CollectDisplay(new MorHorn(), new SharedOneColumnFile(fileNameRoot+"morisitahorn")));
				}else if (globaldata->Estimators[i] == "braycurtis") { 
					cDisplays.push_back(new CollectDisplay(new BrayCurtis(), new SharedOneColumnFile(fileNameRoot+"braycurtis")));
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectSharedCommand class Function CollectSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectSharedCommand class function CollectSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

CollectSharedCommand::~CollectSharedCommand(){
	delete order;
	delete input;
	delete cCurve;
	delete read;
	delete util;
}

//**********************************************************************************************************************

int CollectSharedCommand::execute(){
	try {
		int count = 1;
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }
		
		if (format == "sharedfile") {
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			order = input->getSharedOrderVector();
		}else {
			//you are using a list and a groupfile
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			input = globaldata->ginput;
			SharedList = globaldata->gSharedList;
			order = SharedList->getSharedOrderVector();
		}
		set<string> orderList;
		
		//set users groups
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups, "collect");
		util->updateGroupIndex(globaldata->Groups, globaldata->gGroupmap->groupIndex);

		while(order != NULL){
			orderList.insert(order->getLabel());
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				convert(globaldata->getFreq(), freq);
				cCurve->getSharedCurve(freq);
			
				delete cCurve;
			
				cout << order->getLabel() << '\t' << count << endl;
			}
			
			//get next line to process
			if (format == "sharedfile") {
				order = input->getSharedOrderVector();
			}else {
				//you are using a list and a groupfile
				SharedList = input->getSharedListVector(); //get new list vector to process
				if (SharedList != NULL) {
					order = SharedList->getSharedOrderVector(); //gets new order vector with group info.
				}else {
					break;
				}
			}
			
			count++;
		}
		set<string>::iterator i;
		for(i = globaldata->labels.begin(); i != globaldata->labels.end(); ++i)
			if(orderList.count(*i) == 0)
				cout << "'" << *i << "'" << " is not a valid label.\n";
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}	
		
		//reset groups parameter
		globaldata->Groups.clear();  globaldata->setGroups("");
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectSharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectSharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************/

