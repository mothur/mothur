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
#include "whittaker.h"



//**********************************************************************************************************************

CollectSharedCommand::CollectSharedCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		Estimators.clear();
		Groups.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"freq","label","calc","groups"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters=parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getSharedFile() == "") {
				if (globaldata->getListFile() == "") { mothurOut("You must read a list and a group, or a shared before you can use the collect.shared command."); mothurOutEndLine(); abort = true; }
				else if (globaldata->getGroupFile() == "") { mothurOut("You must read a list and a group, or a shared before you can use the collect.shared command."); mothurOutEndLine(); abort = true; }
			}

			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking..
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
			if (calc == "not found") { calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			else { 
				 if (calc == "default")  {  calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			}
			splitAtDash(calc, Estimators);
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
			}
			globaldata->Groups = Groups;
			
			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, freq); 
						
			if (abort == false) {
			
				string fileNameRoot = getRootName(globaldata->inputFileName);
				format = globaldata->getFormat();
				int i;
				
				validCalculator = new ValidCalculators();
				util = new SharedUtil();
				
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("shared", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedchao") { 
							cDisplays.push_back(new CollectDisplay(new SharedChao1(), new SharedOneColumnFile(fileNameRoot+"shared.chao")));
						}else if (Estimators[i] == "sharedsobs") { 
							cDisplays.push_back(new CollectDisplay(new SharedSobsCS(), new SharedOneColumnFile(fileNameRoot+"shared.sobs")));
						}else if (Estimators[i] == "sharedace") { 
							cDisplays.push_back(new CollectDisplay(new SharedAce(), new SharedOneColumnFile(fileNameRoot+"shared.ace")));
						}else if (Estimators[i] == "jabund") { 	
							cDisplays.push_back(new CollectDisplay(new JAbund(), new SharedOneColumnFile(fileNameRoot+"jabund")));
						}else if (Estimators[i] == "sorabund") { 
							cDisplays.push_back(new CollectDisplay(new SorAbund(), new SharedOneColumnFile(fileNameRoot+"sorabund")));
						}else if (Estimators[i] == "jclass") { 
							cDisplays.push_back(new CollectDisplay(new Jclass(), new SharedOneColumnFile(fileNameRoot+"jclass")));
						}else if (Estimators[i] == "sorclass") { 
							cDisplays.push_back(new CollectDisplay(new SorClass(), new SharedOneColumnFile(fileNameRoot+"sorclass")));
						}else if (Estimators[i] == "jest") { 
							cDisplays.push_back(new CollectDisplay(new Jest(), new SharedOneColumnFile(fileNameRoot+"jest")));
						}else if (Estimators[i] == "sorest") { 
							cDisplays.push_back(new CollectDisplay(new SorEst(), new SharedOneColumnFile(fileNameRoot+"sorest")));
						}else if (Estimators[i] == "thetayc") { 
							cDisplays.push_back(new CollectDisplay(new ThetaYC(), new SharedOneColumnFile(fileNameRoot+"thetayc")));
						}else if (Estimators[i] == "thetan") { 
							cDisplays.push_back(new CollectDisplay(new ThetaN(), new SharedOneColumnFile(fileNameRoot+"thetan")));
						}else if (Estimators[i] == "kstest") { 
							cDisplays.push_back(new CollectDisplay(new KSTest(), new SharedOneColumnFile(fileNameRoot+"kstest")));
						}else if (Estimators[i] == "whittaker") { 
							cDisplays.push_back(new CollectDisplay(new Whittaker(), new SharedOneColumnFile(fileNameRoot+"whittaker")));
						}else if (Estimators[i] == "sharednseqs") { 
							cDisplays.push_back(new CollectDisplay(new SharedNSeqs(), new SharedOneColumnFile(fileNameRoot+"shared.nseqs")));
						}else if (Estimators[i] == "ochiai") { 
							cDisplays.push_back(new CollectDisplay(new Ochiai(), new SharedOneColumnFile(fileNameRoot+"ochiai")));
						}else if (Estimators[i] == "anderberg") { 
							cDisplays.push_back(new CollectDisplay(new Anderberg(), new SharedOneColumnFile(fileNameRoot+"anderberg")));
						}else if (Estimators[i] == "skulczynski") { 
							cDisplays.push_back(new CollectDisplay(new Kulczynski(), new SharedOneColumnFile(fileNameRoot+"kulczynski")));
						}else if (Estimators[i] == "kulczynskicody") { 
							cDisplays.push_back(new CollectDisplay(new KulczynskiCody(), new SharedOneColumnFile(fileNameRoot+"kulczynskicody")));
						}else if (Estimators[i] == "lennon") { 
							cDisplays.push_back(new CollectDisplay(new Lennon(), new SharedOneColumnFile(fileNameRoot+"lennon")));
						}else if (Estimators[i] == "morisitahorn") { 
							cDisplays.push_back(new CollectDisplay(new MorHorn(), new SharedOneColumnFile(fileNameRoot+"morisitahorn")));
						}else if (Estimators[i] == "braycurtis") { 
							cDisplays.push_back(new CollectDisplay(new BrayCurtis(), new SharedOneColumnFile(fileNameRoot+"braycurtis")));
						}
					}
				}	
			}
		}

	}
	catch(exception& e) {
		errorOut(e, "CollectSharedCommand", "CollectSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void CollectSharedCommand::help(){
	try {
		mothurOut("The collect.shared command can only be executed after a successful read.otu command.\n");
		mothurOut("The collect.shared command parameters are label, freq, calc and groups.  No parameters are required \n");
		mothurOut("The collect.shared command should be in the following format: \n");
		mothurOut("collect.shared(label=yourLabel, freq=yourFreq, calc=yourEstimators, groups=yourGroups).\n");
		mothurOut("Example collect.shared(label=unique-.01-.03, freq=10, groups=B-C, calc=sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan).\n");
		mothurOut("The default values for freq is 100 and calc are sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan.\n");
		mothurOut("The default value for groups is all the groups in your groupfile.\n");
		validCalculator->printCalc("shared", cout);
		mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile).\n\n");
		
	}
	catch(exception& e) {
		errorOut(e, "CollectSharedCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

CollectSharedCommand::~CollectSharedCommand(){
	if (abort == false) {
		delete input; globaldata->ginput = NULL;
		delete read;
		delete util;
		delete validCalculator;
		globaldata->gorder = NULL;
	}
}

//**********************************************************************************************************************

int CollectSharedCommand::execute(){
	try {
		
		if (abort == true) {	return 0;	}
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }
		
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		order = input->getSharedOrderVector();
		string lastLabel = order->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
				
		//set users groups
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups, "collect");
		util->updateGroupIndex(globaldata->Groups, globaldata->gGroupmap->groupIndex);

		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {

			if(allLines == 1 || labels.count(order->getLabel()) == 1){
				
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				cCurve->getSharedCurve(freq);
				delete cCurve;
			
				mothurOut(order->getLabel()); mothurOutEndLine();
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			}
			
			//you have a label the user want that is smaller than this label and the last label has not already been processed
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				delete order;
				order = input->getSharedOrderVector(lastLabel);
				
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				cCurve->getSharedCurve(freq);
				delete cCurve;
			
				mothurOut(order->getLabel()); mothurOutEndLine();
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			}
			
			
			lastLabel = order->getLabel();			
			
			//get next line to process
			delete order;
			order = input->getSharedOrderVector();
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
			if (order != NULL) {  delete order;  }
			order = input->getSharedOrderVector(lastLabel);

			cCurve = new Collect(order, cDisplays);
			cCurve->getSharedCurve(freq);
			delete cCurve;
			
			mothurOut(order->getLabel()); mothurOutEndLine();
			delete order;
		}
		
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}	
		
		//reset groups parameter
		globaldata->Groups.clear(); 
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "CollectSharedCommand", "execute");
		exit(1);
	}
}

/***********************************************************/
