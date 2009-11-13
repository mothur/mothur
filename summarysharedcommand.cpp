/*
 *  summarysharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarysharedcommand.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
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

SummarySharedCommand::SummarySharedCommand(string option){
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
			string Array[] =  {"label","calc","groups","all"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getSharedFile() == "") {
				 mothurOut("You must read a list and a group, or a shared before you can use the summary.shared command."); mothurOutEndLine(); abort = true; 
			}
			
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
			if (calc == "not found") { calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			else { 
				 if (calc == "default")  {  calc = "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan";  }
			}
			splitAtDash(calc, Estimators);
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			string temp = validParameter.validFile(parameters, "all", false);				if (temp == "not found") { temp = "true"; }
			all = isTrue(temp);
			
			if (abort == false) {
			
				validCalculator = new ValidCalculators();
				int i;
				
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("sharedsummary", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedsobs") { 
							sumCalculators.push_back(new SharedSobsCS());
						}else if (Estimators[i] == "sharedchao") { 
							sumCalculators.push_back(new SharedChao1());
						}else if (Estimators[i] == "sharedace") { 
							sumCalculators.push_back(new SharedAce());
						}else if (Estimators[i] == "jabund") { 	
							sumCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							sumCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							sumCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							sumCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							sumCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							sumCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							sumCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							sumCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "kstest") { 
							sumCalculators.push_back(new KSTest());
						}else if (Estimators[i] == "sharednseqs") { 
							sumCalculators.push_back(new SharedNSeqs());
						}else if (Estimators[i] == "ochiai") { 
							sumCalculators.push_back(new Ochiai());
						}else if (Estimators[i] == "anderberg") { 
							sumCalculators.push_back(new Anderberg());
						}else if (Estimators[i] == "kulczynski") { 
							sumCalculators.push_back(new Kulczynski());
						}else if (Estimators[i] == "kulczynskicody") { 
							sumCalculators.push_back(new KulczynskiCody());
						}else if (Estimators[i] == "lennon") { 
							sumCalculators.push_back(new Lennon());
						}else if (Estimators[i] == "morisitahorn") { 
							sumCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							sumCalculators.push_back(new BrayCurtis());
						}else if (Estimators[i] == "whittaker") { 
							sumCalculators.push_back(new Whittaker());
						}
					}
				}
				
				outputFileName = ((getRootName(globaldata->inputFileName)) + "shared.summary");
				openOutputFile(outputFileName, outputFileHandle);
				mult = false;
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "SummarySharedCommand", "SummarySharedCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void SummarySharedCommand::help(){
	try {
		mothurOut("The summary.shared command can only be executed after a successful read.otu command.\n");
		mothurOut("The summary.shared command parameters are label, calc and all.  No parameters are required.\n");
		mothurOut("The summary.shared command should be in the following format: \n");
		mothurOut("summary.shared(label=yourLabel, calc=yourEstimators, groups=yourGroups).\n");
		mothurOut("Example summary.shared(label=unique-.01-.03, groups=B-C, calc=sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan).\n");
		validCalculator->printCalc("sharedsummary", cout);
		mothurOut("The default value for calc is sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan\n");
		mothurOut("The default value for groups is all the groups in your groupfile.\n");
		mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		mothurOut("The all parameter is used to specify if you want the estimate of all your groups together.  This estimate can only be made for sharedsobs and sharedchao calculators. The default is true.\n");
		mothurOut("If you use sharedchao and run into memory issues, set all to false. \n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "SummarySharedCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

SummarySharedCommand::~SummarySharedCommand(){
	if (abort == false) {
		delete read;
		delete validCalculator;
	}
}

//**********************************************************************************************************************

int SummarySharedCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
	
		//if the users entered no valid calculators don't execute command
		if (sumCalculators.size() == 0) { return 0; }
		//check if any calcs can do multiples
		else{
			if (all){ 
				for (int i = 0; i < sumCalculators.size(); i++) {
					if (sumCalculators[i]->getMultiple() == true) { mult = true; }
				}
			}
		}
		
		//read first line
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//output estimator names as column headers
		outputFileHandle << "label" <<'\t' << "comparison" << '\t'; 
		for(int i=0;i<sumCalculators.size();i++){
			outputFileHandle << '\t' << sumCalculators[i]->getName();
		}
		outputFileHandle << endl;
		
		//create file and put column headers for multiple groups file
		if (mult == true) {
			outAllFileName = ((getRootName(globaldata->inputFileName)) + "sharedmultiple.summary");
			openOutputFile(outAllFileName, outAll);
			
			outAll << "label" <<'\t' << "comparison" << '\t'; 
			for(int i=0;i<sumCalculators.size();i++){
				if (sumCalculators[i]->getMultiple() == true) { 
					outAll << '\t' << sumCalculators[i]->getName();
				}
			}
			outAll << endl;
		}
		
		if (lookup.size() < 2) { 
			mothurOut("I cannot run the command without at least 2 valid groups."); 
			for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
			
			//close files and clean up
			outputFileHandle.close();  remove(outputFileName.c_str());
			if (mult == true) {  outAll.close();  remove(outAllFileName.c_str());  }
			return 0;
		//if you only have 2 groups you don't need a .sharedmultiple file
		}else if ((lookup.size() == 2) && (mult == true)) { 
			mult = false;
			outAll.close();  
			remove(outAllFileName.c_str());
		}
					
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
			
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
		
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup[0]->getLabel();
					
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					process(lookup);
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					
					//restore real lastlabel to save below
					lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();			
				
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
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
				for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i];	} } 
				lookup = input->getSharedRAbundVectors(lastLabel);

				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				process(lookup);
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		

		//reset groups parameter
		globaldata->Groups.clear();  
		
		//close files
		outputFileHandle.close();
		if (mult == true) {  outAll.close();  }
		
		for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
		
		delete input;  globaldata->ginput = NULL;

		return 0;
	}
	catch(exception& e) {
		errorOut(e, "SummarySharedCommand", "execute");
		exit(1);
	}
}

/***********************************************************/
void SummarySharedCommand::process(vector<SharedRAbundVector*> thisLookup) {
	try {
				//loop through calculators and add to file all for all calcs that can do mutiple groups
				if (mult == true) {
					//output label
					outAll << thisLookup[0]->getLabel() << '\t';
					
					//output groups names
					string outNames = "";
					for (int j = 0; j < thisLookup.size(); j++) {
						outNames += thisLookup[j]->getGroup() +  "-";
					}
					outNames = outNames.substr(0, outNames.length()-1); //rip off extra '-';
					outAll << outNames << '\t';
					
					for(int i=0;i<sumCalculators.size();i++){
						if (sumCalculators[i]->getMultiple() == true) { 
							sumCalculators[i]->getValues(thisLookup);
							outAll << '\t';
							sumCalculators[i]->print(outAll);
						}
					}
					outAll << endl;
				}
	
				int n = 1; 
				vector<SharedRAbundVector*> subset;
				for (int k = 0; k < (thisLookup.size() - 1); k++) { // pass cdd each set of groups to commpare
					for (int l = n; l < thisLookup.size(); l++) {
						
						outputFileHandle << thisLookup[0]->getLabel() << '\t';
						
						subset.clear(); //clear out old pair of sharedrabunds
						//add new pair of sharedrabunds
						subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
						
						//sort groups to be alphanumeric
						if (thisLookup[k]->getGroup() > thisLookup[l]->getGroup()) {
							outputFileHandle << (thisLookup[l]->getGroup() +'\t' + thisLookup[k]->getGroup()) << '\t'; //print out groups
						}else{
							outputFileHandle << (thisLookup[k]->getGroup() +'\t' + thisLookup[l]->getGroup()) << '\t'; //print out groups
						}
						
						for(int i=0;i<sumCalculators.size();i++) {

							sumCalculators[i]->getValues(subset); //saves the calculator outputs
							outputFileHandle << '\t';
							sumCalculators[i]->print(outputFileHandle);
						}
						outputFileHandle << endl;
					}
					n++;
				}

	}
	catch(exception& e) {
		errorOut(e, "SummarySharedCommand", "process");
		exit(1);
	}
}

/***********************************************************/
