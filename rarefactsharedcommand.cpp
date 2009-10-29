/*
 *  rarefactsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactsharedcommand.h"
#include "sharedsobs.h"
#include "sharednseqs.h"

//**********************************************************************************************************************

RareFactSharedCommand::RareFactSharedCommand(string option){
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
			string Array[] =  {"iters","label","calc","groups", "jumble"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
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
			if (calc == "not found") { calc = "sharedobserved";  }
			else { 
				 if (calc == "default")  {  calc = "sharedobserved";  }
			}
			splitAtDash(calc, Estimators);
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
			}
			globaldata->Groups = Groups;
			
			string temp;
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, nIters); 
			
			temp = validParameter.validFile(parameters, "jumble", false);			if (temp == "not found") { temp = "T"; }
			if (isTrue(temp)) { jumble = true; }
			else { jumble = false; }
			globaldata->jumble = jumble;
			
			if (abort == false) {
			
				string fileNameRoot = getRootName(globaldata->inputFileName);
//				format = globaldata->getFormat();

				
				validCalculator = new ValidCalculators();
				
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("sharedrarefaction", Estimators[i]) == true) { 
						if (Estimators[i] == "sharedobserved") { 
							rDisplays.push_back(new RareDisplay(new SharedSobs(), new SharedThreeColumnFile(fileNameRoot+"shared.rarefaction", "")));
						}else if (Estimators[i] == "sharednseqs") { 
							rDisplays.push_back(new RareDisplay(new SharedNSeqs(), new SharedThreeColumnFile(fileNameRoot+"shared.r_nseqs", "")));
						}
					}
				}
			}
				
		}

	}
	catch(exception& e) {
		errorOut(e, "RareFactSharedCommand", "RareFactSharedCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void RareFactSharedCommand::help(){
	try {
		mothurOut("The rarefaction.shared command can only be executed after a successful read.otu command.\n");
		mothurOut("The rarefaction.shared command parameters are label, iters, groups, jumble and calc.  No parameters are required.\n");
		mothurOut("The rarefaction command should be in the following format: \n");
		mothurOut("rarefaction.shared(label=yourLabel, iters=yourIters, calc=yourEstimators, jumble=yourJumble, groups=yourGroups).\n");
		mothurOut("Example rarefaction.shared(label=unique-0.01-0.03,  iters=10000, groups=B-C, jumble=T, calc=sharedobserved).\n");
		mothurOut("The default values for iters is 1000, freq is 100, and calc is sharedobserved which calculates the shared rarefaction curve for the observed richness.\n");
		mothurOut("The default value for groups is all the groups in your groupfile, and jumble is true.\n");
		validCalculator->printCalc("sharedrarefaction", cout);
		mothurOut("The label parameter is used to analyze specific labels in your input.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "RareFactSharedCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

RareFactSharedCommand::~RareFactSharedCommand(){
	if (abort == false) {
		delete input;   globaldata->ginput = NULL;
		delete read;
		delete validCalculator;
	}
}

//**********************************************************************************************************************

int RareFactSharedCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		//if the users entered no valid calculators don't execute command
		if (rDisplays.size() == 0) { return 0; }

		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();

		if (lookup.size() < 2) { 
			mothurOut("I cannot run the command without at least 2 valid groups."); 
			for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
			return 0;
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
	
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){
				
				rCurve = new Rarefact(lookup, rDisplays);
				rCurve->getSharedCurve(freq, nIters);
				delete rCurve;
			
				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup[0]->getLabel();
			
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					rCurve = new Rarefact(lookup, rDisplays);
					rCurve->getSharedCurve(freq, nIters);
					delete rCurve;

					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					
					//restore real lastlabel to save below
					lookup[0]->setLabel(saveLabel);
			}
				
			
			lastLabel = lookup[0]->getLabel();
			
			//get next line to process
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
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i]; }  } 
			lookup = input->getSharedRAbundVectors(lastLabel);

			mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
			rCurve = new Rarefact(lookup, rDisplays);
			rCurve->getSharedCurve(freq, nIters);
			delete rCurve;
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
		for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	
		
		//reset groups parameter
		globaldata->Groups.clear();  

		return 0;
	}
	catch(exception& e) {
		errorOut(e, "RareFactSharedCommand", "execute");
		exit(1);
	}
}


//**********************************************************************************************************************
