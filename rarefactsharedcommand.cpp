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
		lines.clear();
		labels.clear();
		Estimators.clear();
		Groups.clear();
				
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"iters","line","label","calc","groups"};
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
		mothurOut("The rarefaction.shared command parameters are label, line, iters, groups and calc.  No parameters are required, but you may not use \n");
		mothurOut("both the line and label parameters at the same time. The rarefaction command should be in the following format: \n");
		mothurOut("rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, calc=yourEstimators, groups=yourGroups).\n");
		mothurOut("Example rarefaction.shared(label=unique-.01-.03, line=0-5-10, iters=10000, groups=B-C, calc=sharedobserved).\n");
		mothurOut("The default values for iters is 1000, freq is 100, and calc is sharedobserved which calculates the shared rarefaction curve for the observed richness.\n");
		mothurOut("The default value for groups is all the groups in your groupfile.\n");
		validCalculator->printCalc("sharedrarefaction", cout);
		mothurOut("The label and line parameters are used to analyze specific lines in your input.\n");
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
		
		int count = 1;
		
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
		set<int> userLines = lines;
	
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(allLines == 1 || lines.count(count) == 1 || labels.count(lookup[0]->getLabel()) == 1){
				
				rCurve = new Rarefact(lookup, rDisplays);
				rCurve->getSharedCurve(freq, nIters);
				delete rCurve;
			
				mothurOut(lookup[0]->getLabel() + "\t" + toString(count)); mothurOutEndLine();
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					mothurOut(lookup[0]->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					rCurve = new Rarefact(lookup, rDisplays);
					rCurve->getSharedCurve(freq, nIters);
					delete rCurve;

					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
			}
				
			
			lastLabel = lookup[0]->getLabel();
			
			//get next line to process
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
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
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors(lastLabel);

			mothurOut(lookup[0]->getLabel() + "\t" + toString(count)); mothurOutEndLine();
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
