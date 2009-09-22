/*
 *  heatmapsimcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmapsimcommand.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"


//**********************************************************************************************************************

HeatMapSimCommand::HeatMapSimCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		Groups.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label", "calc"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getSharedFile() == "") {
				 mothurOut("You must read a list and a group, or a shared before you can use the heatmap.sim command."); mothurOutEndLine(); abort = true; 
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
			if (label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
			
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jest-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jest-thetayc";  }
			}
			splitAtDash(calc, Estimators);
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			if (abort == false) {
				validCalculator = new ValidCalculators();
				heatmap = new HeatMapSim();
			
				int i;
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("heat", Estimators[i]) == true) { 
						if (Estimators[i] == "jabund") { 	
							heatCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							heatCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							heatCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							heatCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							heatCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							heatCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							heatCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							heatCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "morisitahorn") { 
							heatCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							heatCalculators.push_back(new BrayCurtis());
						}
					}
				}
				
			}
		}

				

	}
	catch(exception& e) {
		errorOut(e, "HeatMapSimCommand", "HeatMapSimCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void HeatMapSimCommand::help(){
	try {
		mothurOut("The heatmap.sim command can only be executed after a successful read.otu command.\n");
		mothurOut("The heatmap.sim command parameters are groups, calc and label.  No parameters are required.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap.\n");
		mothurOut("The group names are separated by dashes. The label parameter allows you to select what distance levels you would like a heatmap created for, and is also separated by dashes.\n");
		mothurOut("The heatmap.sim command should be in the following format: heatmap.sim(groups=yourGroups, calc=yourCalc, label=yourLabels).\n");
		mothurOut("Example heatmap.sim(groups=A-B-C, calc=jabund).\n");
		mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		validCalculator->printCalc("heat", cout);
		mothurOut("The default value for calc is jclass-thetayc.\n");
		mothurOut("The heatmap.sim command outputs a .svg file for each calculator you choose at each label you specify.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");

	}
	catch(exception& e) {
		errorOut(e, "HeatMapSimCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

HeatMapSimCommand::~HeatMapSimCommand(){
	if (abort == false) {
		delete input;  globaldata->ginput = NULL;
		delete read;
		delete heatmap;
		delete validCalculator;
	}
}

//**********************************************************************************************************************

int HeatMapSimCommand::execute(){
	try {
	
		if (abort == true)  { return 0; }
		
		//if the users entered no valid calculators don't execute command
		if (heatCalculators.size() == 0) { mothurOut("No valid calculators."); mothurOutEndLine(); return 0; }
		
		//you have groups
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		if (lookup.size() < 2) { mothurOut("You have not provided enough valid groups.  I cannot run the command."); mothurOutEndLine(); return 0;}
				
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
		
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
	
				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				heatmap->getPic(lookup, heatCalculators);
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
				
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
			
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);				

				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				heatmap->getPic(lookup, heatCalculators);
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
				
			//prevent memory leak
			 
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
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i]; } } 
			lookup = input->getSharedRAbundVectors(lastLabel);				

			mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
			heatmap->getPic(lookup, heatCalculators);
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
			
		//reset groups parameter
		globaldata->Groups.clear();  
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "HeatMapSimCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
