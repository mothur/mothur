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
		lines.clear();
		labels.clear();
		Groups.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","line","label", "calc"};
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
				 cout << "You must read a list and a group, or a shared before you can use the heatmap.sim command." << endl; abort = true; 
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
			if ((line != "") && (label != "")) { cout << "You cannot use both the line and label parameters at the same time. " << endl; abort = true; }
			//if the user has not specified any line or labels use the ones from read.otu
			else if ((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
			
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
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
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapSimCommand class Function HeatMapSimCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapSimCommand class function HeatMapSimCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

void HeatMapSimCommand::help(){
	try {
		cout << "The heatmap.sim command can only be executed after a successful read.otu command." << "\n";
		cout << "The heatmap.sim command parameters are groups, calc, line and label.  No parameters are required, but you may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like a heatmap created for, and are also separated by dashes." << "\n";
		cout << "The heatmap.sim command should be in the following format: heatmap.sim(groups=yourGroups, calc=yourCalc, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example heatmap.sim(groups=A-B-C, line=1-3-5, calc=jabund)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and all lines in your inputfile will be used." << "\n";
		validCalculator->printCalc("heat", cout);
		cout << "The default value for calc is jclass-thetayc." << "\n";
		cout << "The heatmap.sim command outputs a .svg file for each calculator you choose at each line or label you specify." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapSimCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapSimCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

HeatMapSimCommand::~HeatMapSimCommand(){
	delete input;
	delete read;
	delete heatmap;
	delete validCalculator;
}

//**********************************************************************************************************************

int HeatMapSimCommand::execute(){
	try {
	
		if (abort == true)  { return 0; }
		
		int count = 1;	
		
		//if the users entered no valid calculators don't execute command
		if (heatCalculators.size() == 0) { cout << "No valid calculators." << endl; return 0; }
		
		//you have groups
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		vector<SharedRAbundVector*> lastLookup = lookup;
		
		if (lookup.size() < 2) { cout << "You have not provided enough valid groups.  I cannot run the command." << endl; return 0;}
				
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;

		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
			if(allLines == 1 || lines.count(count) == 1 || labels.count(lookup[0]->getLabel()) == 1){			
	
				cout << lookup[0]->getLabel() << '\t' << count << endl;
				heatmap->getPic(lookup, heatCalculators);
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				userLines.erase(count);
			}
				
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLookup[0]->getLabel()) != 1)) {
				cout << lastLookup[0]->getLabel() << '\t' << count << endl;
				heatmap->getPic(lastLookup, heatCalculators);
					
				processedLabels.insert(lastLookup[0]->getLabel());
				userLabels.erase(lastLookup[0]->getLabel());
			}
				
			//prevent memory leak
			if (count != 1) { for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  } }
			lastLookup = lookup;			

			//get next line to process
			lookup = input->getSharedRAbundVectors();				
			count++;
		}
			
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastLookup[0]->getLabel()) != 1) {
				cout << ". I will use " << lastLookup[0]->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastLookup[0]->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastLookup[0]->getLabel() << '\t' << count << endl;
			heatmap->getPic(lastLookup, heatCalculators);
		}
		
		for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  }
			
		//reset groups parameter
		globaldata->Groups.clear();  
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapSimCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapSimCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
