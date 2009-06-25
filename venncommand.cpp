/*
 *  venncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "venncommand.h"
#include "ace.h"
#include "sobs.h"
#include "chao1.h"
//#include "jackknife.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"


//**********************************************************************************************************************

VennCommand::VennCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		lines.clear();
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","line","label","calc", "abund"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getSharedFile() == "")) {
				mothurOut("You must read a list, or a list and a group, or a shared before you can use the venn command."); mothurOutEndLine(); abort = true; 
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
			else if ((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			format = globaldata->getFormat();
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { 
				if(format == "list") { calc = "sobs"; }
				else { calc = "sharedsobs"; }
			}
			else { 
				 if (calc == "default")  {  
					if(format == "list") { calc = "sobs"; }
					else { calc = "sharedsobs"; }
				}
			}
			splitAtDash(calc, Estimators);
			
			string temp;
			temp = validParameter.validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 

			if (abort == false) {
				validCalculator = new ValidCalculators();
		
				int i;
				
				if (format == "list") {
					for (i=0; i<Estimators.size(); i++) {
						if (validCalculator->isValidCalculator("vennsingle", Estimators[i]) == true) { 
							if (Estimators[i] == "sobs") { 
								vennCalculators.push_back(new Sobs());
							}else if (Estimators[i] == "chao") { 
								vennCalculators.push_back(new Chao1());
							}else if (Estimators[i] == "ace") {
								if(abund < 5)
									abund = 10;
								vennCalculators.push_back(new Ace(abund));
							}
						}
					}
				}else {
					for (i=0; i<Estimators.size(); i++) {
						if (validCalculator->isValidCalculator("vennshared", Estimators[i]) == true) { 
							if (Estimators[i] == "sharedsobs") { 
								vennCalculators.push_back(new SharedSobsCS());
							}else if (Estimators[i] == "sharedchao") { 
								vennCalculators.push_back(new SharedChao1());
							}else if (Estimators[i] == "sharedace") { 
								vennCalculators.push_back(new SharedAce());
							}
						}
					}
				}
				
				venn = new Venn();
			}
		}

		
				
	}
	catch(exception& e) {
		errorOut(e, "VennCommand", "VennCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void VennCommand::help(){
	try {
		mothurOut("The venn command can only be executed after a successful read.otu command.\n");
		mothurOut("The venn command parameters are groups, calc, abund, line and label.  No parameters are required, but you may not use line and label at the same time.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included in your venn diagram, you may only use a maximum of 4 groups.\n");
		mothurOut("The group names are separated by dashes. The line and label allow you to select what distance levels you would like a venn diagram created for, and are also separated by dashes.\n");
		mothurOut("The venn command should be in the following format: venn(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels, abund=yourAbund).\n");
		mothurOut("Example venn(groups=A-B-C, line=1-3-5, calc=sharedsobs-sharedchao, abund=20).\n");
		mothurOut("The default value for groups is all the groups in your groupfile up to 4, and all lines in your inputfile will be used.\n");
		mothurOut("The default value for calc is sobs if you have only read a list file or if you have selected only one group, and sharedsobs if you have multiple groups.\n");
		mothurOut("The default available estimators for calc are sobs, chao and ace if you have only read a list file, and sharedsobs, sharedchao and sharedace if you have read a list and group file or a shared file.\n");
		mothurOut("The only estmiator available four 4 groups is sharedsobs.\n");
		mothurOut("The venn command outputs a .svg file for each calculator you specify at each distance you choose.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "VennCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

VennCommand::~VennCommand(){
	if (abort == false) {
		delete input; globaldata->ginput = NULL;
		delete read;
		delete venn;
		globaldata->sabund = NULL;
	}
	
}

//**********************************************************************************************************************

int VennCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		int count = 1;
		string lastLabel;
		
		//if the users entered no valid calculators don't execute command
		if (vennCalculators.size() == 0) { return 0; }
		
		if (format == "sharedfile") {
			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup[0]->getLabel();
		}else if (format == "list") {
			//you are using just a list file and have only one group
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			sabund = globaldata->sabund;
			lastLabel = sabund->getLabel();
			input = globaldata->ginput;
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;
		
		if (format != "list") {	
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {

				if(allLines == 1 || lines.count(count) == 1 || labels.count(lookup[0]->getLabel()) == 1){			
					mothurOut(lookup[0]->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					userLines.erase(count);
					
					if (lookup.size() > 4) {
						mothurOut("Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile."); mothurOutEndLine();
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					}
					venn->getPic(lookup, vennCalculators);
				}
				
				if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					mothurOut(lookup[0]->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());

					if (lookup.size() > 4) {
						mothurOut("Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile."); mothurOutEndLine();
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					}				
					venn->getPic(lookup, vennCalculators);
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
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());

					if (lookup.size() > 4) {
						mothurOut("Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile."); mothurOutEndLine();
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					}				
					venn->getPic(lookup, vennCalculators);
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			}
		

			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else{
		
			while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
				if(allLines == 1 || lines.count(count) == 1 || labels.count(sabund->getLabel()) == 1){			
	
					mothurOut(sabund->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					venn->getPic(sabund, vennCalculators);
					
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
					userLines.erase(count);
				}
				
				if ((anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					delete sabund;
					sabund = input->getSAbundVector(lastLabel);
					
					mothurOut(sabund->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					venn->getPic(sabund, vennCalculators);
					
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
				}		
				
				lastLabel = sabund->getLabel();		
				
				delete sabund;
				sabund = input->getSAbundVector();
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
				delete sabund;
				sabund = input->getSAbundVector(lastLabel);
					
				mothurOut(sabund->getLabel() + "\t" + toString(count)); mothurOutEndLine();
				venn->getPic(sabund, vennCalculators);
				delete sabund;
					
			}
			
		}
		
		for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "VennCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
