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
#include "nseqs.h"

//**********************************************************************************************************************
vector<string> VennCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(pshared);	
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pcalc("calc", "String", "", "", "", "", "",false,false); parameters.push_back(pcalc);
		CommandParameter pabund("abund", "Number", "", "10", "", "", "",false,false); parameters.push_back(pabund);
		CommandParameter pnseqs("nseqs", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pnseqs);
		CommandParameter ppermute("permute", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(ppermute);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string VennCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The venn command parameters are list, shared, groups, calc, abund, nseqs, permute and label.   shared, relabund, list, rabund or sabund is required unless you have a valid current file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included in your venn diagram, you may only use a maximum of 4 groups.\n";
		helpString += "The group names are separated by dashes. The label allows you to select what distance levels you would like a venn diagram created for, and are also separated by dashes.\n";
		helpString += "The venn command should be in the following format: venn(groups=yourGroups, calc=yourCalcs, label=yourLabels, abund=yourAbund).\n";
		helpString += "Example venn(groups=A-B-C, calc=sharedsobs-sharedchao, abund=20).\n";
		helpString += "The default value for groups is all the groups in your groupfile up to 4, and all labels in your inputfile will be used.\n";
		helpString += "The default value for calc is sobs if you have only read a list file or if you have selected only one group, and sharedsobs if you have multiple groups.\n";
		helpString += "The default available estimators for calc are sobs, chao and ace if you have only read a list file, and sharedsobs, sharedchao and sharedace if you have read a shared file.\n";
		helpString += "The nseqs parameter will output the number of sequences represented by the otus in the picture, default=F.\n";
		helpString += "If you have more than 4 groups, the permute parameter will find all possible combos of 4 of your groups and create pictures for them, default=F.\n";
		helpString += "The only estimators available four 4 groups are sharedsobs and sharedchao.\n";
		helpString += "The venn command outputs a .svg file for each calculator you specify at each distance you choose.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
VennCommand::VennCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["svg"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "VennCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

VennCommand::VennCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; m->setListFile(listfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; m->setSharedFile(sharedfile); }
			
			if ((sharedfile == "") && (listfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					listfile = m->getListFile(); 
					if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list or shared file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile);		}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
			
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
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
			
			string temp;
			temp = validParameter.validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			convert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "nseqs", false);		if (temp == "not found"){	temp = "f";				}
			nseqs = m->isTrue(temp); 

			temp = validParameter.validFile(parameters, "permute", false);			if (temp == "not found"){	temp = "f";				}
			perm = m->isTrue(temp); 

		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "VennCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int VennCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
	
		ValidCalculators validCalculator;
					
		if (format == "list") {
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("vennsingle", Estimators[i]) == true) { 
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
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("vennshared", Estimators[i]) == true) { 
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
			
		//if the users entered no valid calculators don't execute command
		if (vennCalculators.size() == 0) { m->mothurOut("No valid calculators given, please correct."); m->mothurOutEndLine(); return 0;  }
		
		venn = new Venn(outputDir, nseqs, inputfile); 
		input = new InputData(inputfile, format);
		
		string lastLabel;
		
		if (format == "sharedfile") {
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup[0]->getLabel();
			
			if ((lookup.size() > 4) && (perm)) { combosOfFour = findCombinations(lookup.size()); }
		}else if (format == "list") {
			sabund = input->getSAbundVector();
			lastLabel = sabund->getLabel();
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		if (format != "list") {	
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->control_pressed) {
					for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					m->clearGroups(); delete venn; delete input;
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
					return 0;
				}

				if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					
					if ((lookup.size() > 4) && (!perm)){
						m->mothurOut("Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile. If you set perm=t, I will find all possible combos of 4 groups."); m->mothurOutEndLine();
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					
						vector<string> outfilenames = venn->getPic(lookup, vennCalculators);
						for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  } }

					}else if ((lookup.size() > 4) && (perm)) {
						set< set<int> >::iterator it3;
						set<int>::iterator it2;
						for (it3 = combosOfFour.begin(); it3 != combosOfFour.end(); it3++) {  
			
							set<int> poss = *it3;
							vector<SharedRAbundVector*> subset;
							for (it2 = poss.begin(); it2 != poss.end(); it2++) {   subset.push_back(lookup[*it2]);   }
							
							vector<string> outfilenames = venn->getPic(subset, vennCalculators);
							for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]); }  }
						}		
					}else {
						vector<string> outfilenames = venn->getPic(lookup, vennCalculators);
						for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }
					}					
				}
				
				if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup[0]->getLabel();
					
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());

					if ((lookup.size() > 4) && (!perm)){
						m->mothurOut("Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile. If you set perm=t, I will find all possible combos of 4 groups."); m->mothurOutEndLine();
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					
						vector<string> outfilenames = venn->getPic(lookup, vennCalculators);
						for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }

					}else if ((lookup.size() > 4) && (perm)) {
						set< set<int> >::iterator it3;
						set<int>::iterator it2;
						for (it3 = combosOfFour.begin(); it3 != combosOfFour.end(); it3++) {  
			
							set<int> poss = *it3;
							vector<SharedRAbundVector*> subset;
							for (it2 = poss.begin(); it2 != poss.end(); it2++) {   subset.push_back(lookup[*it2]);   }
							
							vector<string> outfilenames = venn->getPic(subset, vennCalculators);
							for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }
						}		
					}else {
						vector<string> outfilenames = venn->getPic(lookup, vennCalculators);
						for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }
					}
										
					//restore real lastlabel to save below
					lookup[0]->setLabel(saveLabel);
				}
				
				
				lastLabel = lookup[0]->getLabel();	
						
				//get next line to process
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors();
			}
			
			if (m->control_pressed) {
					for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
					m->clearGroups(); delete venn; delete input; 
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
					return 0;
			}

			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				m->mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
					needToRun = true;
				}else {
					m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
				}
			}
		
			//run last label if you need to
			if (needToRun == true)  {
					for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i]; }  } 
					lookup = input->getSharedRAbundVectors(lastLabel);

					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());

					if ((lookup.size() > 4) && (!perm)){
						m->mothurOut("Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile. If you set perm=t, I will find all possible combos of 4 groups."); m->mothurOutEndLine();
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					
						vector<string> outfilenames = venn->getPic(lookup, vennCalculators);
						for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]); }  }

					}else if ((lookup.size() > 4) && (perm)) {
						set< set<int> >::iterator it3;
						set<int>::iterator it2;
						for (it3 = combosOfFour.begin(); it3 != combosOfFour.end(); it3++) {  
			
							set<int> poss = *it3;
							vector<SharedRAbundVector*> subset;
							for (it2 = poss.begin(); it2 != poss.end(); it2++) {   subset.push_back(lookup[*it2]);   }
							
							vector<string> outfilenames = venn->getPic(subset, vennCalculators);
							for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]); }  }
						}		
					}else {
						vector<string> outfilenames = venn->getPic(lookup, vennCalculators);
						for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }
					}
					
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			}
		

			//reset groups parameter
			m->clearGroups();  
			
			if (m->control_pressed) {
					m->clearGroups(); delete venn; delete input;
					for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
					return 0;
			}

			
		}else{
		
			while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->control_pressed) {
					for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
					delete sabund; delete venn; delete input;
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
					return 0;
				}
		
				if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
	
					m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
					vector<string> outfilenames = venn->getPic(sabund, vennCalculators);
					for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }

					
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
				}
				
				if ((m->anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = sabund->getLabel();
				
					delete sabund;
					sabund = input->getSAbundVector(lastLabel);
					
					m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
					vector<string> outfilenames = venn->getPic(sabund, vennCalculators);
					for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }

					
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
					
					//restore real lastlabel to save below
					sabund->setLabel(saveLabel);
				}		
				
				lastLabel = sabund->getLabel();		
				
				delete sabund;
				sabund = input->getSAbundVector();
			}
			
			if (m->control_pressed) {
					for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
					delete venn; delete input;
					return 0;
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				m->mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
					needToRun = true;
				}else {
					m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
				}
			}
		
			//run last label if you need to
			if (needToRun == true)  {
				if (sabund != NULL) {	delete sabund;	}
				sabund = input->getSAbundVector(lastLabel);
					
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				vector<string> outfilenames = venn->getPic(sabund, vennCalculators);
				for(int i = 0; i < outfilenames.size(); i++) { if (outfilenames[i] != "control" ) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]);  }  }

				delete sabund;
					
			}
			
			if (m->control_pressed) {
					delete venn; delete input;
					for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
					return 0;
			}
		}
		
		for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
		delete venn; delete input;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
//returns a vector of sets containing the 4 group combinations
set< set<int> > VennCommand::findCombinations(int lookupSize){
	try {
		set< set<int> > combos;
		
		set<int> possibles;
		for (int i = 0; i < lookupSize; i++) {  possibles.insert(i);  }
		
		getCombos(possibles, combos);
		
		return combos;
		
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "findCombinations");
		exit(1);
	}
}
//**********************************************************************************************************************
//recusively finds combos of 4
int VennCommand::getCombos(set<int> possibles, set< set<int> >& combos){
	try {
		
		if (possibles.size() == 4) { //done
			if (combos.count(possibles) == 0) { //no dups
				combos.insert(possibles);
			}
		}else { //we still have work to do
			set<int>::iterator it;
			set<int>::iterator it2;
			for (it = possibles.begin(); it != possibles.end(); it++) {  
				
				set<int> newPossibles;
				for (it2 = possibles.begin(); it2 != possibles.end(); it2++) {  //all possible combos of one length smaller
					if (*it != *it2) { 
						newPossibles.insert(*it2);
					}
				}
				getCombos(newPossibles, combos);
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "VennCommand", "getCombos");
		exit(1);
	}
}

//**********************************************************************************************************************
