/*
 *  matrixoutputcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "matrixoutputcommand.h"
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

MatrixOutputCommand::MatrixOutputCommand(string option){
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
			string Array[] =  {"label","calc","groups","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters  = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getSharedFile() == "") {
				if (globaldata->getListFile() == "") { mothurOut("You must read a list and a group, or a shared before you can use the dist.shared command."); mothurOutEndLine(); abort = true; }
				else if (globaldata->getGroupFile() == "") { mothurOut("You must read a list and a group, or a shared before you can use the dist.shared command."); mothurOutEndLine(); abort = true; }
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
				
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			splitAtDash(calc, Estimators);

			if (abort == false) {
			
				validCalculator = new ValidCalculators();
				
				int i;
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("matrix", Estimators[i]) == true) { 
						if (Estimators[i] == "jabund") { 	
							matrixCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							matrixCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							matrixCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							matrixCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							matrixCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							matrixCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							matrixCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							matrixCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "morisitahorn") { 
							matrixCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							matrixCalculators.push_back(new BrayCurtis());
						}
					}
				}
				
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "MatrixOutputCommand", "MatrixOutputCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void MatrixOutputCommand::help(){
	try {
		mothurOut("The dist.shared command can only be executed after a successful read.otu command.\n");
		mothurOut("The dist.shared command parameters are groups, calc and label.  No parameters are required.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n");
		mothurOut("The group names are separated by dashes. The label parameter allows you to select what distance levels you would like distance matrices created for, and is also separated by dashes.\n");
		mothurOut("The dist.shared command should be in the following format: dist.shared(groups=yourGroups, calc=yourCalcs, label=yourLabels).\n");
		mothurOut("Example dist.shared(groups=A-B-C, calc=jabund-sorabund).\n");
		mothurOut("The default value for groups is all the groups in your groupfile.\n");
		mothurOut("The default value for calc is jclass and thetayc.\n");
		validCalculator->printCalc("matrix", cout);
		mothurOut("The dist.shared command outputs a .dist file for each calculator you specify at each distance you choose.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "MatrixOutputCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

MatrixOutputCommand::~MatrixOutputCommand(){
	if (abort == false) {
		delete input; globaldata->ginput = NULL;
		delete read;
		delete validCalculator;
	}
}

//**********************************************************************************************************************

int MatrixOutputCommand::execute(){
	try {
		
		if (abort == true) {	return 0;	}
			
		//if the users entered no valid calculators don't execute command
		if (matrixCalculators.size() == 0) { mothurOut("No valid calculators."); mothurOutEndLine();  return 0; }

		//you have groups
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
					
		if (lookup.size() < 2) { mothurOut("You have not provided enough valid groups.  I cannot run the command."); mothurOutEndLine(); return 0;}
		
		numGroups = lookup.size();
				
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
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {  delete lookup[i]; }  } 
			lookup = input->getSharedRAbundVectors(lastLabel);

			mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
			process(lookup);
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
				
		//reset groups parameter
		globaldata->Groups.clear();  

		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MatrixOutputCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
void MatrixOutputCommand::printSims(ostream& out) {
	try {
		
		//output column headers
		out << simMatrix.size() << endl;
		
		for (int m = 0; m < simMatrix.size(); m++)	{
			out << lookup[m]->getGroup() << '\t';
			for (int n = 0; n < m; n++)	{
				out << simMatrix[m][n] << '\t'; 
			}
			out << endl;
		}

	}
	catch(exception& e) {
		errorOut(e, "MatrixOutputCommand", "printSims");
		exit(1);
	}
}
/***********************************************************/
void MatrixOutputCommand::process(vector<SharedRAbundVector*> thisLookup){
	try {
	
				EstOutput data;
				vector<SharedRAbundVector*> subset;

				//for each calculator												
				for(int i = 0 ; i < matrixCalculators.size(); i++) {
					
					//initialize simMatrix
					simMatrix.clear();
					simMatrix.resize(numGroups);
					for (int m = 0; m < simMatrix.size(); m++)	{
						for (int j = 0; j < simMatrix.size(); j++)	{
							simMatrix[m].push_back(0.0);
						}
					}
				
					for (int k = 0; k < thisLookup.size(); k++) { 
						for (int l = k; l < thisLookup.size(); l++) {
							if (k != l) { //we dont need to similiarity of a groups to itself
								//get estimated similarity between 2 groups
								
								subset.clear(); //clear out old pair of sharedrabunds
								//add new pair of sharedrabunds
								subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
								
								data = matrixCalculators[i]->getValues(subset); //saves the calculator outputs
								//save values in similarity matrix
								simMatrix[k][l] = 1.0 - data[0];  //convert similiarity to distance
								simMatrix[l][k] = 1.0 - data[0];  //convert similiarity to distance
							}
						}
					}
					
					exportFileName = outputDir + getRootName(getSimpleName(globaldata->inputFileName)) + matrixCalculators[i]->getName() + "." + thisLookup[0]->getLabel() + ".dist";
					openOutputFile(exportFileName, out);
					printSims(out);
					out.close();
					
				}

	
		
	}
	catch(exception& e) {
		errorOut(e, "MatrixOutputCommand", "process");
		exit(1);
	}
}
/***********************************************************/



