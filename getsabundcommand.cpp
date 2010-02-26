/*
 *  getsabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getsabundcommand.h"

//**********************************************************************************************************************

GetSAbundCommand::GetSAbundCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"label","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			string outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}

			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "")) { m->mothurOut("You must read a list or rabund before you can use the get.sabund command."); m->mothurOutEndLine(); abort = true; }
			
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
				
			if (abort == false) {
				filename = outputDir + getRootName(getSimpleName(globaldata->inputFileName)) + "sabund";
				openOutputFile(filename, out);
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "GetSAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetSAbundCommand::help(){
	try {
		m->mothurOut("The get.sabund command can only be executed after a successful read.otu of a listfile or rabundfile.\n");
		m->mothurOut("The get.sabund command parameters is label.  No parameters are required.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like included in your .sabund file, and are separated by dashes.\n");
		m->mothurOut("The get.sabund command should be in the following format: get.sabund(label=yourLabels).\n");
		m->mothurOut("Example get.sabund().\n");
		m->mothurOut("The default value for label is all labels in your inputfile.\n");
		m->mothurOut("The get.sabund command outputs a .sabund file containing the labels you selected.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetSAbundCommand::~GetSAbundCommand(){
}

//**********************************************************************************************************************

int GetSAbundCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
	
		//using order vector so you don't have to distinguish between the list and rabund files
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		order = globaldata->gorder;
		string lastLabel = order->getLabel();
		input = globaldata->ginput;
						
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(order->getLabel()) == 1){
					m->mothurOut(order->getLabel());  m->mothurOutEndLine();
					sabund = new SAbundVector();
					*sabund = (order->getSAbundVector());
					sabund->print(out);
					delete sabund;

					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
			}
			
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = order->getLabel();
					
					delete order;		
					order = (input->getOrderVector(lastLabel));
					
					m->mothurOut(order->getLabel());  m->mothurOutEndLine();
					sabund = new SAbundVector();
					*sabund = (order->getSAbundVector());
					sabund->print(out);
					delete sabund;

					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					//restore real lastlabel to save below
					order->setLabel(saveLabel);
			}
			
			
			lastLabel = order->getLabel();	
			
			delete order;		
			order = (input->getOrderVector());
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + ".");  m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + ".");  m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (order != NULL) {	delete order;	}
			order = (input->getOrderVector(lastLabel));
			
			m->mothurOut(order->getLabel());  m->mothurOutEndLine();
			sabund = new SAbundVector();
			*sabund = (order->getSAbundVector());
			sabund->print(out);
			delete sabund;
			delete order;
		}
		globaldata->gorder = NULL;

		out.close();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(filename); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


