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

GetSAbundCommand::GetSAbundCommand(string option){
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
			string Array[] =  {"line","label"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "")) { cout << "You must read a list or rabund before you can use the get.sabund command." << endl; abort = true; }
			
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
			else if((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
				
			if (abort == false) {
				filename = getRootName(globaldata->inputFileName) + "sabund";
				openOutputFile(filename, out);
			}
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetSAbundCommand class Function GetSAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetSAbundCommand class function GetSAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}
//**********************************************************************************************************************

void GetSAbundCommand::help(){
	try {
		cout << "The get.sabund command can only be executed after a successful read.otu of a listfile." << "\n";
		cout << "The get.sabund command parameters are line and label.  No parameters are required, and you may not use line and label at the same time." << "\n";
		cout << "The line and label allow you to select what distance levels you would like included in your .sabund file, and are separated by dashes." << "\n";
		cout << "The get.sabund command should be in the following format: get.sabund(line=yourLines, label=yourLabels)." << "\n";
		cout << "Example get.sabund(line=1-3-5)." << "\n";
		cout << "The default value for line and label are all lines in your inputfile." << "\n";
		cout << "The get.sabund command outputs a .sabund file containing the lines you selected." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetSAbundCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetSAbundCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
	
		int count = 1;
		
		//using order vector so you don't have to distinguish between the list and rabund files
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		order = globaldata->gorder;
		lastOrder = order;
		input = globaldata->ginput;
						
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;

		
		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(allLines == 1 || lines.count(count) == 1 || labels.count(order->getLabel()) == 1){
					cout << order->getLabel() << '\t' << count << endl;
					sabund = new SAbundVector();
					*sabund = (order->getSAbundVector());
					sabund->print(out);
					delete sabund;

					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastOrder->getLabel()) != 1)) {
					cout << lastOrder->getLabel() << '\t' << count << endl;
					sabund = new SAbundVector();
					*sabund = (lastOrder->getSAbundVector());
					sabund->print(out);
					delete sabund;

					processedLabels.insert(lastOrder->getLabel());
					userLabels.erase(lastOrder->getLabel());
			}
			
			if (count != 1) { delete lastOrder; }
			lastOrder = order;	
					
			order = (input->getOrderVector());
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastOrder->getLabel()) != 1) {
				cout << ". I will use " << lastOrder->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastOrder->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastOrder->getLabel() << '\t' << count << endl;
			sabund = new SAbundVector();
			*sabund = (lastOrder->getSAbundVector());
			sabund->print(out);
			delete sabund;
		}
		delete lastOrder;  globaldata->gorder = NULL;

		out.close();
		return 0;		
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetSAbundCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetSAbundCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************


