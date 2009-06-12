/*
 *  getrabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "getrabundcommand.h"

//**********************************************************************************************************************

GetRAbundCommand::GetRAbundCommand(string option){
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
			
			parser = new OptionParser();
			parser->parse(option, parameters);  delete parser;
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter->isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getListFile() == "") { cout << "You must read a listfile before you can use the get.rabund command." << endl; abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter->validFile(parameters, "line", false);				
			if (line == "not found") { line = "";  }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			label = validParameter->validFile(parameters, "label", false);			
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
				
			delete validParameter;
			
			if (abort == false) {
				filename = getRootName(globaldata->inputFileName) + "rabund";
				openOutputFile(filename, out);
			}
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRAbundCommand class Function GetRAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRAbundCommand class function GetRAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}
//**********************************************************************************************************************

void GetRAbundCommand::help(){
	try {
		cout << "The get.rabund command can only be executed after a successful read.otu of a listfile." << "\n";
		cout << "The get.rabund command parameters are line and label.  No parameters are required, and you may not use line and label at the same time." << "\n";
		cout << "The line and label allow you to select what distance levels you would like included in your .rabund file, and are separated by dashes." << "\n";
		cout << "The get.rabund command should be in the following format: get.rabund(line=yourLines, label=yourLabels)." << "\n";
		cout << "Example get.rabund(line=1-3-5)." << "\n";
		cout << "The default value for line and label are all lines in your inputfile." << "\n";
		cout << "The get.rabund command outputs a .rabund file containing the lines you selected." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRAbundCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRAbundCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

GetRAbundCommand::~GetRAbundCommand(){
}

//**********************************************************************************************************************

int GetRAbundCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		int count = 1;
		
		//read first line
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		list = globaldata->gListVector;
		ListVector* lastList = list;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;

		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(allLines == 1 || lines.count(count) == 1 || labels.count(list->getLabel()) == 1){
					cout << list->getLabel() << '\t' << count << endl;
					rabund = new RAbundVector();
					*rabund = (list->getRAbundVector());
					rabund->print(out);
					delete rabund;
															
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
					cout << lastList->getLabel() << '\t' << count << endl;
					rabund = new RAbundVector();
					*rabund = (lastList->getRAbundVector());
					rabund->print(out);
					delete rabund;

					processedLabels.insert(lastList->getLabel());
					userLabels.erase(lastList->getLabel());
			}
			
			if (count != 1) { delete lastList; }
			lastList = list;			
			
			list = input->getListVector();
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastList->getLabel()) != 1) {
				cout << ". I will use " << lastList->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastList->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastList->getLabel() << '\t' << count << endl;
			rabund = new RAbundVector();
			*rabund = (lastList->getRAbundVector());
			rabund->print(out);
			delete rabund;
		}
		delete lastList;

		out.close(); 
		return 0;		
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRAbundCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRAbundCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************


