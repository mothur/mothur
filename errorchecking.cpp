/*
 *  errorchecking.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "errorchecking.h"

/*******************************************************/

/******************************************************/

ErrorCheck::ErrorCheck() {
	globaldata = GlobalData::getInstance();
	validCommand = new ValidCommands();
	validParameter = new ValidParameters();
	validCalculator = new ValidCalculators();
	columnfile = globaldata->getColumnFile();
	phylipfile = globaldata->getPhylipFile();
	listfile = globaldata->getListFile();
	rabundfile = globaldata->getRabundFile();
	sabundfile = globaldata->getSabundFile();
	namefile = globaldata->getNameFile();
	groupfile = globaldata->getGroupFile();
	orderfile = globaldata->getOrderFile();
	cutoff = globaldata->getCutOff();
	format = globaldata->getFormat();
	method = globaldata->getMethod();

}
/*******************************************************/

/******************************************************/

ErrorCheck::~ErrorCheck() {}

/*******************************************************/

/******************************************************/

bool ErrorCheck::checkInput(string input) {
		errorFree = true;
		clear();

		//get command name and parameters
		int openParen = input.find_first_of('(');
		int closeParen = input.find_last_of(')');

		if(openParen != -1 && closeParen != -1){			
			commandName = input.substr(0, openParen);   //commandName contains everything before "("
			optionText = input.substr(openParen+1, closeParen-openParen-1); //optionString contains everything between "(" and ")".
		}else if (openParen == -1) { //there is no parenthesis
			cout << input << " is not a valid command. You are missing the ()." << endl;
			return false;
		}
		
		//is it a valid command
		if (validCommand->isValidCommand(commandName) != true) { return false; }
		
		string parameter, value;
		//reads in parameters and values
		if((optionText != "") && (commandName != "help")){
			while((optionText.find_first_of(',') != -1) && (errorFree)) {  //while there are parameters
				globaldata->splitAtComma(value, optionText);
				globaldata->splitAtEquals(parameter, value);
				
				//is it a valid parameter
				if (validParameter->isValidParameter(parameter) != true) { return false; }
				
				if (parameter == "phylipfile" )		{ phylipfile = value; }
				if (parameter == "columnfile" )		{ columnfile = value; }
				if (parameter == "listfile" )		{ listfile = value; }
				if (parameter == "rabundfile" )		{ rabundfile = value; }
				if (parameter == "sabundfile" )		{ sabundfile = value; }
				if (parameter == "namefile" )		{ namefile = value; }
				if (parameter == "orderfile" )		{ orderfile = value; }
				if (parameter == "groupfile" )		{ groupfile = value; }
				if (parameter == "cutoff" )			{ cutoff = value; }
				if (parameter == "precision" )		{ precision = value; }
				if (parameter == "iters" )			{ iters = value; }
				if (parameter == "jumble" )			{ jumble = value; }
				if (parameter == "freq" )			{ freq = value; }
				if (parameter == "method" )			{ method = value; }
				if (parameter == "fileroot" )		{ fileroot = value; }
				if (parameter == "line" )			{ line = value; }
				if (parameter == "label" )			{ label = value; }

				if (parameter == "single") {//stores estimators in a vector
					singleEsimators.clear(); //clears out old values
					globaldata->splitAtDash(value, singleEsimators);
					for (int i = 0; i < singleEsimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, singleEsimators[i]) != true) { return false; }
					}
				}
				if (parameter == "rarefaction") {//stores estimators in a vector
					rareEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, rareEstimators);
					for (int i = 0; i < rareEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, rareEstimators[i]) != true) { return false; }
					}
				}
				if (parameter == "shared") {//stores estimators in a vector
					sharedEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, sharedEstimators);
					for (int i = 0; i < sharedEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, sharedEstimators[i]) != true) { return false; }
					}
				}
				if (parameter == "summary") { //stores summaries to be used in a vector
					summaryEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, summaryEstimators);
					for (int i = 0; i < summaryEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, summaryEstimators[i]) != true) { return false; }
					}
				}
				if (parameter == "sharedrarefaction") { //stores summaries to be used in a vector
					sharedRareEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, sharedRareEstimators);
					for (int i = 0; i < sharedRareEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, sharedRareEstimators[i]) != true) { return false; }
					}
				}
			}
			
			//gets the last parameter and value
			if (errorFree)  { //gets the last parameter and value
				value = optionText;
				globaldata->splitAtEquals(parameter, value);
				//is it a valid parameter
				if (validParameter->isValidParameter(parameter) != true) { return false; }
				
				if (parameter == "phylipfile" )		{ phylipfile = value; }
				if (parameter == "columnfile" )		{ columnfile = value; }				
				if (parameter == "listfile" )		{ listfile = value; }
				if (parameter == "rabundfile" )		{ rabundfile = value; }
				if (parameter == "sabundfile" )		{ sabundfile = value; }
				if (parameter == "namefile" )		{ namefile = value; }
				if (parameter == "orderfile" )		{ orderfile = value; }
				if (parameter == "groupfile" )		{ groupfile = value; }
				if (parameter == "cutoff" )			{ cutoff = value; }
				if (parameter == "precision" )		{ precision = value; }
				if (parameter == "iters" )			{ iters = value; }
				if (parameter == "jumble" )			{ jumble = value; }
				if (parameter == "freq" )			{ freq = value; }
				if (parameter == "method" )			{ method = value; }
				if (parameter == "fileroot" )		{ fileroot = value; }
				if (parameter == "line" )			{ line = value; }
				if (parameter == "label" )			{ label = value; }

				if (parameter == "single") {//stores estimators in a vector
					singleEsimators.clear(); //clears out old values
					globaldata->splitAtDash(value, singleEsimators);
					for (int i = 0; i < singleEsimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, singleEsimators[i]) != true) { return false; }
					}
				}
				if (parameter == "rarefaction") {//stores estimators in a vector
					rareEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, rareEstimators);
					for (int i = 0; i < rareEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, rareEstimators[i]) != true) { return false; }
					}
				}
				if (parameter == "shared") {//stores estimators in a vector
					sharedEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, sharedEstimators);
					for (int i = 0; i < sharedEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, sharedEstimators[i]) != true) { return false; }
					}
				}
				if (parameter == "summary") { //stores summaries to be used in a vector
					summaryEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, summaryEstimators);
					for (int i = 0; i < summaryEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, summaryEstimators[i]) != true) { return false; }
					}
				}
				if (parameter == "sharedrarefaction") { //stores summaries to be used in a vector
					sharedRareEstimators.clear(); //clears out old values
					globaldata->splitAtDash(value, sharedRareEstimators);
					for (int i = 0; i < sharedRareEstimators.size(); i++) { //loop through estimators
						//is it a valid calculator
						if (validCalculator->isValidCalculator(parameter, sharedRareEstimators[i]) != true) { return false; }
					}
				}

			}
		}
		
		//make sure the user does not use both the line and label parameters
		if ((line != "") && (label != "")) { cout << "You may use either the line or label parameters, but not both." << endl; return false; }
	
		
		if (commandName == "read.dist") { 
			validateReadFiles();
			validateReadDist();
		}else if (commandName == "read.otu") { 
			validateReadFiles();
			validateReadPhil();	
		}else if (commandName == "read.list") { 
			validateParseFiles(); //checks the listfile and groupfile parameters
		}
		
		//are you trying to cluster before you have read something			
		if ((commandName == "cluster") && (globaldata->getSparseMatrix() == NULL) ||
			(commandName == "cluster") && (globaldata->getListVector() == NULL)) {
				cout << "Before you use the cluster command, you first need to read in a distance matrix." << endl; 
				errorFree = false;
		} 
		
		//check for valid method
		if (commandName == "cluster") {
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else {cout << "Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average." << endl; return false; }
		}
		
		if ((commandName == "collect.single") || (commandName == "rarefaction.single") || (commandName == "summary.single") ){ 
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { cout << "You must read a listfile, sabundfile or rabundfile before you can use the collect.single, rarefaction.single or summary.single commands." << endl; return false; }
		}
		
		if ((commandName == "collect.shared") || (commandName == "rarefaction.shared") || (commandName == "summary.shared") || (commandName == "shared") ){ 
			if (globaldata->getListFile() == "") { cout << "You must read a listfile and a groupfile before you can use the collect.shared, rarefaction.shared, summary.shared or shared commands." << endl; return false; }
			else if (globaldata->getGroupFile() == "") { cout << "You must read a listfile and a groupfile before you can use the collect.shared, rarefaction.shared, summary.shared or shared commands." << endl; return false; }
		}
 
		
		return errorFree;
}

/*******************************************************/

/******************************************************/
//This function checks to make sure the user entered a file to 
// read and that the file exists and can be opened.
void ErrorCheck::validateReadFiles() {
	try {
		//Validating files for read
		ifstream filehandle;
		int ableToOpen;
	
		//are we reading a phylipfile
		if (phylipfile != "") {
			ableToOpen = openInputFile(phylipfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) { errorFree = false; }
			else { globaldata->inputFileName = phylipfile; }
		//are we reading a phylipfile
		}else if (columnfile != "") {
			ableToOpen = openInputFile(columnfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) { errorFree = false; }
			else { globaldata->inputFileName = columnfile; }
		//are we reading a listfile
		}else if (listfile!= "") {
			ableToOpen = openInputFile(listfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
			else { globaldata->inputFileName = listfile; }
		//are we reading a rabundfile
		}else if (rabundfile != "") {
			ableToOpen = openInputFile(rabundfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
			else { globaldata->inputFileName = rabundfile; }
		//are we reading a sabundfile
		}else if (sabundfile != "") {
			ableToOpen = openInputFile(sabundfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
			else { globaldata->inputFileName = sabundfile; }
		}else{ //no file given
			errorFree = false;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateReadFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateReadFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}
/*******************************************************/

/******************************************************/
//This function checks to make sure the user entered appropriate
// format parameters on a distfile read
void ErrorCheck::validateReadDist() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		if ((phylipfile == "") && (columnfile == "")) { cout << "When executing a read.dist you must enter a phylipfile or a columnfile." << endl; errorFree = false; }
		else if ((phylipfile != "") && (columnfile != "")) { cout << "When executing a read.dist you must enter ONLY ONE of the following: phylipfile or columnfile." << endl; errorFree = false; }
		
		if (columnfile != "") {
			if (namefile == "") {
				cout << "You need to provide a namefile name if you are going to use the column format." << endl;
				errorFree = false; 
			}else {
				ableToOpen = openInputFile(namefile, filehandle);
				filehandle.close();
				//unable to open
				if (ableToOpen == 1) { errorFree = false; }
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateReadDist. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateReadDist. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************/

/******************************************************/
//This function checks to make sure the user entered appropriate
// format parameters on a parselistcommand
void ErrorCheck::validateParseFiles() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		//checks for valid files
	
		if (listfile == "") { cout << "When executing a read.list you must enter a listfile and a groupfile." << endl; errorFree = false; }
		else if (groupfile == "") { cout << "When executing a read.list you must enter a listfile and a groupfile." << endl; errorFree = false; }
	
		//checks parameters on the read command
		if (listfile != "") {
			ableToOpen = openInputFile(listfile, filehandle);
			filehandle.close();
			if (ableToOpen == 1) { //unable to open
				errorFree = false;
			}
			if (groupfile != "") {
				ableToOpen = openInputFile(groupfile, filehandle);
				filehandle.close();
				if (ableToOpen == 1) { //unable to open
					errorFree = false;;
				}
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateReadPhil. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateReadPhil. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*******************************************************/

/******************************************************/
//This function checks to make sure the user entered appropriate
// format parameters on a distfile read
void ErrorCheck::validateReadPhil() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		//checks to make sure only one file type is given
		if (listfile != "") { 
			if ((rabundfile != "") || (sabundfile != "")) { 
				cout << "When executing a read.otu you must enter ONLY ONE of the following: listfile, rabundfile or sabundfile." << endl; errorFree = false; }
		}else if (rabundfile != "") { 
			if ((listfile != "") || (sabundfile != "")) { 
				cout << "When executing a read.otu you must enter ONLY ONE of the following: listfile, rabundfile or sabundfile." << endl; errorFree = false; }
		}else if (sabundfile != "") { 
			if ((listfile != "") || (rabundfile != "")) { 
				cout << "When executing a read.otu you must enter ONLY ONE of the following: listfile, rabundfile or sabundfile." << endl; errorFree = false; }
		}else if ((listfile == "") && (rabundfile == "") && (sabundfile == "")) {
			    cout << "When executing a read.otu you must enter one of the following: listfile, rabundfile or sabundfile." << endl; errorFree = false; 
		}
		
		//checks parameters on the read command
		if (orderfile != "") {
			ableToOpen = openInputFile(orderfile, filehandle);
			filehandle.close();
			if (ableToOpen == 1) { //unable to open
				errorFree = false;
			}
		}	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateReadPhil. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateReadPhil. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************/

/******************************************************/

void ErrorCheck::clear() {
	//option definitions should go here...
	phylipfile		=	"";
	columnfile		=	"";
	listfile		=	"";
	rabundfile		=	"";
	sabundfile		=	"";
	namefile		=	"";
	groupfile		=	""; 
	orderfile		=	"";
	line			=	"";
	label			=	"";
	method			=   "furthest";
}
/*******************************************************/

/******************************************************/

