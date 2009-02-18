/*
 *  errorchecking.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "errorchecking.h"
#include <math.h>

/*******************************************************/

/******************************************************/

ErrorCheck::ErrorCheck() {
	globaldata = GlobalData::getInstance();
	validCommand = new ValidCommands();
	validParameter = new ValidParameters();
}
/*******************************************************/

/******************************************************/

void ErrorCheck::refresh() {
	columnfile = globaldata->getColumnFile();
	phylipfile = globaldata->getPhylipFile();
	listfile = globaldata->getListFile();
	rabundfile = globaldata->getRabundFile();
	sabundfile = globaldata->getSabundFile();
	namefile = globaldata->getNameFile();
	groupfile = globaldata->getGroupFile();
	orderfile = globaldata->getOrderFile();
	fastafile = globaldata->getFastaFile();
	treefile = globaldata->getTreeFile();
	cutoff = globaldata->getCutOff();
	format = globaldata->getFormat();
	method = globaldata->getMethod();

	
	string p[] = {
		"phylip",              //0
		"column",             //1
		"list",               //2
		"rabund",             //3
		"sabund",             //4
		"name",               //5
		"order",              //6
		"group",              //7
		"fasta",              //8
		"treefile",           //9
		"cutoff",             //10
		"precision",          //11
		"iters",              //12
		"jumble",             //13
		"freq",               //14
		"method",             //15
		"fileroot",           //16
		"line",               //17
		"label",              //18
		"single",             //19
		"rarefaction",        //20
		"shared",             //21
		"summary",            //22
		"sharedrarefaction",  //23
		"sharedsummary",      //24
		"comparegroups",      //25
		"abund",              //26
		};
	
	string c0[] = {p[0],p[5],p[10],p[11]};
	string c1[] = {p[2],p[6],p[7]}; 
	string c2[] = {p[10],p[11],p[15]}; 
	string c3[] = {p[8]};  
	string c4[] = {p[14],p[17],p[18],p[19],p[26]};
	string c5[] = {p[13],p[14],p[17],p[18],p[21],p[25]};
	string c6[] = {""}; 
	string c7[] = {""}; 
	string c8[] = {""}; 
	string c9[] = {p[12],p[14],p[17],p[18],p[20],p[26]};
	string c10[] = {p[12],p[13],p[17],p[18],p[23]};
	string c11[] = {p[17],p[18],p[22],p[26]};   
	string c12[] =  {p[13],p[17],p[18],p[24]}; 
	string c13[] = {""}; 	
	
	vector<string> v0 (c0, c0+sizeof(c0)/sizeof(string)); 
	vector<string> v1 (c1, c1+sizeof(c1)/sizeof(string));
	vector<string> v2 (c2, c2+sizeof(c2)/sizeof(string));
	vector<string> v3 (c3, c3+sizeof(c3)/sizeof(string));
	vector<string> v4 (c4, c4+sizeof(c4)/sizeof(string));
	vector<string> v5 (c5, c5+sizeof(c5)/sizeof(string));
	vector<string> v6 (c6, c6+sizeof(c6)/sizeof(string));
	vector<string> v7 (c7, c7+sizeof(c7)/sizeof(string));
	vector<string> v8 (c8, c8+sizeof(c8)/sizeof(string));
	vector<string> v9 (c9, c9+sizeof(c9)/sizeof(string));
	vector<string> v10 (c10, c10+sizeof(c10)/sizeof(string));
	vector<string> v11 (c11, c11+sizeof(c11)/sizeof(string));
	vector<string> v12 (c12, c12+sizeof(c12)/sizeof(string));
	vector<string> v13 (c13, c13+sizeof(c13)/sizeof(string));
	
	vector<vector<string> > allCommands;
	allCommands.push_back(v0);
	allCommands.push_back(v1);
	allCommands.push_back(v2);
	allCommands.push_back(v3);
	allCommands.push_back(v4);
	allCommands.push_back(v5);
	allCommands.push_back(v6);
	allCommands.push_back(v7);
	allCommands.push_back(v8);
	allCommands.push_back(v9);
	allCommands.push_back(v10);
	allCommands.push_back(v11);
	allCommands.push_back(v12);
	allCommands.push_back(v13);
	
	string commands[] = {
	"read.dist",          //0
	"read.otu",           //1
	"cluster",            //2
	"deconvolute",        //3
	"collect.single",     //4
	"collect.shared",     //5
	"get.group",          //6
	"get.label",          //7
	"get.line",           //8
	"rarefaction.single", //9
	"rarefaction.shared", //10
	"summary.single",     //11
	"summary.shared",     //12
	"quit"                //13
	};
	
	for(int i = 0; i < allCommands.size(); i++)
		commandParameters[commands[i]] = allCommands.at(i);
	
	 //{Lowerbound(piSent if no lowerbound), Upperbound(piSent if no upperbound), 1 if only the first 2 values, 0 if greater than, 0 if less than};
	piSent = 3.14159;
	double ip0[] = {10, piSent, 0, 1, 0};
	double ip1[] = {10, piSent, 0, 1, 0};
	double ip2[] = {0, 1, 1, 0, 0};
	double ip3[] =  {1, piSent, 0, 0, 0};
	double ip4[] = {1, piSent, 0, 1, 0};
	double ip5[] = {5, piSent, 0, 1, 0};
	
	vector<double> ipv0 (ip0, ip0+sizeof(ip0)/sizeof(double)); 
	vector<double> ipv1 (ip1, ip1+sizeof(ip1)/sizeof(double)); 
	vector<double> ipv2 (ip2, ip2+sizeof(ip2)/sizeof(double)); 
	vector<double> ipv3 (ip3, ip3+sizeof(ip3)/sizeof(double)); 
	vector<double> ipv4 (ip4, ip4+sizeof(ip4)/sizeof(double)); 
	vector<double> ipv5 (ip5, ip5+sizeof(ip5)/sizeof(double));

	intParams[p[11]] = ipv0;
	intParams[p[12]] = ipv1;
	intParams[p[13]] = ipv2;
	intParams[p[14]] = ipv3;
	intParams[p[17]] = ipv4;
	intParams[p[26]] = ipv5;
	
	randomtree = globaldata->getRandomTree();
	sharedfile = globaldata->getSharedFile();
}

/*******************************************************/

/******************************************************/

ErrorCheck::~ErrorCheck() {}

/*******************************************************/

/******************************************************/

bool ErrorCheck::checkInput(string input) {
		errorFree = true;
		clear();
		
		//refresh variable
		refresh();
		
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
				splitAtComma(value, optionText);
				splitAtEquals(parameter, value);
				
				//is it a valid parameter
				if (validParameter->isValidParameter(parameter) != true) { return false; }
				if(!validCommandParameter(parameter,commandName)) { 
					cout << "'" << parameter << "' is not a valid parameter for the " << commandName << " command.\n";
					return false; 
				}
				if(!validParameterValue(value, parameter)) {
					if(parameter.compare("precision") == 0)
						cout << "The precision parameter can only take powers of 10 as a value (e.g. 10,1000,1000, etc.)\n";
					else {
					vector<double> bounds = intParams[parameter];
					double a = bounds.at(0);
					double b = bounds.at(1);
					double c = bounds.at(2);
					double d = bounds.at(3);
					double e = bounds.at(4);
					cout << "The '" << parameter << "' parameter needs to be ";
					if(c == 1)
							cout << "either '" << a << "' or '" << b << "'.\n";
					else
					{
						if(a != piSent)
						{
							cout << ">";
							if(d != 0)
								cout << "=";
							cout << " '" << a << "'";
						}
						if(b == piSent)
							cout << ".\n";
						else if(a != piSent)
							cout << " and ";
						if(b != piSent)
						{
							cout << "<";
							if(e != 0)
								cout << "=";
							cout << " '" << b << ".\n";
						}
					}
					}
					return false;
				}

				if (parameter == "phylip" )		{ phylipfile = value; }
				if (parameter == "column" )		{ columnfile = value; }
				if (parameter == "list" )		{ listfile = value; }
				if (parameter == "rabund" )		{ rabundfile = value; }
				if (parameter == "sabund" )		{ sabundfile = value; }
				if (parameter == "name" )		{ namefile = value; }
				if (parameter == "order" )		{ orderfile = value; }
				if (parameter == "fasta" )		{ fastafile = value; }
				if (parameter == "tree" )		{ treefile = value; }
				if (parameter == "group" )		{ groupfile = value; }
				if (parameter == "shared" )		{ sharedfile = value; }
				if (parameter == "cutoff" )			{ cutoff = value; }
				if (parameter == "precision" )		{ precision = value; }
				if (parameter == "iters" )			{ iters = value; }
				if (parameter == "jumble" )			{ jumble = value; }
				if (parameter == "freq" )			{ freq = value; }
				if (parameter == "method" )			{ method = value; }
				if (parameter == "fileroot" )		{ fileroot = value; }
				if (parameter == "line" )			{ line = value; }
				if (parameter == "label" )			{ label = value; }
				if (parameter == "abund" )          { abund = value; }
				if (parameter == "random" )			{ randomtree = value;	}
			}
			
			//gets the last parameter and value
			if (errorFree)  { //gets the last parameter and value
				value = optionText;
				splitAtEquals(parameter, value);
				//is it a valid parameter
				if (validParameter->isValidParameter(parameter) != true) { return false; }
				if(!validCommandParameter(parameter,commandName)) { 
					cout << "'" << parameter << "' is not a valid parameter for the " << commandName << " command.\n";
					return false; 
				}
				if(!validParameterValue(value, parameter)) {
					if(parameter.compare("precision") == 0)
						cout << "The precision parameter can only take powers of 10 as a value (e.g. 10,1000,1000, etc.)\n";
					else {
					vector<double> bounds = intParams[parameter];
					double a = bounds.at(0);
					double b = bounds.at(1);
					double c = bounds.at(2);
					double d = bounds.at(3);
					double e = bounds.at(4);
					cout << "The '" << parameter << "' parameter needs to be ";
					if(c == 1)
							cout << "either '" << a << "' or '" << b << "'.\n";
					else
					{
						if(a != piSent)
						{
							cout << ">";
							if(d != 0)
								cout << "=";
							cout << " '" << a << "'";
						}
						if(b == piSent)
							cout << ".\n";
						else if(a != piSent)
							cout << " and ";
						if(b != piSent)
						{
							cout << "<";
							if(e != 0)
								cout << "=";
							cout << " '" << b << ".\n";
						}
					}
					}
					return false;
				}
				if (parameter == "phylip" )		{ phylipfile = value; }
				if (parameter == "column" )		{ columnfile = value; }				
				if (parameter == "list" )		{ listfile = value; }
				if (parameter == "rabund" )		{ rabundfile = value; }
				if (parameter == "sabund" )		{ sabundfile = value; }
				if (parameter == "name" )		{ namefile = value; }
				if (parameter == "order" )		{ orderfile = value; }
				if (parameter == "group" )		{ groupfile = value; }
				if (parameter == "shared" )		{ sharedfile = value; }
				if (parameter == "fasta" )		{ fastafile = value; }
				if (parameter == "tree" )		{ treefile = value; }
				if (parameter == "cutoff" )			{ cutoff = value; }
				if (parameter == "precision" )		{ precision = value; }
				if (parameter == "iters" )			{ iters = value; }
				if (parameter == "jumble" )			{ jumble = value; }
				if (parameter == "freq" )			{ freq = value; }
				if (parameter == "method" )			{ method = value; }
				if (parameter == "fileroot" )		{ fileroot = value; }
				if (parameter == "line" )			{ line = value; }
				if (parameter == "label" )			{ label = value; }
				if (parameter == "random" )			{ randomtree = value;	}
				if (parameter == "abund" )          { abund = value; }
			}
		}
		
		//make sure the user does not use both the line and label parameters
		if ((line != "") && (label != "")) { cout << "You may use either the line or label parameters, but not both." << endl; return false; }
		
		if (commandName == "read.dist") { 
			validateReadFiles();
			validateReadDist();
		}else if (commandName == "read.otu") { 
			//you want to do shared commands
			if ((listfile != "") && (groupfile != ""))	{
				validateParseFiles(); //checks the listfile and groupfile parameters
			}else { //you want to do single commands
				validateReadFiles();
				validateReadPhil();
			}
		}else if (commandName == "read.shared") { 
			//you want to do shared commands with just the shared file
			validateReadFiles();
		}else if (commandName == "read.tree") { 
			validateTreeFiles(); //checks the treefile and groupfile parameters
		}else if (commandName == "deconvolute") {
			if (fastafile == "") { cout << "You must enter a fastafile with the deconvolute() command." << endl; return false; }
			validateReadFiles();
		}
		
		//are you trying to cluster before you have read something			
		if ((commandName == "cluster") && (globaldata->getSparseMatrix() == NULL) ||
			(commandName == "cluster") && (globaldata->getListVector() == NULL)) {
				cout << "Before you use the cluster command, you first need to read in a distance matrix." << endl; 
				errorFree = false;
		} 
		
		if (commandName == "parsimony") {
			//are you trying to use parsimony without reading a tree or saying you want random distribution
			if (randomtree == "")  {
				if (globaldata->gTree.size() == 0) {
					cout << "You must read a treefile and a groupfile or set the randomtree parameter to the output filename you wish, before you may execute the parsimony command." << endl; return false;  }
			}
		}
		
		if ((commandName == "unifrac.weighted") || (commandName == "unifrac.unweighted")) {
			if (globaldata->gTree.size() == 0) {//no trees were read
				cout << "You must execute the read.tree command, before you may execute the unifrac.weighted or unifrac.unweighted command." << endl; return false;  }
		}
		
		//check for valid method
		if(commandName == "get.group") {
			if ((globaldata->getGroupFile() == "")) { cout << "You must read a group before you can use the get.group command." << endl; return false; }
		}
		if (commandName == "get.label" || commandName == "get.line") {
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { cout << "You must read a list, sabund or rabund before you can use the get.label or get.line command." << endl; return false; }
		}
		if (commandName == "cluster") {
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else {cout << "Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average." << endl; return false; }
		}
		
		if ((commandName == "collect.single") || (commandName == "rarefaction.single") || (commandName == "summary.single") ){ 
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { cout << "You must read a list, sabund or rabund before you can use the collect.single, rarefaction.single or summary.single commands." << endl; return false; }
		}
		
		if ((commandName == "collect.shared") || (commandName == "rarefaction.shared") || (commandName == "summary.shared") ){ 
			if (globaldata->getSharedFile() == "") {
				if (globaldata->getListFile() == "") { cout << "You must read a list and a group, or a shared before you can use the collect.shared, rarefaction.shared or summary.shared commands." << endl; return false; }
				else if (globaldata->getGroupFile() == "") { cout << "You must read a list and a group, or a shared before you can use the collect.shared, rarefaction.shared or summary.shared commands." << endl; return false; }
			}
		}

		globaldata->clearAbund();

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
		}else if (fastafile != "") {
			ableToOpen = openInputFile(fastafile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
			else { globaldata->inputFileName = fastafile; }
		}else if (sharedfile != "") {
			ableToOpen = openInputFile(sharedfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
			else { globaldata->inputFileName = sharedfile; }
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
//This function checks to see if the given paramter
//is a valid paramter for the given command.
bool ErrorCheck::validCommandParameter(string parameter, string commandName) {
	try {
		for(int i = 0; i < commandParameters[commandName].size(); i++)
			if(parameter.compare(commandParameters[commandName][i]) == 0)
				return true;
		return false;
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
//This function checks to see if the given paramter value
//is convertable into an int if that parameter requires it.
bool ErrorCheck::validParameterValue(string value, string parameter) {
	try {
		int pVal;
		if(intParams.count(parameter) == 1)
		{
			vector<double> bounds = intParams[parameter];
			bool valid = convertTest(value, pVal);
			if(!valid)
				return false;
			if(parameter.compare("precision") == 0)
			{
				double logNum = log10((double)pVal);
				double diff = (double)((int)logNum - logNum);
				if(diff != 0)
					return false;
			}
			double a = bounds.at(0);
			double b = bounds.at(1);
			double c = bounds.at(2);
			double d = bounds.at(3);
			double e = bounds.at(4);
			bool a0 = pVal > a;
			bool a1 = pVal >= a;
			bool b0 = pVal < b;
			bool b1 = pVal <= b;
			
			if(c != 1)
			{
				if(a == piSent && b == piSent)
					return true;
				if(a != piSent && b == piSent)
				{
					if(d == 0)
						return a0;
					else
						return a1;
				}
				else if(a == piSent && b != piSent)
				{
					if(e == 0)
						return b0;
					else
						return b1;
				}
				else
				{
					if(d == 0 && e == 0)
						return (a0 && b0);
					else if(d == 0 && e == 1)
						return (a0 && b1);
					else if(d == 1 && e == 0)
						return (a1 && b0);
					else
						return (a1 && b1);
				}
			}
			else
			{
				if(a == piSent && b == piSent)
					return true;
				if(a != piSent && b == piSent)
					return (pVal == a);
				else if(a == piSent && b != piSent)
					return (pVal == b);
				else
					return (pVal == a || pVal == b);
			}
		}
		return true;
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
		
		if ((phylipfile == "") && (columnfile == "")) { cout << "When executing a read.dist you must enter a phylip or a column." << endl; errorFree = false; }
		else if ((phylipfile != "") && (columnfile != "")) { cout << "When executing a read.dist you must enter ONLY ONE of the following: phylip or column." << endl; errorFree = false; }
		
		if (columnfile != "") {
			if (namefile == "") {
				cout << "You need to provide a namefile if you are going to use the column format." << endl;
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
	
		if (listfile == "") { cout << "When executing a read.otu for groups you must enter a list and a group." << endl; errorFree = false; }
		else if (groupfile == "") { cout << "When executing a read.otu for groups you must enter a list and a group." << endl; errorFree = false; }
	
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
// format parameters on a parselistcommand
void ErrorCheck::validateTreeFiles() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		//checks for valid files
	
		if (treefile == "") { cout << "When executing a read.tree you must enter a treefile and a groupfile." << endl; errorFree = false; }
		else if (groupfile == "") { cout << "When executing a read.tree you must enter a treefile and a groupfile." << endl; errorFree = false; }
	
		//checks parameters on the read command
		if (treefile != "") {
			ableToOpen = openInputFile(treefile, filehandle);
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
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateTreeFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateTreeFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
				cout << "When executing a read.otu you must enter ONLY ONE of the following: list, rabund or sabund." << endl; errorFree = false; }
		}else if (rabundfile != "") { 
			if ((listfile != "") || (sabundfile != "")) { 
				cout << "When executing a read.otu you must enter ONLY ONE of the following: list, rabund or sabund." << endl; errorFree = false; }
		}else if (sabundfile != "") { 
			if ((listfile != "") || (rabundfile != "")) { 
				cout << "When executing a read.otu you must enter ONLY ONE of the following: list, rabund or sabund." << endl; errorFree = false; }
		}else if ((listfile == "") && (rabundfile == "") && (sabundfile == "")) {
			    cout << "When executing a read.otu you must enter one of the following: list, rabund or sabund." << endl; errorFree = false; 
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
	sharedfile		=	"";
	line			=	"";
	label			=	"";
	method			=   "furthest";
}
/*******************************************************/

/******************************************************/

