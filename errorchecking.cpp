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
}
/*******************************************************/

/******************************************************/

void ErrorCheck::refresh() {

	//columnfile = globaldata->getColumnFile();
	//phylipfile = globaldata->getPhylipFile();
	//listfile = globaldata->getListFile();
	//rabundfile = globaldata->getRabundFile();
	//sabundfile = globaldata->getSabundFile();
	//namefile = globaldata->getNameFile();
	//groupfile = globaldata->getGroupFile();
	//orderfile = globaldata->getOrderFile();
	//fastafile = globaldata->getFastaFile();
	//treefile = globaldata->getTreeFile();
	//cutoff = globaldata->getCutOff();
	//format = globaldata->getFormat();
	//method = globaldata->getMethod();
	//randomtree = globaldata->getRandomTree();
	//sharedfile = globaldata->getSharedFile();

}

/*******************************************************/

/******************************************************/

ErrorCheck::~ErrorCheck() {
	delete validCommand;
	delete validParameter;
}

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
				if (validParameter->isValidParameter(parameter, commandName, value) != true) { return false; }
				
				if (parameter == "phylip" )		{ phylipfile = value; }
				if (parameter == "column" )		{ columnfile = value; }
				if (parameter == "list" )		{ listfile = value; }
				if (parameter == "rabund" )		{ rabundfile = value; }
				if (parameter == "sabund" )		{ sabundfile = value; }
				if (parameter == "name" )		{ namefile = value; }
				if (parameter == "order" )		{ orderfile = value; }
				if (parameter == "fasta" )		{ fastafile = value; }
				if (parameter == "nexus" )		{ nexusfile = value; }
				if (parameter == "clustal" )	{ clustalfile = value; }
				if (parameter == "tree" )		{ treefile = value; }
				if (parameter == "group" )			{ groupfile = value; }
				if (parameter == "shared" )			{ sharedfile = value; }
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
				if (parameter == "random" )			{ randomtree = value; }
				if (parameter == "sorted" )			{ sorted = value; }
				if (parameter == "trump" )          { trump = value; }
				if (parameter == "soft" )			{ soft = value; }
				if (parameter == "filter" )         { filter = value; }
				if (parameter == "scale" )			{ scale = value;	}
				if (parameter == "ends" )			{ ends = value; }
				if (parameter == "processors" )		{ processors = value;	}

			}
			
			//gets the last parameter and value
			if (errorFree)  { //gets the last parameter and value
				value = optionText;
				splitAtEquals(parameter, value);
				//is it a valid parameter
				if (validParameter->isValidParameter(parameter, commandName, value) != true) { return false; }
	
				
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
				if (parameter == "nexus" )		{ nexusfile = value; }
				if (parameter == "clustal" )	{ clustalfile = value; }
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
				if (parameter == "sorted" )			{ sorted = value;	}
				if (parameter == "trump" )          { trump = value; }
				if (parameter == "soft" )			{ soft = value; }
				if (parameter == "filter" )         { filter = value; }
				if (parameter == "scale" )			{ scale = value;	}
				if (parameter == "ends" )			{ ends = value; }
				if (parameter == "processors" )		{ processors = value;	}

			}
		}
		
		//make sure the user does not use both the line and label parameters
		if ((line != "") && (label != "")) { cout << "You may use either the line or label parameters, but not both." << endl; return false; }
		
		//check for valid files 
		if (commandName == "read.dist") { 
			validateReadFiles();
			validateReadDist();
		}else if (commandName == "read.otu") { 
			//you want to do shared commands
			if ((listfile != "") && (groupfile != ""))	{
				validateParseFiles(); //checks the listfile and groupfile parameters
			//you want to do single commands
			}else if ((listfile != "") || (rabundfile != "") || (sabundfile != "")){ 
				validateReadFiles();
				validateReadPhil();
			//you have not given a file
			}else if ((listfile == "") && (sharedfile == "") && (rabundfile == "") && (sabundfile == "")) {
				cout << "You must enter either a listfile, rabundfile, sabundfile or a sharedfile with the read.otu command. " << endl; return false; 
			//you want to do shared commands with a shared file
			}else if (sharedfile != "") {//you are reading a shared file
				validateReadFiles();
			}
		}else if (commandName == "read.tree") { 
			validateTreeFiles(); //checks the treefile and groupfile parameters
		}else if (commandName == "deconvolute") {
			if (fastafile == "") { cout << "You must enter a fastafile with the deconvolute() command." << endl; return false; }
			validateReadFiles();
		}
		
		//are you trying to cluster before you have read something	
		if (((commandName == "cluster") && (globaldata->gSparseMatrix == NULL)) ||
			((commandName == "cluster") && (globaldata->gListVector == NULL))) {
				cout << "Before you use the cluster command, you first need to read in a distance matrix." << endl; 
				errorFree = false;
		} 
		
		if ((commandName == "libshuff") && ((globaldata->gMatrix == NULL) || (globaldata->gGroupmap == NULL))) {
			 cout << "You must read in a matrix and groupfile using the read.dist command, before you use the libshuff command. " << endl; return false; 
		}
		
		if (commandName == "parsimony") {
			//are you trying to use parsimony without reading a tree or saying you want random distribution
			if (randomtree == "")  {
				if (globaldata->gTree.size() == 0) {
					cout << "You must read a treefile and a groupfile or set the randomtree parameter to the output filename you wish, before you may execute the parsimony command." << endl; return false;  }
			}
		}
		
		if ((commandName == "unifrac.weighted") || (commandName == "unifrac.unweighted") || (commandName == "concensus")) {
			if (globaldata->gTree.size() == 0) {//no trees were read
				cout << "You must execute the read.tree command, before you may execute the unifrac.weighted, unifrac.unweighted or concensus command." << endl; return false;  }
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
		
		if ((commandName == "collect.shared") || (commandName == "rarefaction.shared") || (commandName == "summary.shared") || (commandName == "tree.shared") || (commandName == "bootstrap.shared")){ 
			if (globaldata->getSharedFile() == "") {
				if (globaldata->getListFile() == "") { cout << "You must read a list and a group, or a shared before you can use the collect.shared, rarefaction.shared, summary.shared, tree.shared or bootstrap.shared commands." << endl; return false; }
				else if (globaldata->getGroupFile() == "") { cout << "You must read a list and a group, or a shared before you can use the collect.shared, rarefaction.shared, summary.shared, tree.shared or bootstrap.shared commands." << endl; return false; }
			}
		}
		
		if ((commandName == "heatmap") || (commandName == "venn")) { 
			if ((globaldata->getListFile() == "") && (globaldata->getSharedFile() == "")) {
				 cout << "You must read a list, or a list and a group, or a shared before you can use the heatmap or venn commands." << endl; return false; 
			}
		}
		
		if ((commandName == "filter.seqs") || (commandName == "dist.seqs")) { 
			if ((fastafile == "") && (nexusfile == "") && (clustalfile == "") && (phylipfile == "")) {
				 cout << "You must read either a fasta, nexus, clustal, or phylip file before you can use the filter.seqs command." << endl; return false; 
			}
			validateSeqsFiles();
		}
		
		if ((commandName == "bin.seqs")) { 
			if ((globaldata->getListFile() == "")) { cout << "You must read a list file before you can use the bin.seqs command." << endl; return false; }
			validateBinFiles();
		}
		
		if ((commandName == "get.oturep")) { 
			if ((globaldata->gSparseMatrix == NULL) || (globaldata->gListVector == NULL)) {
				cout << "Before you use the get.oturep command, you first need to read in a distance matrix." << endl; 
				errorFree = false;
			}
			if (listfile == "") { cout << "list is a required parameter for the get.oturep command." << endl; errorFree = false; }
			validateBinFiles();
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
		//are we reading a columnfile
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
//This function checks to make sure the user entered appropriate
// format parameters on a distfile read
void ErrorCheck::validateReadDist() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		if (groupfile != "") {
			ableToOpen = openInputFile(groupfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
		}
		
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
		}else if ((listfile == "") && (rabundfile == "") && (sabundfile == "") && (sharedfile == "")) {
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
//This function checks to make sure the user entered appropriate
// format parameters on a distfile read
void ErrorCheck::validateSeqsFiles() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		//checks to make sure only one file type is given
		if (phylipfile != "") { 
			if ((nexusfile != "") || (fastafile != "") || (clustalfile != "")) { 
				cout << "You may enter ONLY ONE of the following: phylip, fasta, nexus or clustal." << endl; errorFree = false; }
			else {
				ableToOpen = openInputFile(phylipfile, filehandle);
				filehandle.close();
				if (ableToOpen == 1) { //unable to open
					errorFree = false;
				}
			}
		}else if (nexusfile != "") { 
			if ((phylipfile != "") || (fastafile != "") || (clustalfile != "")) { 
				cout << "You may enter ONLY ONE of the following: phylip, fasta, nexus or clustal." << endl; errorFree = false; }
			else {
				ableToOpen = openInputFile(nexusfile, filehandle);
				filehandle.close();
				if (ableToOpen == 1) { //unable to open
					errorFree = false;
				}
			}
		}else if (fastafile != "") { 
			if ((phylipfile != "") || (nexusfile != "") || (clustalfile != "")) { 
				cout << "You may enter ONLY ONE of the following: phylip, fasta, nexus or clustal." << endl; errorFree = false; }
			else {
				ableToOpen = openInputFile(fastafile, filehandle);
				filehandle.close();
				if (ableToOpen == 1) { //unable to open
					errorFree = false;
				}
			}
		}else if (clustalfile != "") { 
			if ((phylipfile != "") || (nexusfile != "") || (fastafile != "")) { 
				cout << "You may enter ONLY ONE of the following: phylip, fasta, nexus or clustal." << endl; errorFree = false; }
			else {
				ableToOpen = openInputFile(clustalfile, filehandle);
				filehandle.close();
				if (ableToOpen == 1) { //unable to open
					errorFree = false;
				}
			}

		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateSeqsFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateSeqsFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*******************************************************/

/******************************************************/
//This function checks to make sure the user entered appropriate
// format parameters on a bin.seq command
void ErrorCheck::validateBinFiles() {
	try {
		ifstream filehandle;
		int ableToOpen;
		
		if (fastafile == "") {
				cout << "fasta is a required parameter for bin.seqs and get.oturep commands." << endl; errorFree = false; 
		}else if (fastafile != "") {
			//is it a valid filename'
			ableToOpen = openInputFile(fastafile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
		}else if (listfile != "") {
			//is it a valid filename'
			ableToOpen = openInputFile(listfile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
		}else if (globaldata->getNameFile() != "") {
			//is it a valid filename'
			ifstream filehandle;
			int ableToOpen = openInputFile(globaldata->getNameFile(), filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
		}else if (namefile != "") {
			//is it a valid filename'
			ifstream filehandle;
			int ableToOpen = openInputFile(namefile, filehandle);
			filehandle.close();
			//unable to open
			if (ableToOpen == 1) {  errorFree = false; }
		}


	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ErrorCheck class Function validateBinFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ErrorCheck class function validateBinFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
	fastafile       =   "";
	nexusfile       =   "";
	clustalfile     =   "";
	line			=	"";
	label			=	"";
	method			=   "furthest";
}
/*******************************************************/

/******************************************************/

