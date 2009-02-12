#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace std;

#include "globaldata.hpp"
#include "sparsematrix.hpp"
#include "tree.h"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include <exception>
#include <iostream>

/*******************************************************/

/******************************************************/
GlobalData* GlobalData::getInstance() {
	if( _uniqueInstance == 0 ) {
		_uniqueInstance = new GlobalData();
	}
	return _uniqueInstance;
}
/*******************************************************/

/******************************************************/

ListVector* GlobalData::getListVector()		{	return gListVector;		}
/*******************************************************/

/******************************************************/
void GlobalData::setListVector(ListVector* lv){
	try {
		if(gListVector != NULL){	delete gListVector;	}
		gListVector = new ListVector(*lv);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function setListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function setListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*******************************************************/

/******************************************************/

SparseMatrix* GlobalData::getSparseMatrix()	{	return gSparseMatrix;	}
/*******************************************************/

/******************************************************/
void GlobalData::setSparseMatrix(SparseMatrix* sm){
	try{
		if(gSparseMatrix != NULL){	delete gSparseMatrix;	}
		gSparseMatrix = new SparseMatrix(*sm);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function setSparseMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function setSparseMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/
//This function parses through the option string of the command to remove its parameters
void GlobalData::parseGlobalData(string commandString, string optionText){
	try {
		allLines = 1;
		commandName = commandString; //save command name to be used by other classes
		
		//clears out data from previous read
		if ((commandName == "read.dist") || (commandName == "read.otu") || (commandName == "read.tree")) { 
			clear();
		}
		
		//saves help request
		if (commandName =="help") {
			helpRequest = optionText;
		}
		
		string key, value;		
		//reads in parameters and values
		if((optionText != "") && (commandName != "help")){
			while((optionText.find_first_of(',') != -1)) {  //while there are parameters
				splitAtComma(value, optionText);
				splitAtEquals(key, value);
				
				if (key == "phylip" )	{ phylipfile = value; inputFileName = value; fileroot = value; format = "phylip";	}
				if (key == "column" )	{ columnfile = value; inputFileName = value; fileroot = value; format = "column";	}
				if (key == "list" )		{ listfile = value; inputFileName = value; fileroot = value; format = "list";		}
				if (key == "rabund" )	{ rabundfile = value; inputFileName = value; fileroot = value; format = "rabund";	}
				if (key == "sabund" )	{ sabundfile = value; inputFileName = value; fileroot = value; format = "sabund";	} 
				if (key == "fasta" )	{ fastafile = value; inputFileName = value; fileroot = value; format = "fasta";		} 
				if (key == "tree" )		{ treefile = value; inputFileName = value; fileroot = value; format = "tree";		}
				if (key == "name" )		{ namefile = value;		}
				if (key == "order" )	{ orderfile = value;	}
				if (key == "group" )	{ groupfile = value;	}
				if (key == "cutoff" )		{ cutoff = value;		}
				if (key == "precision" )	{ precision = value;	}
				if (key == "iters" )		{ iters = value;		}
				if (key == "jumble" )		{ jumble = value;		}
				if (key == "freq" )			{ freq = value;			}
				if (key == "method" )		{ method = value;		}
				if (key == "fileroot" )		{ fileroot = value;		}
				if (key == "randomtree" )	{ randomtree = value;	}
				if (key == "groups" )		{ groups = value;	}
				
				if (key == "single") {//stores estimators in a vector
					singleEstimators.clear(); //clears out old values
					if (value == "default") { value = "sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson-rarefraction"; }
					splitAtDash(value, singleEstimators);
				}
				if (key == "rarefaction") {//stores estimators in a vector
					rareEstimators.clear(); //clears out old values
					if (value == "default") { value = "rarefraction"; }
					splitAtDash(value, rareEstimators);
				}
				if (key == "shared") {//stores estimators in a vector
					sharedEstimators.clear(); //clears out old values
					if (value == "default") { value = "sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN"; }
					splitAtDash(value, sharedEstimators);
				}
				if (key == "summary") { //stores summaries to be used in a vector
					summaryEstimators.clear();
					if (value == "default") { value = "summary-chao-ace-jack-bootstrap-shannon-npshannon-simpson"; }
					splitAtDash(value, summaryEstimators);
				}
				if (key == "sharedsummary") { //stores sharedSummaries to be used in a vector
					sharedSummaryEstimators.clear();
					if (value == "default") { value = "sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN"; }
					splitAtDash(value, sharedSummaryEstimators);
				}
				if (key == "sharedrarefaction") { //stores sharedrarefaction to be used in a vector
					sharedRareEstimators.clear();
					if (value == "default") { value = "sharedobserved"; }
					splitAtDash(value, sharedRareEstimators);
				}
				if (key == "line") {//stores lines to be used in a set
					lines.clear();
					line = value;
					label = "";
					splitAtDash(value, lines);
					allLines = 0;
				}
				if (key == "label") {//stores labels to be used in a set
					labels.clear();
					label = value;
					line = "";
					splitAtDash(value, labels);
					allLines = 0;
				}
				if (key == "groups") {//stores groups to be used in a vector
					Groups.clear();
					groups = value;
					splitAtDash(value, Groups);
				}

			}
			
			//saves the last parameter
			value = optionText;
			splitAtEquals(key, value);
			if (key == "phylip" )	{ phylipfile = value; inputFileName = value; fileroot = value; format = "phylip";	}
			if (key == "column" )	{ columnfile = value; inputFileName = value; fileroot = value; format = "column";	}
			if (key == "list" )		{ listfile = value; inputFileName = value; fileroot = value; format = "list";		}
			if (key == "rabund" )	{ rabundfile = value; inputFileName = value; fileroot = value; format = "rabund";	}
			if (key == "sabund" )	{ sabundfile = value; inputFileName = value; fileroot = value; format = "sabund";	}
			if (key == "fasta" )	{ fastafile = value; inputFileName = value; fileroot = value; format = "fasta";		}
			if (key == "tree" )		{ treefile = value; inputFileName = value; fileroot = value; format = "tree";		}  
			if (key == "name" )		{ namefile = value;		}
			if (key == "order" )	{ orderfile = value;	}
			if (key == "group" )	{ groupfile = value;	}
			if (key == "cutoff" )		{ cutoff = value;		}
			if (key == "precision" )	{ precision = value;	}
			if (key == "iters" )		{ iters = value;		}
			if (key == "jumble" )		{ jumble = value;		}
			if (key == "freq" )			{ freq = value;			}
			if (key == "method" )		{ method = value;		}
			if (key == "fileroot" )		{ fileroot = value;		}
			if (key == "randomtree" )	{ randomtree = value;	}
			if (key == "groups" )		{ groups = value;	}

			
			if (key == "single") {//stores estimators in a vector
				singleEstimators.clear(); //clears out old values
				if (value == "default") { value = "sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson-rarefraction"; }
				splitAtDash(value, singleEstimators);
			}
			if (key == "rarefaction") {//stores estimators in a vector
				rareEstimators.clear(); //clears out old values
				if (value == "default") { value = "rarefraction"; }
				splitAtDash(value, rareEstimators);
			}
			if (key == "shared") {//stores estimators in a vector
				sharedEstimators.clear(); //clears out old values
				if (value == "default") { value = "sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN"; }
				splitAtDash(value, sharedEstimators);
			}
			if (key == "summary") { //stores summaries to be used in a vector
				summaryEstimators.clear();
				if (value == "default") { value = "summary-chao-ace-jack-bootstrap-shannon-npshannon-simpson"; }
				splitAtDash(value, summaryEstimators);
			}
			if (key == "sharedsummary") { //stores sharedSummaries to be used in a vector
				sharedSummaryEstimators.clear();
				if (value == "default") { value = "sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN"; }
				splitAtDash(value, sharedSummaryEstimators);
			}
			if (key == "sharedrarefaction") { //stores sharedrarefaction to be used in a vector
				sharedRareEstimators.clear();
				if (value == "default") { value = "sharedobserved"; }
				splitAtDash(value, sharedRareEstimators);
			}
			
			if (key == "line") {//stores lines to be used in a vector
				lines.clear();
				line = value;
				label = "";
				splitAtDash(value, lines);
				allLines = 0;
			}
			if (key == "label") {//stores lines to be used in a vector
				labels.clear();
				label = value;
				line = "";
				splitAtDash(value, labels);
				allLines = 0;
			}
			if (key == "groups") {//stores groups to be used in a vector
					Groups.clear();
					groups = value;
					splitAtDash(value, Groups);
			}

		}
		
		//set format for shared
		if ((listfile != "") && (groupfile != "")) { format = "shared"; }
				
		//input defaults
		if (commandName == "collect.single") {
			if (singleEstimators.size() == 0) { splitAtDash(single, singleEstimators); }
		}
		if (commandName == "rarefaction.single") {
			if (rareEstimators.size() == 0) { splitAtDash(rarefaction, rareEstimators);  }	
		}
		if (commandName == "collect.shared") {
			if (sharedEstimators.size() == 0) { splitAtDash(shared, sharedEstimators); }	
		}
		if (commandName == "summary.single") {
			if (summaryEstimators.size() == 0) { splitAtDash(summary, summaryEstimators); }
		}
		if (commandName == "summary.shared") {
			if (sharedSummaryEstimators.size() == 0) { splitAtDash(sharedsummary, sharedSummaryEstimators); }
		}
		if (commandName == "rarefaction.shared") {
			if (sharedRareEstimators.size() == 0) { splitAtDash(sharedrarefaction, sharedRareEstimators); }
		}


		//if you have done a read.otu with a groupfile but don't want to use it anymore because you want to do single commands
		if ((commandName == "collect.single") || (commandName == "rarefaction.single") || (commandName == "summary.single")) {
			if (listfile != "") { format = "list"; }
			else if (sabundfile != "") { format = "sabund"; }
			else if (rabundfile != "") { format = "rabund"; }
		}
				
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function parseGlobalData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function parseGlobalData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/
// These functions give you the option parameters of the commands
string GlobalData::getPhylipFile()		{	return phylipfile;	}
string GlobalData::getColumnFile()		{	return columnfile;	}
string GlobalData::getListFile()		{	return listfile;	}
string GlobalData::getRabundFile()		{	return rabundfile;	}
string GlobalData::getSabundFile()		{	return sabundfile;	}
string GlobalData::getNameFile()		{	return namefile;	}
string GlobalData::getGroupFile()		{	return groupfile;	}
string GlobalData::getOrderFile()		{	return orderfile;	}
string GlobalData::getTreeFile()		{	return treefile;	}
string GlobalData::getFastaFile()		{	return fastafile;	}
string GlobalData::getCutOff()			{	return cutoff;		}
string GlobalData::getFormat()			{	return format;		}
string GlobalData::getPrecision()		{	return precision;	}
string GlobalData::getMethod()			{	return method;		}
string GlobalData::getFileRoot()		{	return fileroot;	}
string GlobalData::getIters()			{	return iters;		}
string GlobalData::getJumble()			{	return jumble;		}
string GlobalData::getFreq()			{	return freq;		}
string GlobalData::getRandomTree()		{	return randomtree;	}
void GlobalData::setListFile(string file)	{	listfile = file;	inputFileName = file;}
void GlobalData::setRabundFile(string file)	{	rabundfile = file;	inputFileName = file;}
void GlobalData::setSabundFile(string file)	{	sabundfile = file;	inputFileName = file;}
void GlobalData::setPhylipFile(string file)	{	phylipfile = file;    inputFileName = file;}
void GlobalData::setColumnFile(string file)	{	columnfile = file;    inputFileName = file;}
//void GlobalData::setGroupFile(string file)	{	groupfile = file;		}
void GlobalData::setNameFile(string file)		{	namefile = file;		}
void GlobalData::setFormat(string Format)		{	format = Format;		}
void GlobalData::setRandomTree(string Random)	{	randomtree = Random;	}


/*******************************************************/

/******************************************************/

GlobalData::GlobalData() {
	//option definitions should go here...
	helpRequest = "";
	clear();
}
/*******************************************************/

/******************************************************/

void GlobalData::clear() {
	//option definitions should go here...
	phylipfile		=	"";
	columnfile		=	"";
	listfile		=	"";
	rabundfile		=	"";
	sabundfile		=	"";
	namefile		=	"";
	groupfile		=	""; 
	orderfile		=	"";
	fastafile		=   "";
	treefile		=	"";
	cutoff			=	"10.00";
	format			=	"";
	precision		=	"100";
	iters			=	"1000"; 
	line			=   "";
	label			=	"";
	groups			=	"";
	jumble			=	"1";	//0 means don't jumble, 1 means jumble.
	randomtree		=	"0";  //0 means user will enter some user trees, 1 means they just want the random tree distribution.
	freq			=	"100";
	method			=	"furthest";
	fileroot		=	"";
	single			=	"sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson-rarefraction";
	rarefaction		=	"sobs";
	shared			=	"sharedSobs-sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN";
	sharedsummary	=   "sharedSobs-sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN";
	summary			=	"summary-chao-ace-jack-bootstrap-shannon-npshannon-simpson";
	sharedrarefaction = "sharedobserved";
}
/*******************************************************/

/******************************************************/

GlobalData::~GlobalData() {
	_uniqueInstance = 0;
	if(gListVector != NULL)		{	delete gListVector;		}
	if(gSparseMatrix != NULL)	{	delete gSparseMatrix;	}
	if(gorder != NULL)			{	delete gorder;		}
}
/*******************************************************/

/******************************************************/
//This function parses the estimator options and puts them in a vector
void GlobalData::splitAtDash(string& estim, vector<string>& container) {
	try {
		string individual;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				container.push_back(individual);
			}
		}
		//get last one
		container.push_back(estim);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/
//This function parses the label options and puts them in a set
void GlobalData::splitAtDash(string& estim, set<string>& container) {
	try {
		string individual;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				container.insert(individual);
			}
		}
		//get last one
		container.insert(estim);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/
//This function parses the line options and puts them in a set
void GlobalData::splitAtDash(string& estim, set<int>& container) {
	try {
		string individual;
		int lineNum;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				convert(individual, lineNum); //convert the string to int
				container.insert(lineNum);
			}
		}
		//get last one
		convert(estim, lineNum); //convert the string to int
		container.insert(lineNum);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/

//This function splits up the various option parameters
void GlobalData::splitAtComma(string& prefix, string& suffix){
	try {
		prefix = suffix.substr(0,suffix.find_first_of(','));
		if ((suffix.find_first_of(',')+2) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
			suffix = suffix.substr(suffix.find_first_of(',')+2, suffix.length());
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function splitAtComma. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function splitAtComma. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/
//This function separates the key value from the option value i.e. dist=96_...
void GlobalData::splitAtEquals(string& key, string& value){		
	try {
		if(value.find_first_of('=') != -1){
			key = value.substr(0,value.find_first_of('='));
			if ((value.find_first_of('=')+1) <= value.length()) {
				value = value.substr(value.find_first_of('=')+1, value.length());
			}
		}else{
			key = value;
			value = 1;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function splitAtEquals. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function splitAtEquals. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/*******************************************************/

/******************************************************/
