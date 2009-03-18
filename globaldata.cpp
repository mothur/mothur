#include "globaldata.hpp"
#include "sparsematrix.hpp"
#include "tree.h"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"

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
//This function parses through the option string of the command to remove its parameters
void GlobalData::parseGlobalData(string commandString, string optionText){
	try {
		commandName = commandString; //save command name to be used by other classes
		
		//set all non filename paramters to default
		reset();
		
		//clears out data from previous read
		if ((commandName == "read.dist") || (commandName == "read.otu") || (commandName == "read.tree")) { 
			clear();
			gGroupmap = NULL;
			gTree.clear();
			labels.clear(); lines.clear(); groups.clear();
			allLines = 1;
		}
		
		//saves help request
		if (commandName =="help") {
			helpRequest = optionText;
		}
		
		if (commandName == "libshuff") {
			iters = "10000";
			cutoff = "1.0";
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
				if (key == "shared" )	{ sharedfile = value; inputFileName = value; fileroot = value; format = "sharedfile";	}
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
				if (key == "abund" )        { abund = value;        }
				if (key == "random" )		{ randomtree = value;	}
				if (key == "calc")			{ calc = value;			}
				if (key == "step")			{ step = value;			}
				if (key == "form")			{ form = value;			}
				

				
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
			if (key == "shared" )	{ sharedfile = value; inputFileName = value; fileroot = value; format = "sharedfile";	} 
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
			if (key == "abund" )        { abund = value;        }
			if (key == "random" )		{ randomtree = value;	}
			if (key == "calc")			{ calc = value;		}
			if (key == "step")			{ step = value;			}
			if (key == "form")			{ form = value;			}

			if (key == "line") {//stores lines to be used in a vector
				lines.clear();
				line = value;
				label = "";
				if (line != "all") {  splitAtDash(value, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			if (key == "label") {//stores lines to be used in a vector
				labels.clear();
				label = value;
				line = "";
				if (label != "all") {  splitAtDash(value, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if (key == "groups") {//stores groups to be used in a vector
					Groups.clear();
					groups = value;
					splitAtDash(value, Groups);
			}
		}
		
		//set format for shared
		if ((listfile != "") && (groupfile != "")) { format = "shared"; }
		if ((phylipfile != "") && (groupfile != "")) { format = "matrix"; }
				
		//input defaults for calculators
		if (commandName == "collect.single") {
			if ((calc == "default") || (calc == "")) { calc = "sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "rarefaction.single") {
			if ((calc == "default") || (calc == "")) { calc = "sobs"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "collect.shared") {
			if ((calc == "default") || (calc == "")) { calc = "sharedsobs-sharedchao-sharedace-sharedjabund-sharedsorensonabund-sharedjclass-sharedsorclass-sharedjest-sharedsorest-sharedthetayc-sharedthetan"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "summary.single") {
			if ((calc == "default") || (calc == "")) { calc = "sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "summary.shared") {
			if ((calc == "default") || (calc == "")) { calc = "sharedsobs-sharedchao-sharedace-sharedjabund-sharedsorensonabund-sharedjclass-sharedsorclass-sharedjest-sharedsorest-sharedthetayc-sharedthetan"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
		}
		if (commandName == "rarefaction.shared") {
			if ((calc == "default") || (calc == "")) { calc = "sharedobserved"; }
			Estimators.clear();
			splitAtDash(calc, Estimators); 
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
string GlobalData::getSharedFile()		{	return sharedfile;	}
string GlobalData::getFastaFile()		{	return fastafile;	}
string GlobalData::getCutOff()			{	return cutoff;		}
string GlobalData::getFormat()			{	return format;		}
string GlobalData::getPrecision()		{	return precision;	}
string GlobalData::getMethod()			{	return method;		}
string GlobalData::getFileRoot()		{	return fileroot;	}
string GlobalData::getIters()			{	return iters;		}
string GlobalData::getJumble()			{	return jumble;		}
string GlobalData::getFreq()			{	return freq;		}
string GlobalData::getAbund()           {   return abund;       }
string GlobalData::getRandomTree()		{	return randomtree;	}
string GlobalData::getGroups()			{	return groups;		}
string GlobalData::getStep()			{	return step;		}
string GlobalData::getForm()			{	return form;		}
void GlobalData::setListFile(string file)	{	listfile = file;	inputFileName = file;}
void GlobalData::setRabundFile(string file)	{	rabundfile = file;	inputFileName = file;}
void GlobalData::setSabundFile(string file)	{	sabundfile = file;	inputFileName = file;}
void GlobalData::setPhylipFile(string file)	{	phylipfile = file;    inputFileName = file;}
void GlobalData::setColumnFile(string file)	{	columnfile = file;    inputFileName = file;}
void GlobalData::setNameFile(string file)		{	namefile = file;		}
void GlobalData::setFormat(string Format)		{	format = Format;		}
void GlobalData::setRandomTree(string Random)	{	randomtree = Random;	}
void GlobalData::setGroups(string g)			{	groups = g;				}
void GlobalData::setCalc(string Calc)			{	calc = Calc;			}

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
	sharedfile		=	"";
	cutoff			=	"10.00";
	format			=	"";
	precision		=	"100";
	iters			=	"1000"; 
	line			=   "";
	label			=	"";
	groups			=	"";
	jumble			=	"1";	//0 means don't jumble, 1 means jumble.
	randomtree		=	"";  //"" means user will enter some user trees, "outputfile" means they just want the random tree distribution to be outputted to outputfile.
	freq			=	"100";
	method			=	"furthest";
	fileroot		=	"";
	abund           =   "10";
	step			=	"0.01";
	form			=	"integral";
}

//*******************************************************/

/******************************************************/
void GlobalData::reset() {
	cutoff			=	"10.00";
	precision		=	"100";
	iters			=	"1000"; 
	groups			=	"";
	jumble			=	"1";	//0 means don't jumble, 1 means jumble.
	randomtree		=	"";  //"" means user will enter some user trees, "outputfile" means they just want the random tree distribution to be outputted to outputfile.
	freq			=	"100";
	method			=	"furthest";
	calc			=	"";
	abund			=   "10";
	step			=	"0.01";
	form			=	"integral";
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
