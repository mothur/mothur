#include "globaldata.hpp"
#include "tree.h"
#include "sparsematrix.hpp"

/*******************************************************/

/******************************************************/
GlobalData* GlobalData::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new GlobalData();
	}
	return _uniqueInstance;
}
/*******************************************************/

/******************************************************/
// These functions give you the option parameters of the commands
string GlobalData::getPhylipFile()		{	return phylipfile;		}
string GlobalData::getColumnFile()		{	return columnfile;		}
string GlobalData::getListFile()		{	return listfile;		}
string GlobalData::getRabundFile()		{	return rabundfile;		}
string GlobalData::getSabundFile()		{	return sabundfile;		}
string GlobalData::getNameFile()		{	return namefile;		}
string GlobalData::getGroupFile()		{	return groupfile;		}
string GlobalData::getOrderFile()		{	return orderfile;		}
string GlobalData::getTreeFile()		{	return treefile;		}
string GlobalData::getSharedFile()		{	return sharedfile;		}
string GlobalData::getFastaFile()		{	return fastafile;		}	
string GlobalData::getFormat()			{	return format;			}
string GlobalData::getCandidateFile()	{	return candidatefile;	}


void GlobalData::setListFile(string file)		{	listfile = file;	inputFileName = file;					}
void GlobalData::setFastaFile(string file)		{	fastafile = file;	inputFileName = file;					}
void GlobalData::setTreeFile(string file)		{	treefile = file;	inputFileName = file;					}
void GlobalData::setCandidateFile(string file)	{	candidatefile = file;										}
void GlobalData::setRabundFile(string file)		{	rabundfile = file;	inputFileName = file;					}
void GlobalData::setSabundFile(string file)		{	sabundfile = file;	inputFileName = file;					}
void GlobalData::setPhylipFile(string file)		{	phylipfile = file;    inputFileName = file;					}
void GlobalData::setColumnFile(string file)		{	columnfile = file;    inputFileName = file;					}
void GlobalData::setGroupFile(string file)		{	groupfile = file;											}
void GlobalData::setSharedFile(string file)		{	sharedfile = file;	inputFileName = file;					}
void GlobalData::setNameFile(string file)		{	namefile = file;		}
void GlobalData::setOrderFile(string file)		{	orderfile = file;		}
void GlobalData::setFormat(string Format)		{	format = Format;		}


/*******************************************************/

/******************************************************/
GlobalData::GlobalData() {
	//option definitions should go here...
	clear();
	gListVector = NULL;		
	gSparseMatrix = NULL;	
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
	candidatefile	=	"";
}
/*******************************************************/

/******************************************************/
void GlobalData::newRead() {
	try{
			clear();
			gGroupmap = NULL;
			gListVector = NULL;
			gSparseMatrix = NULL;
			gTree.clear();
			Treenames.clear();
			labels.clear(); lines.clear(); Groups.clear();
			allLines = 1;
			runParse = true;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GlobalData class Function newRead. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GlobalData class function newRead. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//******************************************************/

/******************************************************/
GlobalData::~GlobalData() {
	_uniqueInstance = 0;
	if(gListVector != NULL)		{	delete gListVector;		}
	if(gSparseMatrix != NULL)	{	delete gSparseMatrix;	}
	if(gorder != NULL)			{	delete gorder;		}
}
/*******************************************************/

/*******************************************************/


