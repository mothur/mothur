

#include "globaldata.hpp"
#include "sharedlistvector.h"
#include "inputdata.h"
#include "fullmatrix.h"

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
string GlobalData::getFormat()			{	return format;			}

void GlobalData::setListFile(string file)		{	listfile = file;	inputFileName = file;					}
void GlobalData::setTreeFile(string file)		{	treefile = file;	inputFileName = file;					}
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
	m = MothurOut::getInstance();
	//option definitions should go here...
	clear();
	gListVector = NULL;		
	gSparseMatrix = NULL;	
	ginput = NULL;
	gorder = NULL;
	glist = NULL;
	gSharedList = NULL;
	sabund = NULL;
	rabund = NULL;
	gGroupmap = NULL;
	gMatrix = NULL;
	gTreemap = NULL;
	gSequenceDB = NULL;
	nameMap = NULL;
}
/*******************************************************/

/******************************************************/
void GlobalData::clear() {
	//option definitions should go here...
	phylipfile		=	""; //do we need this?
	columnfile		=	""; //do we need this?
	listfile		=	"";
	rabundfile		=	"";
	sabundfile		=	"";
	namefile		=	""; //do we need this?
	groupfile		=	""; //do we need this?
	orderfile		=	"";
//	fastafile		=   ""; //do we need this?
	treefile		=	"";
	sharedfile		=	"";
	format = "";
}


/*******************************************************/

/******************************************************/
void GlobalData::newRead() {
	try{	
			//remove old file names
			clear();
			
			//free memory
			if (gGroupmap != NULL) { delete gGroupmap; gGroupmap = NULL; }

			if (gListVector != NULL) { delete gListVector; gListVector = NULL;}

			if (gSparseMatrix != NULL) { delete gSparseMatrix; gSparseMatrix = NULL; }

			if (ginput != NULL) { delete ginput; ginput = NULL;}

			if (gorder != NULL) { delete gorder; gorder = NULL; }

			if (glist != NULL) { delete glist; glist = NULL;}

			if (gSharedList != NULL) { delete gSharedList; gSharedList = NULL; }

			if (sabund != NULL) { delete sabund; sabund = NULL;}

			if (rabund != NULL) { delete rabund; rabund = NULL; }

			if (gMatrix != NULL) { delete gMatrix; gMatrix = NULL;}

			if (gTreemap != NULL) { delete gTreemap; gTreemap = NULL; }

			if (gSequenceDB != NULL) { delete gSequenceDB; gSequenceDB = NULL;}
			
			if (nameMap != NULL) { delete nameMap; nameMap = NULL; }


			gTree.clear();
			Treenames.clear();
			labels.clear(); Groups.clear();
			allLines = 1;
			runParse = true;
			names.clear();
	}
	catch(exception& e) {
		m->errorOut(e, "GlobalData", "newRead");
		exit(1);
	}
}

//******************************************************/

/******************************************************/
GlobalData::~GlobalData() {
	_uniqueInstance = 0;
	try {
		if (gGroupmap != NULL) { delete gGroupmap; gGroupmap = NULL; }
		if (gListVector != NULL) { delete gListVector; gListVector = NULL;}
		if (gSparseMatrix != NULL) { delete gSparseMatrix; gSparseMatrix = NULL; }
		if (ginput != NULL) { delete ginput; ginput = NULL;}
		if (gorder != NULL) { delete gorder; gorder = NULL; }
		if (glist != NULL) { delete glist; glist = NULL;}
		if (gSharedList != NULL) { delete gSharedList; gSharedList = NULL; }
		if (sabund != NULL) { delete sabund; sabund = NULL;}
		if (rabund != NULL) { delete rabund; rabund = NULL; }
		if (gMatrix != NULL) { delete gMatrix; gMatrix = NULL;}
		if (gTreemap != NULL) { delete gTreemap; gTreemap = NULL; }
		if (gSequenceDB != NULL) { delete gSequenceDB; gSequenceDB = NULL;}
		if (nameMap != NULL) { delete nameMap; nameMap = NULL; }
	}
	catch(exception& e) {
		m->errorOut(e, "GlobalData", "~GlobalData");
		exit(1);
	}
}
/*******************************************************/

/*******************************************************/


