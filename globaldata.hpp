#ifndef GLOBALDATA_HPP
#define GLOBALDATA_HPP

#include "mothur.h"
#include "groupmap.h"
#include "treemap.h"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "tree.h"
#include "sparsematrix.hpp"
#include "sequencedb.h"


class ListVector;
class SharedListVector;
class SparseMatrix;
class FullMatrix;
class Tree;
class OrderVector;
class InputData;
class GroupMap;
class TreeMap;
class SAbundVector;
class RAbundVector;
class SequenceDB;

class GlobalData {
public:
	static GlobalData* getInstance();
	ListVector* gListVector;
	SparseMatrix* gSparseMatrix;
	InputData* ginput;
	OrderVector* gorder;
	ListVector* glist;
	vector<Tree*> gTree;
	SharedListVector* gSharedList;
	SAbundVector* sabund;
	RAbundVector* rabund;
	GroupMap* gGroupmap;
	FullMatrix* gMatrix;
	TreeMap* gTreemap;
	SequenceDB* gSequenceDB;
	string inputFileName, argv;
	bool allLines, runParse;
	vector<string>  Estimators, Groups; //holds estimators to be used
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	vector<string> Treenames;
	
	
	string getPhylipFile();
	string getColumnFile();
	string getListFile();
	string getRabundFile();
	string getSabundFile();
	string getNameFile();	//do we need this?
	string getGroupFile();	//do we need this?
	string getOrderFile();
//	string getFastaFile();
	string getTreeFile();
	string getSharedFile();
	string getFormat();	//do we need this?
//	string getCandidateFile();
//	string getTemplateFile();

	void setListFile(string);
//	void setFastaFile(string);
	void setTreeFile(string);
//	void setCandidateFile(string);
//	void setTemplateFile(string);
	void setGroupFile(string);		//do we need this?
	void setPhylipFile(string);
	void setColumnFile(string);
	void setNameFile(string);	//do we need this?
	void setRabundFile(string);
	void setSabundFile(string);
	void setSharedFile(string);
	void setOrderFile(string file);
	void setFormat(string);	//do we need this?
	
	void clear(); 
	void clearLabels();
	void clearAbund();
	
	void newRead();
	
private:

	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, orderfile, treefile, sharedfile, format;

	static GlobalData* _uniqueInstance;
	GlobalData( const GlobalData& ); // Disable copy constructor
	void operator=( const GlobalData& ); // Disable assignment operator
	GlobalData();
	~GlobalData();
	
	
};

#endif
