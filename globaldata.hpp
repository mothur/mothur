#ifndef GLOBALDATA_HPP
#define GLOBALDATA_HPP

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <exception>

#include "groupmap.h"
#include "treemap.h"

using namespace std;

class ListVector;
class SharedListVector;
class SparseMatrix;
class Tree;
class OrderVector;
class InputData;
class GroupMap;
class TreeMap;
class SAbundVector;

class GlobalData {
public:
	static GlobalData* getInstance();
	ListVector* getListVector();
	SparseMatrix* getSparseMatrix();
	InputData* ginput;
	OrderVector* gorder;
	ListVector* glist;
	vector<Tree*> gTree;
	SharedListVector* gSharedList;
	SAbundVector* sabund;
	GroupMap* gGroupmap;
	TreeMap* gTreemap;
	string inputFileName, helpRequest, commandName;
	bool allLines;
	vector<string>  Estimators, Groups; //holds estimators to be used
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	
	string getPhylipFile();
	string getColumnFile();
	string getListFile();
	string getRabundFile();
	string getSabundFile();
	string getNameFile();
	string getGroupFile();
	string getOrderFile();
	string getFastaFile();
	string getTreeFile();
	string getSharedFile();
	string getCutOff();
	string getFormat();
	string getPrecision();
	string getMethod();
	string getFileRoot();
	string getIters();
	string getJumble();
	string getFreq();
	string getAbund();
	string getRandomTree();

	void setListFile(string);
	void setPhylipFile(string);
	void setColumnFile(string);
	void setNameFile(string);
	void setRabundFile(string);
	void setSabundFile(string);
	void setFormat(string);
	void setRandomTree(string);
	void setCalc(string);

	
	void setListVector(ListVector*);
	void setSparseMatrix(SparseMatrix*);
	void clear(); 
	void clearLabels();
	void clearAbund();
	
	void parseGlobalData(string, string);
		
private:
	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, orderfile, fastafile, treefile, sharedfile, line, label, randomtree, groups;
	string cutoff, format, precision, method, fileroot, iters, jumble, freq, calc, abund;

	static GlobalData* _uniqueInstance;
	GlobalData( const GlobalData& ); // Disable copy constructor
	void operator=( const GlobalData& ); // Disable assignment operator
	GlobalData();
	~GlobalData();
	ListVector* gListVector;
	SparseMatrix* gSparseMatrix;
	void clear();  //clears all parameters
	void reset();	//clears all non filename parameters
	
	
	
};

#endif
