#ifndef GLOBALDATA_HPP
#define GLOBALDATA_HPP

#include <string>
#include <vector>
#include <set>
#include "groupmap.h"

using namespace std;

class ListVector;
class SparseMatrix;
class OrderVector;
class InputData;
class GroupMap;
class SAbundVector;

class GlobalData {
public:
	static GlobalData* getInstance();
	ListVector* getListVector();
	SparseMatrix* getSparseMatrix();
	InputData* ginput;
	OrderVector* gorder;
	ListVector* glist;
	SAbundVector* sabund;
	GroupMap* gGroupmap;
	string inputFileName, helpRequest, commandName;
	bool allLines;
	vector<string> singleEstimators, summaryEstimators, sharedEstimators, rareEstimators, sharedRareEstimators, sharedSummaryEstimators; //holds estimators to be used
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
	string getCutOff();
	string getFormat();
	string getPrecision();
	string getMethod();
	string getFileRoot();
	string getIters();
	string getJumble();
	string getFreq();
	void setListFile(string);
	void setPhylipFile(string);
	void setColumnFile(string);
	void setNameFile(string);
	void setRabundFile(string);
	void setSabundFile(string);
	void setFormat(string);

	
	void setListVector(ListVector*);
	void setSparseMatrix(SparseMatrix*);
	void clear(); 
	
	void parseGlobalData(string, string);
	void splitAtEquals(string&, string&);
	void splitAtComma(string&, string&);
	void splitAtDash(string&, vector<string>&);
	void splitAtDash(string&, set<int>&);
	void splitAtDash(string&, set<string>&);
	
private:
	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, orderfile, fastafile,line, label;
	string cutoff, format, precision, method, fileroot, iters, jumble, freq, single, rarefaction, shared, summary, sharedsummary, sharedrarefaction;
	static GlobalData* _uniqueInstance;
	GlobalData( const GlobalData& ); // Disable copy constructor
	void operator=( const GlobalData& ); // Disable assignment operator
	GlobalData();
	~GlobalData();
	ListVector* gListVector;
	SparseMatrix* gSparseMatrix;
	
};

//**********************************************************************************************************************

#endif
