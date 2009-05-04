#ifndef GLOBALDATA_HPP
#define GLOBALDATA_HPP

#include "mothur.h"
#include "groupmap.h"
#include "treemap.h"

#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"


using namespace std;

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
	string inputFileName, helpRequest, commandName, vertical;
	bool allLines;
	vector<string>  Estimators, Groups; //holds estimators to be used
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	vector<string> Treenames;
	
	string getPhylipFile();
	string getColumnFile();
	string getListFile();
	string getRabundFile();
	string getSabundFile();
	string getNameFile();
	string getGroupFile();
	string getOrderFile();
	string getFastaFile();
	string getNexusFile();
	string getClustalFile();
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
	string getGroups();
	string getStep();
	string getForm();
	string getSorted();

	string getTrump();
	string getSoft();
	string getFilter();
	

	string getScale();


	void setListFile(string);
	void setPhylipFile(string);
	void setColumnFile(string);
	void setNameFile(string);
	void setRabundFile(string);
	void setSabundFile(string);
	void setFormat(string);
	void setRandomTree(string);
	void setGroups(string);
	void setCalc(string);

	void clear(); 
	void clearLabels();
	void clearAbund();
	
	void parseGlobalData(string, string);
	
	void parseTreeFile();		//parses through tree file to find names of nodes and number of them
							//this is required in case user has sequences in the names file that are
							//not included in the tree. 
							//only takes names from the first tree in the tree file and assumes that all trees use the same names.

		
private:

	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, orderfile, fastafile, nexusfile, clustalfile, treefile, sharedfile, line, label, randomtree, groups;
	string cutoff, format, precision, method, fileroot, iters, jumble, freq, calc, abund, step, form, sorted, trump, soft, filter, scale;


	static GlobalData* _uniqueInstance;
	GlobalData( const GlobalData& ); // Disable copy constructor
	void operator=( const GlobalData& ); // Disable assignment operator
	GlobalData();
	~GlobalData();
	void reset();	//clears all non filename parameters
	void readTreeString(ifstream&);
	
	
	
};

#endif
