#ifndef GLOBALDATA_HPP
#define GLOBALDATA_HPP

#include "mothur.h"
#include "groupmap.h"
#include "treemap.h"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"

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
	string inputFileName, helpRequest, commandName, vertical, argv;
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
	string getCountEnds();
	string getProcessors();
	string getSize();
	string getCandidateFile();
	string getSearch();
	string getKSize();
	string getAlign();
	string getMatch();
	string getMismatch();
	string getGapopen();
	string getGapextend();
	string getVertical();
	string getTrump();
	string getSoft();
	string getHard();
	string getScale();
	string getStartPos();
	string getEndPos();
	string getMaxAmbig();
	string getMaxHomoPolymer();
	string getMinLength();
	string getMaxLength();

	void setListFile(string);
	void setGroupFile(string file);	
	void setPhylipFile(string);
	void setColumnFile(string);
	void setNameFile(string);
	void setRabundFile(string);
	void setSabundFile(string);
	void setSharedFile(string);
	void setFormat(string);
	void setRandomTree(string);
	void setGroups(string);
	void setCalc(string);
	void setCountEnds(string);
	void setProcessors(string);

	void clear(); 
	void clearLabels();
	void clearAbund();
	
	void parseGlobalData(string, string);
	
	void parseTreeFile();	//parses through tree file to find names of nodes and number of them
							//this is required in case user has sequences in the names file that are
							//not included in the tree. 
							//only takes names from the first tree in the tree file and assumes that all trees use the same names.

		
private:

	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, orderfile, fastafile, treefile, sharedfile, line, label, randomtree, groups, cutoff, format, precision, method, fileroot, iters, jumble, freq, calc, abund, step, form, sorted, trump, soft, hard, scale, countends, processors, candidatefile, search, ksize, align, match, size, mismatch, gapopen, gapextend, minLength, maxLength, startPos, endPos, maxAmbig, maxHomoPolymer;


	static GlobalData* _uniqueInstance;
	GlobalData( const GlobalData& ); // Disable copy constructor
	void operator=( const GlobalData& ); // Disable assignment operator
	GlobalData();
	~GlobalData();
	void reset();	//clears all non filename parameters
	void readTreeString(ifstream&);
	
	
	
};

#endif
