#ifndef FILEOUTPUT_H
#define FILEOUTPUT_H

#include "mothur.h"
#include "mothurout.h"

/***********************************************************************/

class FileOutput {
	
public:
	FileOutput(){ m = MothurOut::getInstance(); }
	virtual ~FileOutput(){};
	
	virtual void initFile(string) = 0;
	virtual void initFile(string, vector<string>) = 0;
	virtual void output(int, vector<double>) = 0;
	virtual void output(vector<double>) = 0;
	virtual void resetFile() = 0;
	virtual string getFileName() = 0;

protected:
	int renameOk;
	MothurOut* m;

};	
	
/***********************************************************************/

class ThreeColumnFile : public FileOutput {
	
public:
	ThreeColumnFile(string n) : FileOutput(), inName(n), counter(0), outName(n + ".temp") { };
	~ThreeColumnFile();
	void initFile(string);
	void output(int, vector<double>);
	void resetFile();
	string getFileName()	{ return inName;	};
	
	void initFile(string, vector<string>){};
	void output(vector<double>) {};

private:
	string inName;
	string outName;
	ifstream inFile;
	ofstream outFile;
	int counter;
};


/***********************************************************************/
class OneColumnFile : public FileOutput {
	
	
public:
	OneColumnFile(string n) : inName(n), counter(0), outName(n + ".temp") {};
	~OneColumnFile();
	void output(int, vector<double>);
	void initFile(string);
	void resetFile();
	string getFileName()	{ return inName;	};
	
	void initFile(string, vector<string>) {};
	void output(vector<double>) {};


private:
	string outName;
	ifstream inFile;
	string inName;
	ofstream outFile;
	int counter;
};

/***********************************************************************/
class SharedOneColumnFile : public FileOutput {
	
	
public:
	SharedOneColumnFile(string n) : inName(n), counter(0), outName(n + ".temp") {};
	~SharedOneColumnFile();
	void output(int, vector<double>);
	void initFile(string);
	void resetFile();
	string getFileName()	{ return inName;	};
	
	void initFile(string, vector<string>) {};
	void output(vector<double>) {};


private:
	string outName;
	ifstream inFile;
	string inName;
	ofstream outFile;
	int counter;
		
};

/***********************************************************************/

class SharedThreeColumnFile : public FileOutput {
	
public:
	SharedThreeColumnFile(string n, string groups) : FileOutput(), groupLabel(groups), inName(n), counter(0), numGroup(1), outName(n + ".temp") {	};
	~SharedThreeColumnFile();
	void initFile(string);
	void output(int, vector<double>);
	void resetFile();
	string getFileName()	{ return inName;	};
	
	
	void initFile(string, vector<string>) {};
	void output(vector<double>) {};

private:
	string inName, groupLabel;
	string outName;
	ifstream inFile;
	ofstream outFile;
	int counter, numGroup;
};

/***********************************************************************/
//used by parsimony, unifrac.weighted and unifrac.unweighted
class ColumnFile : public FileOutput {
	
public:
	ColumnFile(string n, string i) : FileOutput(), iters(i), inName(n), counter(0), outName(n + ".temp") {};
	~ColumnFile();
	
	//to make compatible with parent class
	void output(int, vector<double>){};
	void initFile(string){};
	
	void initFile(string, vector<string>);
	void output(vector<double>);
	void resetFile();
	string getFileName()	{ return inName;	};
private:
	string inName;
	string outName;
	ifstream inFile;
	ofstream outFile;
	int counter;
	string iters;
};



#endif
