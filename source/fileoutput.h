#ifndef FILEOUTPUT_H
#define FILEOUTPUT_H

#include "mothurout.h"
#include "utils.hpp"

/***********************************************************************/

class FileOutput {
	
public:
	FileOutput(string n){ m = MothurOut::getInstance(); fileHeader = ""; filename = n; firstLabel = true; }
    virtual ~FileOutput(){ printFile(); }
	
    virtual void setLabelName(string) {}
    virtual void updateOutput(int, vector<double>) {}
    virtual void resetFile() { firstLabel = false;  }
    
    virtual void setLabelName(string, vector<string>) {}
    virtual void updateOutput(vector<double>) {}

protected:
	MothurOut* m;
    Utils util;
    string filename, fileHeader;
    bool firstLabel;
    vector< vector<double> > results;
    
    void printFile();
};	
	
/***********************************************************************/

class ThreeColumnFile : public FileOutput {
	
public:
    ThreeColumnFile(string n) : FileOutput(n) { }
    ~ThreeColumnFile() {}
    
	void setLabelName(string);
	void updateOutput(int, vector<double>);

private:
    
	
};

/***********************************************************************/
class OneColumnFile : public FileOutput {
	
	
public:
	OneColumnFile(string n) : FileOutput(n) { }
	~OneColumnFile() {}
    
    void setLabelName(string);
    void updateOutput(int, vector<double>);
	
private:
	
};

/***********************************************************************/
class SharedOneColumnFile : public FileOutput {
	
	
public:
	SharedOneColumnFile(string n) : FileOutput(n) {}
	~SharedOneColumnFile() {}
	
    void setLabelName(string);
    void updateOutput(int, vector<double>);
	
private:	
		
};

/***********************************************************************/

class SharedThreeColumnFile : public FileOutput {
	
public:
    SharedThreeColumnFile(string n, string groups) : FileOutput(n), groupLabel(groups), numGroup(1) { }
    ~SharedThreeColumnFile() {}
    
    void setLabelName(string);
    void updateOutput(int, vector<double>);
		
private:
	string groupLabel;
	int numGroup;
};

/***********************************************************************/
//used by parsimony, unifrac.weighted and unifrac.unweighted
class ColumnFile : public FileOutput {
	
public:
	ColumnFile(string n, string i) : FileOutput(n) {}
    ~ColumnFile() {}
		
	void setLabelName(string, vector<string>);
	void updateOutput(vector<double>);
    
private:
	
    
};
/***********************************************************************/


#endif
