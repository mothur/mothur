#ifndef LIST_H
#define LIST_H

#include "datavector.hpp"

/*	DataStructure for a list file.
	This class is a child to datavector.  It represents OTU information at a certain distance. 
	A list vector can be converted into and ordervector, rabundvector or sabundvector.
	Each member of the internal container "data" represents an individual OTU.
	So data[0] = "a,b,c,d,e,f".
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m  
			 rabundvector	=	6				3			2		1		1
			 sabundvector	=	2		1		1		0		0		1
			 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */

class ListVector : public DataVector {
    
public:
	ListVector();
	ListVector(int);
	ListVector(string, vector<string>, string&);
    ListVector(const ListVector& lv) : DataVector(lv.label), data(lv.data), maxRank(lv.maxRank), numBins(lv.numBins), numSeqs(lv.numSeqs), binLabels(lv.binLabels), otuTag(lv.otuTag), printListHeaders(lv.printListHeaders) {};
	ListVector(ifstream&, string&, string&);
	~ListVector(){};
	
	int getNumBins()							{	return numBins;		}
	int getNumSeqs()							{	return numSeqs;		}
	int getMaxRank()							{	return maxRank;		}

	void set(int, string);	
	string get(int);
    vector<string> getLabels();
    void setLabels(vector<string>);
    bool getPrintedLabels();
    void setPrintedLabels(bool pl) { printListHeaders = pl; }
    
	void push_back(string);
	void resize(int);
	void clear();
	int size();
    void print(ostream&);
    void print(ostream&, bool);
    void print(ostream&, map<string, int>&);
    
    RAbundVector getRAbundVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*);
    
private:
    vector<string> data;  //data[i] is a list of names of sequences in the ith OTU.
    int maxRank;
    int numBins;
    int numSeqs;
    vector<string> binLabels;
    string otuTag;
    bool printListHeaders;
    void printHeaders(ostream&, map<string, int>&, bool);
    
};

#endif
