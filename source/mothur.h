#ifndef MOTHUR_H
#define MOTHUR_H



/*
 *  mothur.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/19/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This file contains all the standard incudes we use in the project as well as some common utilities. */

//#include <cstddef>

//boost libraries
#ifdef USE_BOOST
    #include <boost/iostreams/filtering_stream.hpp>
    #include <boost/iostreams/filter/gzip.hpp>
#endif

//io libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <signal.h>


//exception
#include <stdexcept>
#include <exception>
#include <cstdlib> 


//containers
#include <vector>
#include <set>
#include <map>
#include <string>
#include <list>
#include <string.h>

//math
#include <cmath>
#include <math.h>
#include <algorithm>
#include <numeric>

//misc
#include <cerrno>
#include <ctime>
#include <limits>

#ifdef USE_MPI
	#include "mpi.h"
#endif
/***********************************************************************/

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
	#include <sys/wait.h>
	#include <sys/time.h>
	#include <sys/resource.h>
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
	
	#ifdef USE_READLINE
		#include <readline/readline.h>
		#include <readline/history.h>
	#endif

#else
	#include <conio.h> //allows unbuffered screen capture from stdin
	#include <direct.h> //get cwd
	#include <windows.h>
	#include <psapi.h>
	#include <direct.h>
	#include <tchar.h>

#endif

using namespace std;

#define exp(x) (exp((double) x))
#define sqrt(x) (sqrt((double) x))
#define log10(x) (log10((double) x))
#define log2(x) (log10(x)/log10(2))
#define isnan(x) ((x) != (x))
#define isinf(x) (fabs(x) == std::numeric_limits<double>::infinity())


typedef unsigned long ull;
typedef unsigned short intDist;

struct IntNode {
	int lvalue;
	int rvalue;
	int lcoef;
	int rcoef;
	IntNode* left;
	IntNode* right;
	
	IntNode(int lv, int rv, IntNode* l, IntNode* r) : lvalue(lv), rvalue(rv), left(l), right(r) {};
	IntNode() {};
};

struct ThreadNode {
	int* pid;
	IntNode* left;
	IntNode* right;
};

struct diffPair {
	float	prob;
	float	reverseProb;
	
	diffPair() {
		prob = 0; reverseProb = 0;
	}
	diffPair(float p, float rp) {
		prob = p;
		reverseProb = rp;
	}
};

/**********************************************************/
struct CommonHeader {
	unsigned int magicNumber;
	string version;
	unsigned long long indexOffset;
	unsigned int indexLength;
	unsigned int numReads;
	unsigned short headerLength;
	unsigned short keyLength;
	unsigned short numFlowsPerRead;
	int flogramFormatCode;
	string flowChars; //length depends on number flow reads
	string keySequence; //length depends on key length
	
	CommonHeader(){ magicNumber=0; indexOffset=0; indexLength=0; numReads=0; headerLength=0; keyLength=0; numFlowsPerRead=0; flogramFormatCode='s'; }
	~CommonHeader() { }
};
/**********************************************************/
struct Header {
	unsigned short headerLength;
	unsigned short nameLength;
	unsigned int numBases;
	unsigned short clipQualLeft;
	unsigned short clipQualRight;
	unsigned short clipAdapterLeft;
	unsigned short clipAdapterRight;
	string name; //length depends on nameLength
	string timestamp;
	string region;
	string xy;
	
	Header() { headerLength=0; nameLength=0; numBases=0; clipQualLeft=0; clipQualRight=0; clipAdapterLeft=0; clipAdapterRight=0; }
	~Header() { }
};
/**********************************************************/
struct seqRead {
	vector<unsigned short> flowgram;
	vector<unsigned int> flowIndex;
	string bases;
	vector<unsigned int> qualScores;
	
	seqRead() { }
	~seqRead() { }
};
/**********************************************************/
struct linePair {
    unsigned long long start;
    unsigned long long end;
    linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
    linePair(){ start=0; end=0; }
    ~linePair(){}
};
/***********************************************************************/
struct PDistCell{
	ull index;
	float dist;
	PDistCell() :  index(0), dist(0) {};
	PDistCell(ull c, float d) :  index(c), dist(d) {}
};
/***********************************************************************/
struct consTax{
	string name;
    string taxonomy;
    int abundance;
	consTax() :  name(""), taxonomy("unknown"), abundance(0) {};
	consTax(string n, string t, int a) :  name(n), taxonomy(t), abundance(a) {}
};
/***********************************************************************/
struct listCt{
    string bin;
    int binSize;
    listCt() :  bin(""), binSize(0) {};
    listCt(string b, int a) :  bin(b), binSize(a) {}
};
/***********************************************************************/
struct consTax2{
    string taxonomy;
    int abundance;
    string otuName;
	consTax2() :  otuName("OTUxxx"), taxonomy("unknown"), abundance(0) {};
	consTax2(string n, string t, int a) :  otuName(n), taxonomy(t), abundance(a) {}
};
/************************************************************/
struct clusterNode {
	int numSeq;
	int parent;
	int smallChild; //used to make linkTable work with list and rabund. represents bin number of this cluster node
	clusterNode(int num, int par, int kid) : numSeq(num), parent(par), smallChild(kid) {};
};
/************************************************************/
struct seqDist {
	int seq1;
	int seq2;
	double dist;
	seqDist() {}
	seqDist(int s1, int s2, double d) : seq1(s1), seq2(s2), dist(d) {}
	~seqDist() {}
};
/************************************************************/
struct distlinePair {
	int start;
	int end;
	
};
/************************************************************/
struct oligosPair {
	string forward;
	string reverse;
	
	oligosPair() { forward = ""; reverse = "";  }
	oligosPair(string f, string r) : forward(f), reverse(r) {}
	~oligosPair() {}
};

/************************************************************/
struct seqPriorityNode {
	int numIdentical;
	string seq;
	string name;
	seqPriorityNode() {}
	seqPriorityNode(int n, string s, string nm) : numIdentical(n), seq(s), name(nm) {}
	~seqPriorityNode() {}
};
/************************************************************/
struct compGroup {
	string group1;
	string group2;
	compGroup() {}
	compGroup(string s, string nm) : group1(s), group2(nm) {}
    string getCombo() { return group1+"-"+group2; }
	~compGroup() {}
};
/***************************************************************/
struct spearmanRank {
	string name;
	float score;
	
	spearmanRank(string n, float s) : name(n), score(s) {}
};
//***********************************************************************
inline bool compareIndexes(PDistCell left, PDistCell right){
	return (left.index > right.index);	
}
//********************************************************************************************************************
inline bool compareSpearman(spearmanRank left, spearmanRank right){
	return (left.score < right.score);	
}
//********************************************************************************************************************
inline double max(double left, double right){
    if (left > right) { return left; }
    else { return right; }
}
//********************************************************************************************************************
inline double max(int left, double right){
    double value = left;
    if (left > right) { return value; }
    else { return right; }
}
//********************************************************************************************************************
inline double max(double left, int right){
    double value = right;
    if (left > value) { return left; }
    else { return value; }
}
//********************************************************************************************************************
//sorts highest to lowest
inline bool compareSeqPriorityNodes(seqPriorityNode left, seqPriorityNode right){
	if (left.numIdentical > right.numIdentical) {
        return true;
    }else if (left.numIdentical == right.numIdentical) {
        if (left.seq > right.seq) { return true; }
        else { return false; }
    }
    return false;	
} 
 
/************************************************************/
//sorts lowest to highest
inline bool compareDistLinePairs(distlinePair left, distlinePair right){
	return (left.end < right.end);	
} 
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSequenceDistance(seqDist left, seqDist right){
	return (left.dist < right.dist);	
}
//********************************************************************************************************************
//returns sign of double
inline double sign(double temp){
	//find sign
    if (temp > 0)       { return 1.0;   }
    else if (temp < 0)  { return -1.0;  }
    return 0;
}
/***********************************************************************/

// snagged from http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2
// works for now, but there should be a way to do it without killing the whole program

class BadConversion : public runtime_error {
public:
	BadConversion(const string& s) : runtime_error(s){ }
};

//**********************************************************************************************************************
template<typename T>
void convert(const string& s, T& x, bool failIfLeftoverChars = true){
	
		istringstream i(s);
		char c;
		if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
			throw BadConversion(s);
	
}
//**********************************************************************************************************************
template <typename T> int sgn(T val){ return (val > T(0)) - (val < T(0)); }
//**********************************************************************************************************************

template<typename T>
bool convertTestFloat(const string& s, T& x, bool failIfLeftoverChars = true){
	
		istringstream i(s);
		char c;
		if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
		{
			return false;
		} 
		return true;
	
}

//**********************************************************************************************************************

template<typename T>
bool convertTest(const string& s, T& x, bool failIfLeftoverChars = true){
	
		istringstream i(s);
		char c;
		if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
		{
			return false;
		} 
		return true;
	
}
//**********************************************************************************************************************
template<typename T>
string toString(const T&x){
	
		stringstream output;
		output << x;
		return output.str();
	
}

//**********************************************************************************************************************

template<typename T>
string toHex(const T&x){
	
		stringstream output;
		
		output << hex << x;

		return output.str();
	
}
//**********************************************************************************************************************

template<typename T>
string toString(const T&x, int i){
	
		stringstream output;
		
		output.precision(i);
		output << fixed << x;
		
		return output.str();
	
}
//**********************************************************************************************************************

template<class T>
T fromString(const string& s){
	istringstream stream (s);
	T t;
	stream >> t;
	return t;
}

//**********************************************************************************************************************

#endif

