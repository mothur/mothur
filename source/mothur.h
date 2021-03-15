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
    #include <boost/filesystem.hpp>
#endif

#ifdef USE_HDF5
    #include "H5Cpp.h"
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
#include <random>
#include <chrono>

//misc
#include <cerrno>
#include <ctime>
#include <limits>
#include <thread>
#include <mutex>

/*GSL includes*/
#ifdef USE_GSL
    #include <gsl/gsl_vector.h>
    #include <gsl/gsl_matrix.h>
    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>
    #include <gsl/gsl_math.h>
    #include <gsl/gsl_sf.h>
    #include <gsl/gsl_integration.h>
    #include <gsl/gsl_errno.h>
    #include <gsl/gsl_roots.h>
    #include <gsl/gsl_statistics_double.h>
    #include <gsl/gsl_fft_complex.h>
    #include <gsl/gsl_complex_math.h>
    #include <gsl/gsl_multimin.h>
#endif
/***********************************************************************/

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
	#include <sys/wait.h>
	#include <sys/time.h>
	#include <sys/resource.h>
	#include <sys/types.h>
	#include <sys/stat.h>
    #include <sys/sysctl.h>
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
	#include <tchar.h>

#endif

using namespace std;

#define exp(x) (exp((double) x))
#define sqrt(x) (sqrt((double) x))
#define log10(x) (log10((double) x))
#define log2(x) (log10(x)/log10(2))
#define isnan(x) ((x) != (x))
#define isinf(x) (fabs(x) == std::numeric_limits<double>::infinity())
#define GIG 1073741824
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#define PATH_SEPARATOR "/"
#define EXECUTABLE_EXT ""
#define NON_WINDOWS

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

#undef WINDOWS

#else
#define PATH_SEPARATOR "\\"
#define EXECUTABLE_EXT ".exe"
#define WINDOWS
#undef NON_WINDOWS
#define M_PI 3.14159265358979323846264338327950288

#endif

#define MOTHURMAX 1e6

typedef unsigned long ull;
typedef unsigned short intDist;
const vector<string> nullVector; //used to pass blank vector
const vector<int> nullIntVector; //used to pass blank ints
const vector<char> nullCharVector; //used to pass blank char
const map<int, int> nullIntMap;

/**************************************************************************************************/
struct classifierOTU {
    vector<vector<char> > otuData; //otuData[0] -> vector of first characters from each sequence in the OTU, otuData[1] -> vector of second characters from each sequence in the OTU
    int numSeqs;
    
    /*
     otuData.size = num columns in seq's alignment
     otuData[i].size() = numSeqs in otu
     
     seq1 > atgcaag
     seq2 > gacctga
     seq3 > cctgacg
     
     otuData[0] > {a,g,c}
     otuData[1] > {t,a,c}
     otuData[2] > {g,c,t}
     
     otuData[i] > {charInAllCols} if all chars in otuData[i] are identical. ie, ignore column
     otuData[i] > {a} all seqs contain 'a' in column i of alignment
     */
    classifierOTU(){ numSeqs = 0; }
    classifierOTU(string aligned) {
        for (int i = 0; i < aligned.length(); i++) {
            vector<char> thisSpot;
            thisSpot.push_back(aligned[i]);
            
            otuData.push_back(thisSpot);
        }
        numSeqs = 1;
    }
    classifierOTU(vector<string> otu) { readSeqs(otu); }
    classifierOTU(vector<vector<char> > otu, int num) : otuData(otu), numSeqs(num) {}
    
    void readSeqs(vector<vector<char> > otu, int num) { otuData = otu; numSeqs = num; } //for shortcut files
    
    void readSeqs(vector<string> otu) {
        int alignedLength = 0;
        bool error = false;
        if (otu.size() != 0) { alignedLength = otu[0].length(); }
        for (int j = 0; j < otu.size(); j++) { if (otu[j].length() != alignedLength) { error = true;} }
        
        if (!error) {
            
            for (int i = 0; i < alignedLength; i++) {
                vector<char> thisSpot; set<char> thisChars;
                for (int j = 0; j < otu.size(); j++) {
                    thisSpot.push_back(otu[j][i]);
                    thisChars.insert(otu[j][i]);
                }
                if (thisChars.size() == 1) { thisSpot.clear(); thisSpot.push_back(*thisChars.begin()); }// all same, reduce to 1.
                otuData.push_back(thisSpot);
            }
            numSeqs = otu.size();
        }else {  numSeqs = 0; }
    }
};
//******************************************************
struct mcmcSample {
    double alpha, beta; //dmDash, dV
    double dNu;
    int ns;
    
    mcmcSample() {}
    mcmcSample(double a, double b, double d, int n) : alpha(a), beta(b), dNu(d), ns(n) {}
    
};
typedef struct s_Params
{
    long lSeed;
    
    string szOutFileStub;
    
    double dSigmaX; //dSigmaM, dSigmaA
    
    double dSigmaY;  //dSigmaV, dSigmaB
    
    double dSigmaN;
    
    double dSigmaS;
    
    int nIter;
} t_Params;


typedef struct s_Data
{
    int nNA;
    
    int **aanAbund;
    
    int nL;
    
    int nJ;
}t_Data;

typedef struct s_LNParams
{
    double dMDash;
    
    double dV;
    
    int    n;
} t_LNParams;

typedef struct s_LSParams
{
    double dMDash;
    
    double dV;
    
    double dNu;
    
    double dC;
    
    int n;
    
} t_LSParams;

typedef struct s_IGParams //s_SIParams
{
    int    nS;      /*number of species in community*/
    
    double dAlpha;
    
    double dBeta;
    
    double dC; //dGamma
    
    int n;
    
} t_IGParams;

#ifdef USE_GSL
typedef struct s_MetroInit
{
    t_Params *ptParams;
    
    t_Data   *ptData;
    
    gsl_vector* ptX;
    
    int nAccepted;
    
    long lSeed;
    
    int nThread;
    
} t_MetroInit;

#endif

//***********************************************************************
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

struct countTableItem {
    int abund;
    int group;
    
    countTableItem() { abund = 0; group = -1; }
    countTableItem(int a, int g) : abund(a), group(g) {}
    ~countTableItem() {}
};

struct item {
    int name;
    int group;
    
    item() { name = -1; group = -1; }
    item(int n, int g) : name(n), group(g) {}
    ~item() {}
};


struct weightedSeq {
    long long name;
    long long weight;
    weightedSeq(long long n, long long w)  { name = n; weight = w;     }
};

struct PCell{
    ull row;
    ull column;
    float dist;
    PCell** vectorMap;
    PCell() : row(0), column(0), dist(0), vectorMap(NULL) {};
    PCell(ull r, ull c, float d) : row(r), column(c), dist(d), vectorMap(NULL) {};
};

/* For each distance in a sparse matrix we have a row, column and distance.
 The PDistCell consists of the column and distance.
 We know the row by the distances row in the seqVec matrix.
 SeqVec is square and each row is sorted so the column values are ascending to save time in the search for the smallest distance. */

/***********************************************************************/
struct PDistCellMin{
    ull row;
    ull col;
    
    PDistCellMin(ull r, ull c) :  col(c), row(r) {}
};
/***********************************************************************/

struct colDist {
    int col;
    int row;
    float dist;
    colDist(int r, int c, double d) : row(r), col(c), dist(d) {}
};
/************************************************************/
struct seqPNode {
    int numIdentical;
    string name;
    string sequence;
    vector<int> clusteredIndexes; //indexes of merge nodes. Can use this later to construct names
    int diffs;
    
    seqPNode() { diffs = 0; numIdentical = 0; name = ""; sequence = "";  }
    seqPNode(string na, string seq, int n, vector<int> nm) : numIdentical(n), name(na), sequence(seq), clusteredIndexes(nm) { diffs = 0; }
    ~seqPNode() {}
};

/**********************************************************/
struct linePair {
    double start;
    double end;
    linePair(double i, double j) : start(i), end(j) {}
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
    string label;
    listCt() :  bin(""), binSize(0), label("") {};
    listCt(string b, int a, string l) :  bin(b), binSize(a), label(l) {}
};
/***********************************************************************/
struct consTax2{
    string otuName;
    string taxonomy;
    int abundance;
	consTax2() :  otuName("OTUxxx"), taxonomy("unknown"), abundance(0) {};
	consTax2(string n, string t, int a) :  otuName(n), taxonomy(t), abundance(a) {}
};
/***********************************************************************/
struct Taxon {
    string name;
    float confidence;

    Taxon(string n, float conf) : name(n), confidence(conf) {}
    ~Taxon(){}
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
inline bool compareGroups(countTableItem left, countTableItem right){
	return (left.group > right.group);
}
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
/********************************************************************************************************************/
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
//*******************************************************************************
template <typename T> int sgn(T val){ return (val > T(0)) - (val < T(0)); }
//*******************************************************************************
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
//***********************************************************************
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
//**********************************************************************
template<typename T>
string toString(const T&x){
	
		stringstream output;
		output << x;
		return output.str();
	
}
//***********************************************************************
template<typename T>
void mothurSwap(T&x, T&y){
    T temp = y;
    y = x;
    x = temp;
}
//*************************************************************************
template<typename T>
string toHex(const T&x){
	
		stringstream output;
		
		output << hex << x;

		return output.str();
	
}
//*********************************************************************
template<typename T>
string toString(const T&x, int i){
	
		stringstream output;
		
		output.precision(i);
		output << fixed << x;
		
		return output.str();
}
//***********************************************************************
template<class T>
T fromString(const string& s){
	istringstream stream (s);
	T t;
	stream >> t;
	return t;
}
//****************************************************************************

#endif
