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

//misc
#include <cerrno>
#include <ctime>
#include <limits>

#ifdef USE_MPI
	#include "mpi.h"
#endif
/***********************************************************************/

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
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
	float dist;
	seqDist() {}
	seqDist(int s1, int s2, float d) : seq1(s1), seq2(s2), dist(d) {}
	~seqDist() {}
};
/************************************************************/
struct distlinePair {
	int start;
	int end;
	
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
/***************************************************************/
struct spearmanRank {
	string name;
	float score;
	
	spearmanRank(string n, float s) : name(n), score(s) {}
};
//********************************************************************************************************************
//sorts highest to lowest
inline bool compareSpearman(spearmanRank left, spearmanRank right){
	return (left.score > right.score);	
} 
//********************************************************************************************************************
//sorts highest to lowest
inline bool compareSeqPriorityNodes(seqPriorityNode left, seqPriorityNode right){
	return (left.numIdentical > right.numIdentical);	
} 
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSpearmanReverse(spearmanRank left, spearmanRank right){
	return (left.score < right.score);	
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

