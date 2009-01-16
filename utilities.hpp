#ifndef UTILITIES_H
#define UTILITIES_H

using namespace std;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <stdexcept>

typedef unsigned long long ull;

/***********************************************************************/

// snagged from http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2
// works for now, but there should be a way to do it without killing the whole program

class BadConversion : public runtime_error {
public:
	BadConversion(const string& s) : runtime_error(s){ }
};

//**********************************************************************************************************************

template<typename T>
inline void convert(const string& s, T& x, bool failIfLeftoverChars = true){
	istringstream i(s);
	char c;
	if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
		throw BadConversion(s);
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
string toString(const T&x, int i){
	stringstream output;
	
	output.precision(i);
    output << fixed << x;
	
    return output.str();
}


/***********************************************************************/

inline void gobble(ifstream& f){
	
	char d;
    while(isspace(d=f.get()))		{;}
	f.putback(d);
	
}

/***********************************************************************/

inline float roundDist(float dist, int precision){
	
	return int(dist * precision + 0.5)/float(precision);
	
}

/***********************************************************************/

inline int getNumNames(string names){
	
	int count = 0;
	
	if(names != ""){
		count = 1;
		for(int i=0;i<names.size();i++){
			if(names[i] == ','){
				count++;
			}
		}
	}
	
	return count;
	
}

/**************************************************************************************************/

inline vector<vector<double> > binomial(int maxOrder){
	
	vector<vector<double> > binomial(maxOrder+1);
	
    for(int i=0;i<=maxOrder;i++){
		binomial[i].resize(maxOrder+1);
		binomial[i][0]=1;
		binomial[0][i]=0;
    }
    binomial[0][0]=1;
	
    binomial[1][0]=1;
    binomial[1][1]=1;
	
    for(int i=2;i<=maxOrder;i++){
		binomial[1][i]=0;
    }
	
    for(int i=2;i<=maxOrder;i++){
		for(int j=1;j<=maxOrder;j++){
			if(i==j){	binomial[i][j]=1;									}
			if(j>i)	{	binomial[i][j]=0;									}
			else	{	binomial[i][j]=binomial[i-1][j-1]+binomial[i-1][j];	}
		}
    }
	
	return binomial;
}

/***********************************************************************/

inline string getRootName(string longName){
 
	string rootName = longName;
	
	if(longName.find_last_of(".") != longName.npos){
		int pos = longName.find_last_of('.')+1;
		rootName = longName.substr(0, pos);
	}

	return rootName;
}
/***********************************************************************/

inline string getSimpleName(string longName){
 
	string simpleName = longName;
	
	if(longName.find_last_of("/") != longName.npos){
		int pos = longName.find_last_of('/')+1;
		simpleName = longName.substr(pos, longName.length());
	}

	return simpleName;
}

/***********************************************************************/

inline string getPathName(string longName){
 
	string rootPathName = longName;
	
	if(longName.find_last_of("/") != longName.npos){
		int pos = longName.find_last_of('/')+1;
		rootPathName = longName.substr(0, pos);
	}

	return rootPathName;
}

/***********************************************************************/

inline int openInputFile(string fileName, ifstream& fileHandle){

	fileHandle.open(fileName.c_str());
	if(!fileHandle) {
		cerr << "Error: Could not open " << fileName << endl;
		return 1;
	}
	else {
		return 0;
	}
	
}

/***********************************************************************/

inline int openOutputFile(string fileName, ofstream& fileHandle){
	
	fileHandle.open(fileName.c_str(), ios::trunc);
	if(!fileHandle) {
		cerr << "Error: Could not open " << fileName << endl;
		return 1;
	}
	else {
		return 0;
	}

}

/***********************************************************************/

#endif
