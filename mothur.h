#ifndef MOTHUR_H
#define MOTHUR_H

using namespace std;


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

//math
#include <cmath>
#include <math.h>
#include <algorithm>

typedef unsigned long long ull;

struct IntNode {
	int lvalue;
	int rvalue;
	int lcoef;
	int rcoef;
	IntNode* left;
	IntNode* right;
};
	
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
inline bool convertTest(const string& s, T& x, bool failIfLeftoverChars = true){
	istringstream i(s);
	char c;
	if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
	{
		cout << "'" << s << "' is unable to be converted into an integer.\n";
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


/***********************************************************************/

inline void gobble(istream& f){
	
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

inline int factorial(int num){
	int total = 1;
	
	for (int i = 1; i <= num; i++) {
		total *= i;
	}
	
	return total;
}
/**************************************************************************************************

double min(double x, double y)
{
    if(x<y){	return x;    }
    else   {	return y;    }
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

//This function parses the estimator options and puts them in a vector
inline void splitAtDash(string& estim, vector<string>& container) {
	try {
		string individual;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				container.push_back(individual);
			}
		}
		//get last one
		container.push_back(estim);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/***********************************************************************/
//This function parses the label options and puts them in a set
inline void splitAtDash(string& estim, set<string>& container) {
	try {
		string individual;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				container.insert(individual);
			}
		}
		//get last one
		container.insert(estim);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/***********************************************************************/
//This function parses the line options and puts them in a set
inline void splitAtDash(string& estim, set<int>& container) {
	try {
		string individual;
		int lineNum;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				convert(individual, lineNum); //convert the string to int
				container.insert(lineNum);
			}
		}
		//get last one
		convert(estim, lineNum); //convert the string to int
		container.insert(lineNum);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function splitAtDash. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/***********************************************************************/
//This function parses the a string and puts peices in a vector
inline void splitAtComma(string& estim, vector<string>& container) {
	try {
		string individual;
		
		while (estim.find_first_of(',') != -1) {
			individual = estim.substr(0,estim.find_first_of(','));
			if ((estim.find_first_of(',')+1) <= estim.length()) { //checks to make sure you don't have comma at end of string
				estim = estim.substr(estim.find_first_of(',')+1, estim.length());
				container.push_back(individual);
			}
		}
		//get last one
		container.push_back(estim);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function splitAtComma. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function splitAtComma. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

//This function splits up the various option parameters
inline void splitAtComma(string& prefix, string& suffix){
	try {
		prefix = suffix.substr(0,suffix.find_first_of(','));
		if ((suffix.find_first_of(',')+2) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
			suffix = suffix.substr(suffix.find_first_of(',')+1, suffix.length());
			string space = " ";
			while(suffix.at(0) == ' ')
				suffix = suffix.substr(1, suffix.length());
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function splitAtComma. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function splitAtComma. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/***********************************************************************/

//This function separates the key value from the option value i.e. dist=96_...
inline void splitAtEquals(string& key, string& value){		
	try {
		if(value.find_first_of('=') != -1){
			key = value.substr(0,value.find_first_of('='));
			if ((value.find_first_of('=')+1) <= value.length()) {
				value = value.substr(value.find_first_of('=')+1, value.length());
			}
		}else{
			key = value;
			value = 1;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function splitAtEquals. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function splitAtEquals. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/**************************************************************************************************/

inline bool inUsersGroups(string groupname, vector<string> Groups) {
	try {
		for (int i = 0; i < Groups.size(); i++) {
			if (groupname == Groups[i]) { return true; }
		}
		return false;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the mothur class Function inUsersGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the mothur class function inUsersGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/**************************************************************************************************/


#endif

