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

//misc
#include <cerrno>
#include <ctime>
#include <limits>


#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
	#include <sys/wait.h>
	#include <unistd.h>
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
};

struct ThreadNode {
	int* pid;
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
inline bool convertTestFloat(const string& s, T& x, bool failIfLeftoverChars = true){
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
inline bool convertTest(const string& s, T& x, bool failIfLeftoverChars = true){
	istringstream i(s);
	char c;
	if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
	{
		cout << "unable to be converted into an integer.\n" << endl;
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

inline int openOutputFileAppend(string fileName, ofstream& fileHandle){
	
	fileHandle.open(fileName.c_str(), ios::app);
	if(!fileHandle) {
		cout << "Error: Could not open " <<  fileName << endl; 
		return 1;
	}
	else {
		return 0;
	}

}
/***********************************************************************/

inline string getline(ifstream& fileHandle) {
	try {
	
		string line = "";
		
		while (!fileHandle.eof())	{
			//get next character
			char c = fileHandle.get(); 
			
			//are you at the end of the line
			if ((c == '\n') || (c == '\r') || (c == '\f')){   break;	}	
			else {		line += c;		}
		}
		
		return line;
		
	}
	catch(exception& e) {
		cout << "Error in mothur function getline" << endl;
		exit(1);
	}
}

/**************************************************************************************************/

inline void mothurOut(string message) {
	try{
		ofstream out;
		string logFileName = "mothur.logFile";
		openOutputFileAppend(logFileName, out);
		
		cout << message;
		out << message;
		
		out.close();
	}
	catch(exception& e) {
		cout << "Error in mothur class mothurOut" << endl;
		exit(1);
	}
}
/**************************************************************************************************/

inline void mothurOut(string message, string precision) {
	try{
		ofstream out;
		string logFileName = "mothur.logFile";
		openOutputFileAppend(logFileName, out);
		
		cout << precision << message;
		out << precision << message;
		
		out.close();
	}
	catch(exception& e) {
		cout << "Error in mothur class mothurOut" << endl;
		exit(1);
	}
}

/**************************************************************************************************/

inline void mothurOutEndLine() {
	try {
		ofstream out;
		string logFileName = "mothur.logFile";
		openOutputFileAppend(logFileName, out);
		
		cout << endl;  
		out << endl;
		
		out.close();
	}
	catch(exception& e) {
		cout << "error in mothur mothurOutEndLine" << endl;
		exit(1);
	}
}


/**************************************************************************************************/

inline void errorOut(exception& e, string object, string function) {
	
		mothurOut("Error: ");
		mothurOut(toString(e.what()));
		mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at pschloss@microbio.umass.edu, and be sure to include the mothur.logFile with your inquiry.");
		mothurOutEndLine();
	
}




/***********************************************************************/

inline void gobble(istream& f){
	
	char d;
    while(isspace(d=f.get()))		{;}
	f.putback(d);
	
}
/***********************************************************************/

inline bool isTrue(string f){
	
	if ((f == "TRUE") || (f == "T")	|| (f == "true") || (f == "t")) {	return true;	}
	else {	return false;  }
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
	
	if(longName.find_last_of('/') != longName.npos){
		int pos = longName.find_last_of('/')+1;
		rootPathName = longName.substr(0, pos);
	}

	return rootPathName;
}

/***********************************************************************/

inline string getExtension(string longName){
	
	string extension = longName;
	
	if(longName.find_last_of('.') != longName.npos){
		int pos = longName.find_last_of('.');
		extension = longName.substr(pos, longName.length());
	}
	
	return extension;
}

/***********************************************************************/

inline int openInputFile(string fileName, ifstream& fileHandle){

	fileHandle.open(fileName.c_str());
	if(!fileHandle) {
		mothurOut("Error: Could not open " + fileName);  mothurOutEndLine();
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
		mothurOut("Error: Could not open " + fileName);  mothurOutEndLine();
		return 1;
	}
	else {
		return 0;
	}

}

/***********************************************************************/

inline int getNumSeqs(ifstream& file){
	
	int numSeqs = count(istreambuf_iterator<char>(file),istreambuf_iterator<char>(), '>');
	file.seekg(0);
	return numSeqs;

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
		errorOut(e, "mothur", "splitAtDash");
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
		errorOut(e, "mothur", "splitAtDash");
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
		errorOut(e, "mothur", "splitAtDash");
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
		errorOut(e, "mothur", "splitAtComma");
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
		errorOut(e, "mothur", "splitAtComma");
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
		errorOut(e, "mothur", "splitAtEquals");
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
		errorOut(e, "mothur", "inUsersGroups");
		exit(1);
	}
}

/**************************************************************************************************/

inline void mothurOutJustToLog(string message) {
	try {
		ofstream out;
		string logFileName = "mothur.logFile";
		openOutputFileAppend(logFileName, out);
		
		out << message;
		
		out.close();
	}
	catch(exception& e) {
		errorOut(e, "mothur", "mothurOutJustToLog");
		exit(1);
	}
}


/**************************************************************************************************/

inline void mothurOut(float num) {
	try {
		ofstream out;
		string logFileName = "mothur.logFile";
		openOutputFileAppend(logFileName, out);
		
		cout << num;  
		out << num;
		
		out.close();
	}
	catch(exception& e) {
		cout << "Error in mothur class mothurOut float" << endl;
		exit(1);
	}
}
/***********************************************************************/
inline void mothurOut(double value) {
	try {
		ofstream out;
		string logFileName = "mothur.logFile";
		openOutputFileAppend(logFileName, out);
		
		cout << value;  
		out << value;
		
		out.close();
	}
	catch(exception& e) {
		cout << "Error in mothur class mothurOut double" << endl;
		exit(1);
	}
}

/***********************************************************************/
//this function determines if the user has given us labels that are smaller than the given label.
//if so then it returns true so that the calling function can run the previous valid distance.
//it's a "smart" distance function.  It also checks for invalid labels.
inline bool anyLabelsToProcess(string label, set<string>& userLabels, string errorOff) {
	try {
		set<string>::iterator it;
		vector<float> orderFloat;
		map<string, float> userMap;  //the conversion process removes trailing 0's which we need to put back
		map<string, float>::iterator it2;
		float labelFloat;
		bool smaller = false;
		
		//unique is the smallest line
		if (label == "unique") {  return false;  }
		else { convert(label, labelFloat); }
		
		//go through users set and make them floats
		for(it = userLabels.begin(); it != userLabels.end(); ++it) {
			
			float temp;
			if ((*it != "unique") && (convertTestFloat(*it, temp) == true)){
				convert(*it, temp);
				orderFloat.push_back(temp);
				userMap[*it] = temp;
			}else if (*it == "unique") { 
				orderFloat.push_back(-1.0);
				userMap["unique"] = -1.0;
			}else {
				if (errorOff == "") {  mothurOut(*it + " is not a valid label."); mothurOutEndLine();  }
				userLabels.erase(*it); 
				it--;
			}
		}
		
		//sort order
		sort(orderFloat.begin(), orderFloat.end());
		
		/*************************************************/
		//is this label bigger than any of the users labels
		/*************************************************/
				
		//loop through order until you find a label greater than label
		for (int i = 0; i < orderFloat.size(); i++) {
			if (orderFloat[i] < labelFloat) {
				smaller = true;
				if (orderFloat[i] == -1) { 
					if (errorOff == "") { mothurOut("Your file does not include the label unique."); mothurOutEndLine(); }
					userLabels.erase("unique");
				}
				else {  
					if (errorOff == "") { mothurOut("Your file does not include the label "); mothurOutEndLine(); }
					string s = "";
					for (it2 = userMap.begin(); it2!= userMap.end(); it2++) {  
						if (it2->second == orderFloat[i]) {  
							s = it2->first;  
							//remove small labels
							userLabels.erase(s);
							break;
						}
					}
					if (errorOff == "") { mothurOut(s + ". I will use the next smallest distance. "); mothurOutEndLine(); }
				}
			//since they are sorted once you find a bigger one stop looking
			}else { break; }
		}
		
		return smaller;
						
	}
	catch(exception& e) {
		errorOut(e, "mothur", "anyLabelsToProcess");
		exit(1);
	}
}

/**************************************************************************************************/
#endif

