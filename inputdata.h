#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "ordervector.hpp"
#include "listvector.hpp"


using namespace std;

class InputData {
	
public:
	InputData(string, string);
	InputData(string, string, string);
	~InputData();
	ListVector* getListVector();
	OrderVector* getOrderVector();
	SAbundVector* getSAbundVector();
	
private:
	string format;
	ifstream fileHandle;
	DataVector* input;
	ListVector* list;
	OrderVector* output;
	SAbundVector* sabund;
	map<string,int> orderMap;
};


#endif