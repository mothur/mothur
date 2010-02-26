#ifndef INPUTDATA_H
#define INPUTDATA_H

#include "mothur.h"
#include "ordervector.hpp"
#include "sharedlistvector.h"
#include "sharedordervector.h"
#include "listvector.hpp"


class InputData {
	
public:
	InputData(string, string);
	InputData(string, string, string);
	~InputData();
	ListVector* getListVector();
	ListVector* getListVector(string);  //pass the label you want
	SharedListVector* getSharedListVector();
	SharedListVector* getSharedListVector(string);  //pass the label you want
	OrderVector* getOrderVector(); 
	OrderVector* getOrderVector(string); //pass the label you want
	SharedOrderVector* getSharedOrderVector();
	SharedOrderVector* getSharedOrderVector(string);  //pass the label you want
	SAbundVector* getSAbundVector();
	SAbundVector* getSAbundVector(string);  //pass the label you want
	RAbundVector* getRAbundVector();
	RAbundVector* getRAbundVector(string);  //pass the label you want
	vector<SharedRAbundVector*> getSharedRAbundVectors();
	vector<SharedRAbundVector*> getSharedRAbundVectors(string);  //pass the label you want
	
private:
	string format;
	ifstream fileHandle;
	DataVector* input;
	ListVector* list;
	SharedListVector* SharedList;
	OrderVector* output;
	SharedOrderVector* SharedOrder;
	SAbundVector* sabund;
	RAbundVector* rabund;
	map<string,int> orderMap;
	string filename;
	MothurOut* m;
};


#endif
