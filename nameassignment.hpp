#ifndef NAMEASSIGNMENT_HPP
#define NAMEASSIGNMENT_HPP

#include "mothur.h"
#include "listvector.hpp"

class NameAssignment : public map<string,int> {
public:
	NameAssignment(string);
	NameAssignment(){};
	void readMap();
	ListVector getListVector();
	int get(string);
	string get(int);
	void print(ostream&);
	void push_back(string);
private:
	ifstream fileHandle;
	ListVector list;
	map<int, string> reverse;
	MothurOut* m;
};




#endif
