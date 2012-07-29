#ifndef CALCULATOR_H
#define CALCULATOR_H


#include "mothur.h"
#include "sabundvector.hpp"
#include "sharedsabundvector.h"
#include "rabundvector.hpp"
#include "uvest.h"
#include "mothurout.h"

/* The calculator class is the parent class for all the different estimators implemented in mothur except the tree calculators.
It has 2 pure functions EstOutput getValues(SAbundVector*), which works on a single group, and 
EstOutput getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2), which compares 2 groups. */ 


typedef vector<double> EstOutput;

/***********************************************************************/

class Calculator {

public:
	Calculator(){ m = MothurOut::getInstance(); needsAll = false; }
	virtual ~Calculator(){};
	Calculator(string n, int c, bool f) : name(n), cols(c), multiple(f) { m = MothurOut::getInstance(); needsAll = false; };
	Calculator(string n, int c, bool f, bool a) : name(n), cols(c), multiple(f), needsAll(a) { m = MothurOut::getInstance(); };
	virtual EstOutput getValues(SAbundVector*) = 0;	
	virtual EstOutput getValues(vector<SharedRAbundVector*>) = 0;
	virtual void print(ostream& f)	{ f.setf(ios::fixed, ios::floatfield); f.setf(ios::showpoint);
									  f << data[0]; for(int i=1;i<data.size();i++){	f << '\t' << data[i];	}}
	virtual string getName()		{	return name;	}
	virtual int getCols()		{	return cols;	}
	virtual bool getMultiple()  {   return multiple;   }
	virtual bool getNeedsAll()  {   return needsAll;   }
	virtual string getCitation() = 0;
	void citation() { m->mothurOut(getCitation()); m->mothurOutEndLine(); }
protected:
	MothurOut* m;
	EstOutput data;
	string name;
	int cols;
	bool multiple;
	bool needsAll;

};

/**************************************************************************************************/
/*This Class holds all of the methods that manipulate vectors.
These methods are used in the other classes.
This class must be included if any of the other classes are to be used. */

class VecCalc
{
	// The methods seen in the order here is how they are ordered throughout the class.
	public:
		VecCalc(){};
		//void printElements(vector<double>); //This prints the values of the vector on one line with a space between each value.
		//void printElements(vector<string>); //This prints the values of the vector on one line with a space between each value.
		//int findString(vector<string>, string);//This returns the index of the given string in the given <string> vector, if the string does not exist in the vector it returns -1.
		//double mean(vector<double>); //This returns the mean value of the vector.
		//double stError(vector<double>); //This returns the standard error of the vector.
		int sumElements(vector<int>, int);
		int sumElements(vector<int>);
		double sumElements(vector<double>); //This returns the sum of all the values in the vector.
		double sumElements(vector<double>, int); //This returns the sum of all the values in the vector excluding those whose index is before the given index.  
		//double findMax(vector<double>); //This returns the maximum value in the vector.
		int numNZ(vector<int>); //This returns the number of non-zero values in the vector.
		double numNZ(vector<double>); //This returns the number of non-zero values in the vector.
};

/**************************************************************************************************/
//This Class stores the table of the confidence limits of the Student-T distribution.
class TDTable
{
	public:
		double getConfLimit(int,int);
};
/**************************************************************************************************/
#endif

