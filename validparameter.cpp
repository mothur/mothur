/*
 *  validparameter.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validparameter.h"

/***********************************************************************/

ValidParameters::ValidParameters() {
	try {
		initParameterRanges();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function ValidParameters. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function ValidParameters. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

ValidParameters::~ValidParameters() {}

/***********************************************************************/
bool ValidParameters::isValidParameter(string parameter, vector<string> cParams, string value) {
	try {	
		bool valid = false;
		//vector<string> cParams = commandParameters[command];
		int numParams = cParams.size(); 
		for(int i = 0; i < numParams; i++) {
			if(cParams.at(i).compare(parameter) == 0) {
				valid = true;
				i = numParams;
			}
		}
		if(!valid) {
			cout << "'" << parameter << "' is not a valid parameter." << endl;
			cout << "The valid parameters are: ";
			for(int i = 0; i < numParams-1; i++)
				cout << cParams.at(i) << ", ";
			cout << "and " << cParams.at(numParams-1) << ".\n";
			return false;
		}
		
		if(parameterRanges.count(parameter) != 1)
			return true;
	
		int pVal;
		double piSentinel = 3.14159;
		vector<string> range = parameterRanges[parameter];
		
		vector<string> values;
		splitAtDash(value, values);
		
		for(int i = 0; i < values.size(); i++) {
			value = values.at(i);
			valid = convertTest(value, pVal);
		
			if(!valid)
				return false;
			
			
			
			/********************************************************************************************************
				   Special Cases
			*********************************************************************************************************/
			
			if(parameter.compare("precision") == 0) {
				double logNum = log10((double)pVal);
				double diff = (double)((int)logNum - logNum);
				if(diff != 0) {
					cout << "The precision parameter can only take powers of 10 as a value (e.g. 10,1000,1000, etc.)\n";
					return false;
				}
			}
			
			/************************************************************************************************************/
			
			
			
			double a,b,c,d,e;
			
			if(range.at(1).compare("NA") == 0)
				a = piSentinel;
			else
				a = atoi(range.at(1).c_str()); 
				
			if(range.at(3).compare("NA") == 0)
				b = piSentinel;
			else
				b = atoi(range.at(3).c_str()); 
						
			if(range.at(4).compare("between") == 0)
				c = 0;
			else if(range.at(4).compare("only") == 0)
				c = 1;
			else {
				cout << "The range can only be 'between' or 'only' the bounding numbers.\n";
				return false;
			}
			
			if(range.at(0).compare(">") == 0)
				d = 0;
			else if(range.at(0).compare(">=") == 0 || range[3].compare("=>") == 0)
				d = 1;
			else {
				cout << "The parameter value can only be '>', '>=', or '=>' the lower bounding number.\n";
				return false;
			}
			
			if(range.at(2).compare("<") == 0)
				e = 0;
			else if(range.at(2).compare("<=") == 0 || range[4].compare("=<") == 0)
				e = 1;
			else {
				cout << "The parameter value can only be '<', '<=', or '=<' the upper bounding number.\n";
				return false;
			}
			
			bool a0 = pVal > a;
			bool a1 = pVal >= a;
			bool b0 = pVal < b;
			bool b1 = pVal <= b;
			
			if(c != 1) {
				if(a != piSentinel && b == piSentinel) {
					if(d == 0)
						valid = a0;
					else
						valid = a1;
				}
				else if(a == piSentinel && b != piSentinel) {
					if(e == 0)
						valid = b0;
					else
						valid = b1;
				}
				else {
					if(d == 0 && e == 0)
						valid = (a0 && b0);
					else if(d == 0 && e == 1)
						valid = (a0 && b1);
					else if(d == 1 && e == 0)
						valid = (a1 && b0);
					else
						valid = (a1 && b1);
				}
			}
			else {
				if(a != piSentinel && b == piSentinel)
					valid = (pVal == a);
				else if(a == piSentinel && b != piSentinel)
					valid = (pVal == b);
				else
					valid = (pVal == a || pVal == b);
			}
			
			
			if(!valid) {
				cout << "The '" << parameter << "' parameter needs to be ";
				if(c == 1)
					cout << "either '" << a << "' or '" << b << "'.\n";
				else {
					if(a != piSentinel) {
						cout << ">";
						if(d != 0)
							cout << "=";
						cout << " '" << a << "'";
					}
					if(b == piSentinel)
						cout << "'.\n";
					else if(a != piSentinel)
						cout << " and ";
					if(b != piSentinel) {
						cout << "<";
						if(e != 0)
							cout << "=";
						cout << " '" << b << "'.\n";
					}
				}
				return false;
			}
		}
		return true;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************/

/******************************************************/

string ValidParameters::validFile(map<string, string> container, string parameter, bool isFile) {
	try {
		int ableToOpen;
		ifstream in;
		map<string, string>::iterator it;
		
		it = container.find(parameter);
		if(it != container.end()){ //no parameter given
			if(isFile == true) {
				ableToOpen = openInputFile(it->second, in);
				if (ableToOpen == 1) { return "not open"; }
				in.close();
			}
		}else { return "not found"; }
		
		return it->second;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function validFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function validFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

/***********************************************************************/
void ValidParameters::initParameterRanges() {
	try {	
		int rangeSize = 5;

		/**************************************************************************************************************
			{">=" or "=>" or ">" if the value should be greater than or equal to or just greater than the lower bound,
		    A number representing the lower bound ("NA" if there is no lower bound), 
		   "<=" or "=<" or "<" if the value shoud be less than or equal to or just less than the upper bound,
		    A number representing the upper bound ("NA" if there is no lower bound),
		   "between" if between lower and upper bounds or "only" if exactly one of the bounds};
		   
		   # = parameter
		   # (>, >=) lower bound, # (<, <=) upperbound, # should be (between, only) lower and upper bounds.
		   ***********************************************************************************************************/
		
		string precisionArray[] = {">=","10", "<","NA", "between"};
		parameterRanges["precision"] = addParameters(precisionArray, rangeSize);
		
		string itersArray[] = {">=","10", "<","NA", "between"};
		parameterRanges["iters"] = addParameters(itersArray, rangeSize);

		string freqArray[] = {">=","1", "<","NA", "between"};
		parameterRanges["freq"] = addParameters(freqArray, rangeSize);

		//string lineArray[] = {">=","1", "<","NA", "between"};
		//parameterRanges["line"] = addParameters(lineArray, rangeSize);

		string abundArray[] = {">=","5", "<","NA", "between"};
		parameterRanges["abund"] = addParameters(abundArray, rangeSize);
		
		string softArray[] = {">=","0", "<=","100", "between"};
		parameterRanges["soft"] = addParameters(softArray, rangeSize);
		
		string sizeArray[] = {">=","1", "<","NA", "between"};
		parameterRanges["size"] = addParameters(sizeArray, rangeSize);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

/***********************************************************************/
vector<string> ValidParameters::addParameters(string parameters[], int size) {
	try {	
		vector<string> pVector (parameters, parameters+size); 
		return pVector;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

