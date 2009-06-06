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
		initCommandParameters();		
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
bool ValidParameters::isValidParameter(string parameter, string command, string value) {
	try {	
		bool valid = false;
		vector<string> cParams = commandParameters[command];
		int numParams = cParams.size(); 
		for(int i = 0; i < numParams; i++) {
			if(cParams.at(i).compare(parameter) == 0) {
				valid = true;
				i = numParams;
			}
		}
		if(!valid) {
			cout << "'" << parameter << "' is not a valid parameter for the " << command << " command.\n";
			cout << "The valid paramters for the " << command << " command are: ";
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

/***********************************************************************/

/***********************************************************************/
void ValidParameters::initCommandParameters() {
	try {	
		//{"parameter1","parameter2",...,"last parameter"};
		
		string readdistArray[] = {"phylip","column", "name","cutoff","precision", "group"};
		commandParameters["read.dist"] = addParameters(readdistArray, sizeof(readdistArray)/sizeof(string));

		string readotuArray[] =  {"list","order","shared", "line", "label","group","sabund", "rabund"};
		commandParameters["read.otu"] = addParameters(readotuArray, sizeof(readotuArray)/sizeof(string));
		
		string readtreeArray[] = {"tree","group"};
		commandParameters["read.tree"] = addParameters(readtreeArray, sizeof(readtreeArray)/sizeof(string));
		
		string clusterArray[] =  {"cutoff","precision","method"};
		commandParameters["cluster"] = addParameters(clusterArray, sizeof(clusterArray)/sizeof(string));
		
		string deconvoluteArray[] =  {"fasta"};
		commandParameters["unique.seqs"] = addParameters(deconvoluteArray, sizeof(deconvoluteArray)/sizeof(string));
		
		string collectsingleArray[] =  {"freq","line","label","calc","abund","size"};
		commandParameters["collect.single"] = addParameters(collectsingleArray, sizeof(collectsingleArray)/sizeof(string));

		string collectsharedArray[] =  {"freq","line","label","calc","groups"};
		commandParameters["collect.shared"] = addParameters(collectsharedArray, sizeof(collectsharedArray)/sizeof(string));

		string getgroupArray[] =  {};
		commandParameters["get.group"]	 = addParameters(getgroupArray, sizeof(getgroupArray)/sizeof(string));
		
		string getlabelArray[] =  {};
		commandParameters["get.label"]	= addParameters(getlabelArray, sizeof(getlabelArray)/sizeof(string));

		string getlineArray[] =  {};
		commandParameters["get.line"] = addParameters(getlineArray, sizeof(getlineArray)/sizeof(string));
		
		string getsabundArray[] =  {"line", "label"};
		commandParameters["get.sabund"] = addParameters(getsabundArray, sizeof(getsabundArray)/sizeof(string));
		
		string getrabundArray[] =  {"line", "label"};
		commandParameters["get.rabund"] = addParameters(getrabundArray, sizeof(getrabundArray)/sizeof(string));

		string rarefactionsingleArray[] =  {"iters","freq","line","label","calc","abund"};
		commandParameters["rarefaction.single"] = addParameters(rarefactionsingleArray, sizeof(rarefactionsingleArray)/sizeof(string));

		string rarefactionsharedArray[] =  {"iters","jumble","line","label","calc","groups"};
		commandParameters["rarefaction.shared"] = addParameters(rarefactionsharedArray, sizeof(rarefactionsharedArray)/sizeof(string));
		
		string libshuffArray[] =  {"iters","groups","step","form","cutoff"};
		commandParameters["libshuff"] = addParameters(libshuffArray, sizeof(libshuffArray)/sizeof(string));
		
		string summarysingleArray[] =  {"line","label","calc","abund","size"};
		commandParameters["summary.single"] = addParameters(summarysingleArray, sizeof(summarysingleArray)/sizeof(string));

		string summarysharedArray[] =  {"line","label","calc","groups"};
		commandParameters["summary.shared"] = addParameters(summarysharedArray, sizeof(summarysharedArray)/sizeof(string));

		string parsimonyArray[] =  {"random","groups","iters"};
		commandParameters["parsimony"] = addParameters(parsimonyArray, sizeof(parsimonyArray)/sizeof(string));

		string unifracWeightedArray[] =  {"groups","iters"};
		commandParameters["unifrac.weighted"] = addParameters(unifracWeightedArray, sizeof(unifracWeightedArray)/sizeof(string));

		string unifracUnweightedArray[] =  {"groups","iters"};
		commandParameters["unifrac.unweighted"] = addParameters(unifracUnweightedArray, sizeof(unifracUnweightedArray)/sizeof(string));

		string heatmapArray[] =  {"groups","line","label","sorted","scale"};
		commandParameters["heatmap"] = addParameters(heatmapArray, sizeof(heatmapArray)/sizeof(string));
		
		string filterseqsArray[] =  {"fasta", "trump", "soft", "hard", "vertical"};
		commandParameters["filter.seqs"] = addParameters(filterseqsArray, sizeof(filterseqsArray)/sizeof(string));

		string summaryseqsArray[] =  {"fasta"};
		commandParameters["summary.seqs"] = addParameters(summaryseqsArray, sizeof(summaryseqsArray)/sizeof(string));

		string screenseqsArray[] =  {"fasta", "start", "end", "maxambig", "maxhomop", "minlength", "maxlength", "name", "group"};
		commandParameters["screen.seqs"] = addParameters(screenseqsArray, sizeof(screenseqsArray)/sizeof(string));

		string reverseseqsArray[] =  {"fasta"};
		commandParameters["reverse.seqs"] = addParameters(reverseseqsArray, sizeof(reverseseqsArray)/sizeof(string));

		string vennArray[] =  {"groups","line","label","calc"};
		commandParameters["venn"] = addParameters(vennArray, sizeof(vennArray)/sizeof(string));
		
		string binseqsArray[] =  {"fasta","line","label","name", "group"};
		commandParameters["bin.seqs"] = addParameters(binseqsArray, sizeof(binseqsArray)/sizeof(string));
		
		string distsharedArray[] =  {"line","label","calc","groups"};
		commandParameters["dist.shared"] = addParameters(distsharedArray, sizeof(distsharedArray)/sizeof(string));
		
		string getOTURepArray[] =  {"fasta","list","line","label","name", "group"};
		commandParameters["get.oturep"] = addParameters(getOTURepArray, sizeof(getOTURepArray)/sizeof(string));
		
		string treeGroupsArray[] =  {"line","label","calc","groups", "phylip", "column", "name"};
		commandParameters["tree.shared"] = addParameters(treeGroupsArray, sizeof(treeGroupsArray)/sizeof(string));
		
		string bootstrapArray[] =  {"line","label","calc","groups","iters"};
		commandParameters["bootstrap.shared"] = addParameters(bootstrapArray, sizeof(bootstrapArray)/sizeof(string));
		
		string concensusArray[] =  {};
		commandParameters["concensus"] = addParameters(concensusArray, sizeof(concensusArray)/sizeof(string));
		
		string distanceArray[] =  {"fasta", "calc", "countends", "cutoff", "processors"};
		commandParameters["dist.seqs"] = addParameters(distanceArray, sizeof(distanceArray)/sizeof(string));
		
		string AlignArray[] =  {"fasta", "candidate", "search", "ksize", "align", "match", "mismatch", "gapopen", "gapextend"};
		commandParameters["align.seqs"] = addParameters(AlignArray, sizeof(AlignArray)/sizeof(string));
		
		string quitArray[] = {};
		commandParameters["quit"] = addParameters(quitArray, sizeof(quitArray)/sizeof(string));

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

		string jumbleArray[] = {">","0", "<","1", "only"};
		parameterRanges["jumble"] = addParameters(jumbleArray, rangeSize);

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

