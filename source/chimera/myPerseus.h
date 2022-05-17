#ifndef MOTHURPERSEUS
#define MOTHURPERSEUS

/*
 *  myPerseus.h
 *  
 *
 *  Created by Pat Schloss on 9/5/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */


#include "mothurout.h"

/**************************************************************************************************/
struct seqData {
	
	seqData(string name, string seq, int freq) : seqName(name), sequence(seq), frequency(freq) {	}
	
	bool operator<( seqData const& rhs ) const {
		
		bool verdict = 0;
		
		if(frequency < rhs.frequency){
			verdict = 1;
		}
		else if(frequency == rhs.frequency){
			verdict = (seqName > rhs.seqName);
		}
		
		return verdict; 
	}
	
	string seqName;
	string sequence;
	int frequency;
};
/**************************************************************************************************/
struct pwModel {
	pwModel(double m, double mm, double g): MATCH(m), MISMATCH(mm), GAP_OPEN(g) {;}
	double MATCH;
	double MISMATCH;
	double GAP_OPEN;	
};
/**************************************************************************************************/
struct pwAlign {
	pwAlign(): query(""), reference(""){}
	pwAlign(string q, string r): query(q), reference(r){}
	string query;
	string reference;
	
};
/**************************************************************************************************/
class Perseus {
	
public:
	Perseus() { m = MothurOut::getInstance(); }
	~Perseus() = default;
	
	vector<vector<double> > binomial(int);
	double modeledPairwiseAlignSeqs(string, string, string&, string&, vector<vector<double> >&);
	int getAlignments(int, vector<seqData>, vector<pwAlign>&, vector<vector<int> >& , vector<vector<int> >&, vector<vector<int> >&, vector<vector<int> >&, int&, int&, vector<bool>&);
	int getChimera(vector<seqData>,vector<vector<int> >&, vector<vector<int> >&,int&, int&, int&,vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<bool>);
	string stitchBimera(vector<pwAlign>&, int, int, int, vector<vector<int> >&, vector<vector<int> >&);
	int getTrimera(vector<seqData>&, vector<vector<int> >&, int&, int&, int&, int&, int&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<bool>);
	string stitchTrimera(vector<pwAlign>, int, int, int, int, int, vector<vector<int> >&, vector<vector<int> >&);
	double calcLoonIndex(string, string, string, int, vector<vector<double> >&);
	double classifyChimera(double, double, double, double, double);
	
private:
	MothurOut* m;
	int toInt(char);
	double basicPairwiseAlignSeqs(string, string, string&, string&, pwModel);
	int getDiffs(string, string, vector<int>&, vector<int>&, vector<int>&, vector<int>&);
	int getLastMatch(char, vector<vector<char> >&, int, int, string&, string&);
	int threeWayAlign(string, string, string, string&, string&, string&);
	double calcBestDistance(string, string);

	
};
/**************************************************************************************************/
#endif


