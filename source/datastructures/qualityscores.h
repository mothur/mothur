#ifndef QUALITYSCORES
#define QUALITYSCORES

/*
 *  qualityscores.h
 *  Mothur
 *
 *  Created by Pat Schloss on 7/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

//DataStructure for a quality file.


#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"
#include "utils.hpp"
#include "writer.h"

/**************************************************************************************************/

class QualityScores {
public:
	QualityScores();
    ~QualityScores() = default;
    QualityScores(string n, vector<int> qs);
	QualityScores(ifstream&);
    #ifdef USE_BOOST
    QualityScores(boost::iostreams::filtering_istream&);
    #endif
    int read(ifstream&);
	string getName();
	int getLength(){    return (int)qScores.size();  }
	//vector<int> getQualityScores() { return qScores; }
	void printQScores(ofstream&);
    void printQScores(ostream&);
    void printQScores(OutputWriter*);
	void trimQScores(int, int);
	void flipQScores();
	bool stripQualThreshold(Sequence&, double);
	bool stripQualRollingAverage(Sequence&, double, bool);
	bool stripQualWindowAverage(Sequence&, int, int, double, bool);
	bool cullQualAverage(Sequence&, double, bool);
	void updateQScoreErrorMap(map<char, vector<int> >&, string, int, int, int);
	void updateForwardMap(vector<vector<int> >&, int, int, int);
	void updateReverseMap(vector<vector<int> >&, int, int, int);
    void setName(string n);
    void setScores(vector<int> qs) { qScores = qs; seqLength = (int)qScores.size(); }
    vector<int> getScores() { return qScores; }

private:

	double calculateAverage(bool);
	double calculateExpectedErrors(void);
	MothurOut* m;
	vector<int> qScores;
    Utils util;

	string seqName;
	int seqLength;

    string getSequenceName(ifstream&);
    string getCommentString(ifstream&);

    #ifdef USE_BOOST
    string getCommentString(boost::iostreams::filtering_istream&);
    string getSequenceName(boost::iostreams::filtering_istream&);
    #endif
};

/**************************************************************************************************/

#endif
