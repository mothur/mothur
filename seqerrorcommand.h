#ifndef SEQERRORCOMMAND
#define SEQERRORCOMMAND

/*
 *  seqerrorcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 7/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "sequence.hpp"
#include "referencedb.h"
#include "counttable.h"


class SeqErrorCommand : public Command {
public:
	SeqErrorCommand(string);
	SeqErrorCommand();
	~SeqErrorCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "seq.error";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nhttp://www.mothur.org/wiki/Seq.error"; }
	string getDescription()		{ return "seq.error"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;
	ReferenceDB* rdb;
	
	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
        ~linePair(){}
	};
	
    struct Compare {
        int AA, AT, AG, AC,	TA, TT, TG, TC,	GA, GT, GG, GC,	CA, CT, CG, CC,	NA, NT, NG, NC, Ai, Ti, Gi, Ci, Ni, dA, dT, dG, dC;
        string refName, queryName, sequence;
        double errorRate;
        int weight, matches, mismatches, total;
        
        Compare(){
            AA=0; AT=0; AG=0; AC=0;
            TA=0; TT=0; TG=0; TC=0;
            GA=0; GT=0; GG=0; GC=0;
            CA=0; CT=0; CG=0; CC=0;
            NA=0; NT=0; NG=0; NC=0;
            Ai=0; Ti=0; Gi=0; Ci=0; Ni=0;
            dA=0; dT=0; dG=0; dC=0;
            refName = "";
            queryName = "";
            weight = 1;
            matches = 0;
            mismatches = 0;
            total = 0;
            errorRate = 1.0000;
            sequence = "";
        }
        ~Compare(){};
    };

	vector<int> processIDS;   //processid
	vector<linePair> lines;
	vector<linePair> qLines;
	vector<linePair> rLines;

	void getReferences();
	map<string,int> getWeights();
	int getErrors(Sequence, Sequence, Compare&);
	void printErrorHeader(ofstream&);
	void printErrorData(Compare, int, ofstream&, ofstream&);
	void printSubMatrix();
	void printErrorFRFile(map<char, vector<int> >, map<char, vector<int> >);
	void printErrorQuality(map<char, vector<int> >);
	void printQualityFR(vector<vector<int> >, vector<vector<int> >);
	
	int setLines(string, string, string, vector<unsigned long long>&, vector<unsigned long long>&, vector<unsigned long long>&);
	int driver(string, string, string, string, string, string, linePair, linePair, linePair);
	int createProcesses(string, string, string, string, string, string);

	string queryFileName, referenceFileName, qualFileName, reportFileName, namesFileName, outputDir, countfile;
	double threshold;
	bool ignoreChimeras, save, aligned;
	int numRefs, processors;
	int maxLength, totalBases, totalMatches;
	//ofstream errorSummaryFile, errorSeqFile;
	vector<string> outputNames;
	
	vector<Sequence> referenceSeqs;
	vector<vector<int> > substitutionMatrix;
	vector<vector<int> > qualForwardMap;
	vector<vector<int> > qualReverseMap;
	vector<int> misMatchCounts;
	map<char, vector<int> > qScoreErrorMap;
	map<char, vector<int> > errorForward;
	map<char, vector<int> > errorReverse;
	map<string, int> weights;
	vector<string> megaAlignVector;

};

#endif
