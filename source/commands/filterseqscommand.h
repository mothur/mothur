#ifndef FILTERSEQSCOMMAND_H
#define FILTERSEQSCOMMAND_H

/*
 *  filterseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "filters.h"

class Sequence;
class FilterSeqsCommand : public Command {

public:
	FilterSeqsCommand(string);
	FilterSeqsCommand();
	~FilterSeqsCommand() {};
	
	vector<string> setParameters();
	string getCommandName()			{ return "filter.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Filter.seqs"; }
	string getDescription()		{ return "removes columns from alignments based on a criteria defined by the user"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:

    vector< vector<unsigned long long> >  savedPositions;

	string vertical, filter, fasta, hard, outputDir, filterFileName;
	vector<string> fastafileNames;	
	int alignmentLength, processors;
	vector<int> bufferSizes;
	vector<string> outputNames;

	char trump;
	bool abort, recalced;
	float soft;
	long long numSeqs;
	
	string createFilter();
	int filterSequences();
	long long createProcessesCreateFilter(Filters&, string);
	long long createProcessesRunFilter(string, string, string, vector<linePair>);	
};


/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct filterData {
	Filters F;
    int alignmentLength, threadid;
    unsigned long long start, end;
    long long count;
    MothurOut* m;
    string filename, hard;
    char trump;
    float soft;
    bool vertical;
	
	filterData(){}
	filterData(string fn, MothurOut* mout, unsigned long long st, unsigned long long en, int aLength, char tr, bool vert, float so, string ha, int tid) {
        filename = fn;
		m = mout;
		start = st;
		end = en;
        trump = tr;
        alignmentLength = aLength;
        vertical = vert;
        soft = so;
        hard = ha;
		count = 0;
        threadid = tid;
	}
};
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct filterRunData {
    int alignmentLength;
    unsigned long long start, end;
    long long count;
    MothurOut* m;
    string filename;
    string filter, outputFilename;
	
	filterRunData(){}
	filterRunData(string f, string fn, string ofn, MothurOut* mout, unsigned long long st, unsigned long long en, int aLength) {
        filter = f;
        outputFilename = ofn;
        filename = fn;
		m = mout;
		start = st;
		end = en;
        alignmentLength = aLength;
		count = 0;
	}
};
/**************************************************************************************************/

#endif
