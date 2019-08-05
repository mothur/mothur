#ifndef Mothur_makecontigscommand_h
#define Mothur_makecontigscommand_h

//
//  makecontigscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "command.hpp"
#include "sequence.hpp"
#include "qualityscores.h"
#include "alignment.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"
#include "trimoligos.h"
#include "oligos.h"
#include "fastqread.h"
#include "kmeralign.h"
#include "splitgroupscommand.h"
#include "filefile.hpp"


#        define PROBABILITY(score) (pow(10.0, (-(double)(score)) / 10.0))
#        define PHREDMAX 46
#        define PHREDCLAMP(x) ((x) > PHREDMAX ? PHREDMAX : ((x) < 0 ? 0 : (x)))

struct pairFastqRead {
	FastqRead forward;
    FastqRead reverse;
    FastqRead findex;
    FastqRead rindex;

	pairFastqRead() {};
	pairFastqRead(FastqRead f, FastqRead r) : forward(f), reverse(r){};
    pairFastqRead(FastqRead f, FastqRead r, FastqRead fi, FastqRead ri) : forward(f), reverse(r), findex(fi), rindex(ri) {};
	~pairFastqRead() {};
};
/**************************************************************************************************/

class MakeContigsCommand : public Command {
public:
    MakeContigsCommand(string);
    MakeContigsCommand();
    ~MakeContigsCommand(){}

    vector<string> setParameters();
    string getCommandName()			{ return "make.contigs";			}
    string getCommandCategory()		{ return "Sequence Processing";		}
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden

	string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Make.contigs"; }
    string getDescription()		{ return "description"; }

    int execute();
    void help() { m->mothurOut(getHelpString()); }

private:

#define perfectMatch 2
#define poundMatch  1
#define offByOne  3

    char delim;
    bool abort, allFiles, trimOverlap, createFileGroup, createOligosGroup, makeCount, noneOk, reorient, gz;
    string outputDir, ffastqfile, rfastqfile, align, oligosfile, rfastafile, ffastafile, rqualfile, fqualfile, findexfile, rindexfile, file, format, inputDir;
	float match, misMatch, gapOpen, gapExtend, maxee;
	int processors, longestBase, insert, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs, deltaq, kmerSize, nameType, offByOneTrimLength;
    vector<string> outputNames;
    set<string> badNames;
	map<string, int> groupCounts;
    map<string, string> groupMap;

    unsigned long long processMultipleFileOption(string& outFastaFile, string&);
    unsigned long long processSingleFileOption(string& outFastaFile, string& outScrapFastaFile, string& outQualFile, string& outScrapQualFile, string& outMisMatchFile, string group);

    unsigned long long createProcesses(vector<string>, vector<string>, string, string, string, string, string, vector<vector<string> >, vector<vector<string> >, string, map<int, oligosPair>& pairedPrimers, map<int, oligosPair>& rpairedPrimers, map<int, oligosPair>& pairedBarcodes, map<int, oligosPair>& rpairedBarcodes, vector<string>& barcodeNames, vector<string>& primerNames);
    unsigned long long createProcessesGroups(vector< vector<string> >, string compositeFastaFile, string compositeScrapFastaFile, string compositeQualFile, string compositeScrapQualFile, string compositeMisMatchFile, map<int, string>& file2Groups);

    int createGroupFile(string outputGroupFile, string resultFastafile);
    vector< vector<string> > readFileNames(string, map<int, string>&);
    bool getOligos(map<int, oligosPair>& pairedPrimers, map<int, oligosPair>& rpairedPrimers, map<int, oligosPair>& pairedBarcodes, map<int, oligosPair>& rpairedBarcodes, vector<string>& barcodeNames, vector<string>& primerNames);
    int setLines(vector<string>, vector<string>, vector<linePair>& fastaFilePos, vector<linePair>& qfileFilePos, char delim); //the delim let you know whether this is fasta and qual, or fastq and index. linePair entries will always be in sets of two. One for the forward and one for hte reverse.  (fastaFilePos[0] - ffasta, fastaFilePos[1] - rfasta) - processor1
    bool testGZReadable(vector<string>&, vector<string>&, bool&);

};


/**************************************************************************************************/

#endif
