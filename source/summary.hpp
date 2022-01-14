//
//  summary.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include "mothurout.h"
#include "sequence.hpp"
#include "counttable.h"

class Summary {

public:

#ifdef UNIT_TEST
    friend class TestSummary;
#endif

    Summary(int p) { processors = p; m = MothurOut::getInstance(); total = 0; numUniques = 0; hasNameOrCount = false; nameCountNumUniques = 0; type = "count"; }
    ~Summary() {}

    long long summarizeFasta(string f, string n, string o); //provide fasta file to summarize (paralellized) and optional nameorCountfile and optional outputfile for individual seqs info. To skip nameCount or output file, n="" and / or o=""
    long long summarizeFasta(string f, string o); //provide fasta file to summarize (paralellized) and optional outputfile for individual seqs info. To skip output file, o=""
    long long summarizeFastaSummary(string f); //provide summary of fasta file to summarize (paralellized)
    long long summarizeFastaSummary(string f, string n); //provide summary of fasta file and name or count file to summarize (paralellized)
    long long summarizeContigsSummary(string f); //provide summary of contigs summary file to summarize (paralellized)
    long long summarizeContigsSummary(string f, string n); //provide summary of contigs summary file and name or count file to summarize (paralellized)
    long long summarizeAlignSummary(string f); //provide summary of contigs summary file to summarize (paralellized)
    long long summarizeAlignSummary(string f, string n); //provide summary of contigs summary file and name or count file to summarize (paralellized)

    vector<long long> getDefaults();
    //fasta and summary
    vector<long long> getStart() { return (getValues(startPosition)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getStart(double value) { return (getValue(startPosition, value)); } //2.5 = 2.5% of sequences of sequences start before, 25 = location 25% of sequences start before
    vector<long long> getEnd() { return (getValues(endPosition)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getEnd(double value) { return (getValue(endPosition, value)); } //2.5 = 2.5% of sequences of sequences end after, 25 = location 25% of sequences end after
    vector<long long> getAmbig() { return (getValues(ambigBases)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getAmbig(double value) { return (getValue(ambigBases, value)); } //25 = max abigous bases 25% of sequences contain
    vector<long long> getLength() { return (getValues(seqLength)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getLength(double value) { return (getValue(seqLength, value)); } // 25 = min length of 25% of sequences
    vector<long long> getHomop() { return (getValues(longHomoPolymer)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getHomop(double value) { return (getValue(longHomoPolymer, value)); }

    //contigs
    vector<long long> getOStart() { return (getValues(ostartPosition)); } //contigs overlap start - returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getOStart(double value) { return (getValue(ostartPosition, value)); } //contigs overlap start - 2.5 = 2.5% of sequences of sequences start before, 25 = location 25% of sequences start before
    vector<long long> getOEnd() { return (getValues(oendPosition)); } //contigs overlap end -returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getOEnd(double value) { return (getValue(oendPosition, value)); } //contigs overlap end -2.5 = 2.5% of sequences of sequences end after, 25 = location 25% of sequences end after
    vector<long long> getOLength() { return (getValues(oseqLength)); } //contigs overlap length - returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getOLength(double value) { return (getValue(oseqLength, value)); } //contigs overlap length - 25 = min length of 25% of sequences
    vector<long long> getMisMatches() { return (getValues(misMatches)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getMisMatches(double value) { return (getValue(misMatches, value)); }
    vector<long long> getNumNs() { return (getValues(numNs)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getNumNs(double value) { return (getValue(numNs, value)); } //25 = max abigous bases 25% of sequences contain
    vector<long long> getSims() { return (getValues(sims)); } //contigs overlap length - returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getSims(double value) { return (getValue(sims, value)); } //contigs overlap length - 25 = min length of 25% of sequences
    vector<long long> getScores() { return (getValues(scores)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getScores(double value) { return (getValue(scores, value)); }
    vector<long long> getNumInserts() { return (getValues(inserts)); } //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getNumInserts(double value) { return (getValue(inserts, value)); } //25 = max abigous bases 25% of sequences contain
		int getMaxAbundance();

    long long getTotalSeqs() { return total; }
    long long getUniqueSeqs() { return numUniques; }

private:

    MothurOut* m;
    Utils util;
    int processors;
    long long total, numUniques, nameCountNumUniques;
    bool hasNameOrCount;
    string type;
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    map<int, long long> ostartPosition;
    map<int, long long> oendPosition;
    map<int, long long> oseqLength;
    map<int, long long> misMatches;
    map<int, long long> numNs;
    map<float, long long> sims;
    map<float, long long> scores;
    map<int, long long> inserts;

    map<string, int> nameMap;
    map<int, long long>::iterator it;

    void processNameCount(string n); //determines whether name or count and fills nameMap, ignored if n = ""
    vector<long long> getValues(map<int, long long>& positions);
    long long getValue(map<int, long long>& positions, double);
    vector<long long> getValues(map<float, long long>& positions);
    long long getValue(map<float, long long>& positions, double);
    bool isCountFile(string);



};
/**************************************************************************************************/
struct seqSumData {
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    map<int, long long> ostartPosition;
    map<int, long long> oendPosition;
    map<int, long long> oseqLength;
    map<int, long long> misMatches;
    map<int, long long> numNs;
    map<float, long long> sims;
    map<float, long long> scores;
    map<int, long long> inserts;


    string filename, summaryFile, contigsfile, output;
    double start;
    double end;
    long long count;
    long long total;
    MothurOut* m;
    bool hasNameMap;
    map<string, int> nameMap;
    Utils util;


    seqSumData(){}
    //FastaSummarize - output file created
    seqSumData(string f, string sum, double st, double en, bool na, map<string, int> nam) {
        filename = f;
        m = MothurOut::getInstance();
        start = st;
        end = en;
        hasNameMap = na;
        nameMap = nam;
        count = 0;
        total = 0;
        summaryFile = sum;
    }

    //FastaSummarySummarize - no output files
    seqSumData(string f, double st, double en, bool na, map<string, int> nam) {
        filename = f;
        m = MothurOut::getInstance();
        start = st;
        end = en;
        hasNameMap = na;
        nameMap = nam;
        count = 0;
        total = 0;
    }
};

/**************************************************************************************************/

#endif /* summary_hpp */
