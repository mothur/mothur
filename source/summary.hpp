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
    
    Summary() { m = MothurOut::getInstance(); total = 0; numUniques = 0; hasNameOrCount = false; nameCountNumUniques = 0; type = "count"; }
    Summary(string n); //provide name or count file to include in counts
    ~Summary() {}
    
    long long summarize(string f, string o); //provide fasta file to summarize (paralellized) and output file for individual seqs info. To skip output file, o = ""
    string addSeq(Sequence); //return summary output line
    vector<long long> getDefaults();

    vector<long long> getStart(); //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getStart(double value); //2.5 = 2.5% of sequences of sequences start before, 25 = location 25% of sequences start before
    vector<long long> getEnd(); //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getEnd(double value); //2.5 = 2.5% of sequences of sequences end after, 25 = location 25% of sequences end after
    vector<long long> getAmbig(); //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getAmbig(double value); //25 = max abigous bases 25% of sequences contain
    vector<long long> getLength(); //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getLength(double value); // 25 = min length of 25% of sequences
    vector<long long> getHomop(); //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    long long getHomop(double value);
    
    long long getTotalSeqs() { return total; }
    long long getUniqueSeqs() { return numUniques; }
    
private:
    
    MothurOut* m;
    long long total, numUniques, nameCountNumUniques;
    bool hasNameOrCount;
    string type;
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    map<string, int> nameMap;
    map<string, int>::iterator itFindName;
    map<int, long long>::iterator it;
    
    int driverSummarize(string, string, linePair lines);
    
    
};
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct seqSumData {
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    string filename, summaryFile;
    unsigned long long start;
    unsigned long long end;
    long long count;
    long long total;
    MothurOut* m;
    bool hasNameMap;
    map<string, int> nameMap;
    
    
    seqSumData(){}
    seqSumData(string f, string sum, MothurOut* mout, unsigned long long st, unsigned long long en, bool na, map<string, int> nam) {
        filename = f;
        m = mout;
        start = st;
        end = en;
        hasNameMap = na;
        nameMap = nam;
        count = 0;
        summaryFile = sum;
    }
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MySeqSumThreadFunction(LPVOID lpParam){
    seqSumData* pDataArray;
    pDataArray = (seqSumData*)lpParam;
    
    try {
        ofstream out;
        if (pDataArray->summaryFile != "") { pDataArray->m->openOutputFile(pDataArray->summaryFile, out); }
        
        ifstream in;
        pDataArray->m->openInputFile(pDataArray->filename, in);
        
        //print header if you are process 0
        if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
            in.seekg(0);
            pDataArray->m->zapGremlins(in);
        }else { //this accounts for the difference in line endings.
            in.seekg(pDataArray->start-1); pDataArray->m->gobble(in);
        }
        
        for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
            
            if (pDataArray->m->control_pressed) { in.close(); pDataArray->count = 1; return 1; }
            
            Sequence current(in); pDataArray->m->gobble(in);
            
            if (current.getName() != "") {
                
                long long num = 1;
                if (pDataArray->hasNameMap){
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator it = pDataArray->nameMap.find(current.getName());
                    
                    if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + current.getName() + " is not in your name or count file, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; }
                    else { num = it->second; }
                }
                
                int thisStartPosition = current.getStartPos();
                map<int, long long>::iterator it = pDataArray->startPosition.find(thisStartPosition);
                if (it == pDataArray->startPosition.end()) { pDataArray->startPosition[thisStartPosition] = num; } //first finding of this start position, set count.
                else { it->second += num; } //add counts
                
                int thisEndPosition = current.getEndPos();
                it = pDataArray->endPosition.find(thisEndPosition);
                if (it == pDataArray->endPosition.end()) { pDataArray->endPosition[thisEndPosition] = num; } //first finding of this end position, set count.
                else { it->second += num; } //add counts
                
                int thisSeqLength = current.getNumBases();
                it = pDataArray->seqLength.find(thisSeqLength);
                if (it == pDataArray->seqLength.end()) { pDataArray->seqLength[thisSeqLength] = num; } //first finding of this length, set count.
                else { it->second += num; } //add counts
                
                int thisAmbig = current.getAmbigBases();
                it = pDataArray->ambigBases.find(thisAmbig);
                if (it == pDataArray->ambigBases.end()) { pDataArray->ambigBases[thisAmbig] = num; } //first finding of this ambig, set count.
                else { it->second += num; } //add counts
                
                int thisHomoP = current.getLongHomoPolymer();
                it = pDataArray->longHomoPolymer.find(thisHomoP);
                if (it == pDataArray->longHomoPolymer.end()) { pDataArray->longHomoPolymer[thisHomoP] = num; } //first finding of this homop, set count.
                else { it->second += num; } //add counts
                
                pDataArray->count++;
                
                string output = "";
                output += seq.getName() + '\t';
                output += thisStartPosition + '\t' + thisEndPosition + '\t';
                output += thisSeqLength + '\t' + thisAmbig + '\t';
                output += thisHomoP + '\t' + num;
                
                if (pDataArray->summaryFile != "") { out << output << endl; }
            }
        }
        
        if (pDataArray->summaryFile != "") { out.close(); }
        in.close();
        
        return 0;
        
    }
    catch(exception& e) {
        pDataArray->m->errorOut(e, "SeqSummaryCommand", "MySeqSumThreadFunction");
        exit(1);
    }
} 
#endif




#endif /* summary_hpp */
