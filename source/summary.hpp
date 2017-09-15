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
    map<int, long long> ostartPosition;
    map<int, long long> oendPosition;
    map<int, long long> oseqLength;
    map<int, long long> misMatches;
    map<int, long long> numNs;
    map<float, long long> sims;
    map<float, long long> scores;
    map<int, long long> inserts;

    map<string, int> nameMap;
    map<string, int>::iterator itFindName;
    map<int, long long>::iterator it;
    
    string addSeq(Sequence); //return summary output line
    string addSeq(string name, int length, int olength, int ostart, int oend, int mismatches, int numns);
    string addSeq(string name, int start, int end, int length, int ambigs, int polymer, long long numReps);
    string addSeq(string name, int length, float SimBtwnQueryTemplate, float SearchScore, int LongestInsert);
    //int driverSummarize(string, string, linePair lines); //fastafile, outputfile (optional set to "" to ignore), file positions
    void processNameCount(string n); //determines whether name or count and fills nameMap, ignored if n = ""
    int driverFastaSummarySummarize(string, linePair lines); //summaryfile, file positions
    int driverContigsSummarySummarize(string, linePair lines); //summaryfile, file positions
    int driverAlignSummarySummarize(string, linePair lines); //summaryfile, file positions
    vector<long long> getValues(map<int, long long>& positions);
    long long getValue(map<int, long long>& positions, double);
    vector<long long> getValues(map<float, long long>& positions);
    long long getValue(map<float, long long>& positions, double);


    
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
    unsigned long long start;
    unsigned long long end;
    long long count;
    long long total;
    MothurOut* m;
    bool hasNameMap;
    map<string, int> nameMap;
    
    
    seqSumData(){}
    //FastaSummarize
    seqSumData(string f, string sum, MothurOut* mout, unsigned long long st, unsigned long long en, bool na, map<string, int> nam) {
        filename = f;
        m = mout;
        start = st;
        end = en;
        hasNameMap = na;
        nameMap = nam;
        count = 0;
        total = 0;
        summaryFile = sum;
    }
    
    //FastaSummarySummarize
    seqSumData(string f, MothurOut* mout, unsigned long long st, unsigned long long en, bool na, map<string, int> nam) {
        filename = f;
        m = mout;
        start = st;
        end = en;
        hasNameMap = na;
        nameMap = nam;
        count = 0;
        total = 0;
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
            out << "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs" << endl;
        }else { //this accounts for the difference in line endings.
            in.seekg(pDataArray->start-1); pDataArray->m->gobble(in);
        }
        
        for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
            
            if (pDataArray->m->getControl_pressed()) { in.close(); pDataArray->count = 1; return 1; }
            
            Sequence current(in); pDataArray->m->gobble(in);
            
            if (current.getName() != "") {
                
                long long num = 1;
                if (pDataArray->hasNameMap){
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator it = pDataArray->nameMap.find(current.getName());
                    
                    if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + current.getName() + " is not in your name or count file, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->setControl_pressed(true); }
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
                
                int thisNumNs = current.getNumNs();
                it = pDataArray->numNs.find(thisNumNs);
                if (it == pDataArray->numNs.end()) { pDataArray->numNs[thisNumNs] = num; } //first finding of this homop, set count.
                else { it->second += num; } //add counts
                
                pDataArray->count++;
                
                string output = "";
                output += current.getName() + '\t';
                output += toString(thisStartPosition) + '\t' + toString(thisEndPosition) + '\t';
                output += toString(thisSeqLength) + '\t' + toString(thisAmbig) + '\t';
                output += toString(thisHomoP) + '\t' + toString(num);
                
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
/**************************************************************************************************/

static DWORD WINAPI MySeqFastaSumThreadFunction(LPVOID lpParam){
    seqSumData* pDataArray;
    pDataArray = (seqSumData*)lpParam;
    
    try {
        ifstream in;
        pDataArray->m->openInputFile(pDataArray->summaryFile, in);
        
        //print header if you are process 0
        if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
            in.seekg(0);
            pDataArray->m->zapGremlins(in);
        }else { //this accounts for the difference in line endings.
            in.seekg(pDataArray->start-1); pDataArray->m->gobble(in);
        }
        
        string name;
        int start, end, length, ambigs, polymer, numReps;
        
        for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
            
            if (pDataArray->m->getControl_pressed()) { in.close(); pDataArray->count = 1; return 1; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> start >> end >> length >> ambigs >> polymer >> numReps; pDataArray->m->gobble(in);
            
            if (pDataArray->m->getDebug()) { pDataArray->m->mothurOut("[DEBUG]: " + name + "\t" + toString(start) + "\t" + toString(end) + "\t" + toString(length) + "\n"); }
            
            if (name != "") {
                
                if ((numReps == 1) && pDataArray->hasNameMap){
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator it = pDataArray->nameMap.find(name);
                    
                    if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + name + " is not in your name or count file, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->setControl_pressed(true); }
                    else { numReps = it->second; }
                }
                
                map<int, long long>::iterator it = pDataArray->startPosition.find(start);
                if (it == pDataArray->startPosition.end()) { pDataArray->startPosition[start] = numReps; } //first finding of this start position, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->endPosition.find(end);
                if (it == pDataArray->endPosition.end()) { pDataArray->endPosition[end] = numReps; } //first finding of this end position, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->seqLength.find(length);
                if (it == pDataArray->seqLength.end()) { pDataArray->seqLength[length] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts
;
                it = pDataArray->ambigBases.find(ambigs);
                if (it == pDataArray->ambigBases.end()) { pDataArray->ambigBases[ambigs] = numReps; } //first finding of this ambig, set count.
                else { it->second += numReps; } //add counts

                it = pDataArray->longHomoPolymer.find(polymer);
                if (it == pDataArray->longHomoPolymer.end()) { pDataArray->longHomoPolymer[polymer] = numReps; } //first finding of this homop, set count.
                else { it->second += numReps; } //add counts
                
                pDataArray->count++;
            }
        }
        in.close();
        
        return 0;
        
    }
    catch(exception& e) {
        pDataArray->m->errorOut(e, "SeqSummaryCommand", "MySeqFastaSumThreadFunction");
        exit(1);
    }
}
/**************************************************************************************************/
static DWORD WINAPI MySeqContigsSumThreadFunction(LPVOID lpParam){
    seqSumData* pDataArray;
    pDataArray = (seqSumData*)lpParam;
    
    try {
        ifstream in;
        pDataArray->m->openInputFile(pDataArray->filename, in);
        
        //print header if you are process 0
        if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
            in.seekg(0);
            pDataArray->m->zapGremlins(in); pDataArray->m->gobble(in);
            pDataArray->m->getline(in); pDataArray->m->gobble(in);
        }else { //this accounts for the difference in line endings.
            in.seekg(pDataArray->start-1); pDataArray->m->gobble(in);
        }
        
        string name;
        int length, OLength, thisOStart, thisOEnd, numMisMatches, numns, numReps; //Name	Length	Overlap_Length	Overlap_Start	Overlap_End	MisMatches	Num_Ns
        

        for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
            
            if (pDataArray->m->getControl_pressed()) { in.close(); pDataArray->count = 1; return 1; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numns; pDataArray->m->gobble(in);

            if (pDataArray->m->getDebug()) { pDataArray->m->mothurOut("[DEBUG]: " + name + "\t" + toString(thisOStart) + "\t" + toString(thisOEnd) + "\t" + toString(OLength) + "\n"); }
            
            numReps = 1;
            if (name != "") {
                if (pDataArray->hasNameMap){
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator it = pDataArray->nameMap.find(name);
                    
                    if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + name + " is not in your name or count file, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->setControl_pressed(true); }
                    else { numReps = it->second; }
                }
                
                map<int, long long>::iterator it = pDataArray->ostartPosition.find(thisOStart);
                if (it == pDataArray->ostartPosition.end()) { pDataArray->ostartPosition[thisOStart] = numReps; } //first finding of this start position, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->oendPosition.find(thisOEnd);
                if (it == pDataArray->oendPosition.end()) { pDataArray->oendPosition[thisOEnd] = numReps; } //first finding of this end position, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->seqLength.find(length);
                if (it == pDataArray->seqLength.end()) { pDataArray->seqLength[length] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->oseqLength.find(OLength);
                if (it == pDataArray->oseqLength.end()) { pDataArray->oseqLength[OLength] = numReps; } //first finding of this ambig, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->misMatches.find(numMisMatches);
                if (it == pDataArray->misMatches.end()) { pDataArray->misMatches[numMisMatches] = numReps; } //first finding of this ambig, set count.
                else { it->second += numReps; } //add counts
                
                it = pDataArray->numNs.find(numns);
                if (it == pDataArray->numNs.end()) { pDataArray->numNs[numns] = numReps; } //first finding of this homop, set count.
                else { it->second += numReps; } //add counts
                
                pDataArray->count++;
            }
        }
        in.close();
        
        return 0;
        
    }
    catch(exception& e) {
        pDataArray->m->errorOut(e, "SeqSummaryCommand", "MySeqContigsSumThreadFunction");
        exit(1);
    }
}
/**************************************************************************************************/
static DWORD WINAPI MySeqAlignSumThreadFunction(LPVOID lpParam){
    seqSumData* pDataArray;
    pDataArray = (seqSumData*)lpParam;
    
    try {
        ifstream in;
        pDataArray->m->openInputFile(pDataArray->filename, in);
        
        //print header if you are process 0
        if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
            in.seekg(0);
            pDataArray->m->zapGremlins(in); pDataArray->m->gobble(in);
            pDataArray->m->getline(in); pDataArray->m->gobble(in);
        }else { //this accounts for the difference in line endings.
            in.seekg(pDataArray->start-1); pDataArray->m->gobble(in);
        }
        
        string name, TemplateName, SearchMethod, AlignmentMethod;
        int length, TemplateLength,	 QueryStart,	QueryEnd,	TemplateStart,	TemplateEnd,	PairwiseAlignmentLength,	GapsInQuery,	GapsInTemplate,	LongestInsert, numReps;
        float SearchScore, SimBtwnQueryTemplate;  //QueryName	QueryLength	TemplateName	TemplateLength	SearchMethod	SearchScore	AlignmentMethod	QueryStart	QueryEnd	TemplateStart	TemplateEnd	PairwiseAlignmentLength	GapsInQuery	GapsInTemplate	LongestInsert	SimBtwnQuery&Template
        //checking for minScore, maxInsert, minSim
        
        for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
            
            if (pDataArray->m->getControl_pressed()) { in.close(); pDataArray->count = 1; return 1; }
            
            in >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; pDataArray->m->gobble(in);
            
            if (pDataArray->m->getDebug()) { pDataArray->m->mothurOut("[DEBUG]: " + name + "\t" + toString(length) + "\t" + toString(SearchScore) + "\t" + toString(SimBtwnQueryTemplate) + "\n"); }
            
            numReps = 1;
            if (name != "") {
                if (pDataArray->hasNameMap){
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator it = pDataArray->nameMap.find(name);
                    
                    if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + name + " is not in your name or count file, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->setControl_pressed(true); }
                    else { numReps = it->second; }
                }
                
                map<int, long long>::iterator it = pDataArray->seqLength.find(length);
                if (it == pDataArray->seqLength.end()) { pDataArray->seqLength[length] = numReps; } //first finding of this start position, set count.
                else { it->second += numReps; } //add counts
                
                map<float, long long>::iterator itFloat = pDataArray->sims.find(SimBtwnQueryTemplate);
                if (itFloat == pDataArray->sims.end()) { pDataArray->sims[SimBtwnQueryTemplate] = numReps; } //first finding of this end position, set count.
                else { itFloat->second += numReps; } //add counts
                
                itFloat = pDataArray->scores.find(SearchScore);
                if (itFloat == pDataArray->scores.end()) { pDataArray->scores[SearchScore] = numReps; } //first finding of this length, set count.
                else { itFloat->second += numReps; } //add counts
                
                it = pDataArray->inserts.find(LongestInsert);
                if (it == pDataArray->inserts.end()) { pDataArray->inserts[LongestInsert] = numReps; } //first finding of this ambig, set count.
                else { it->second += numReps; } //add counts
                
                pDataArray->count++;
            }
        }
        in.close();
        
        return 0;
        
    }
    catch(exception& e) {
        pDataArray->m->errorOut(e, "SeqSummaryCommand", "MySeqAlignSumThreadFunction");
        exit(1);
    }
}
/**************************************************************************************************/

#endif



#endif /* summary_hpp */
