//
//  fastqread.h
//  Mothur
//
//  Created by Sarah Westcott on 1/26/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef Mothur_fastqread_h
#define Mothur_fastqread_h

#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"
#include "qualityscores.h"


/* This class is a representation of a fastqread.  If no format is given, defaults to illumina1.8+.
 
 @M00704:50:000000000-A3G0K:1:1101:15777:1541 2:N:0:0
 NCTCTACCAGGCCAAGCATAATGGGCGGGATCGTATCGAAGTAGCCTTGATGGGTAAGGTTGCCTGAGTTTCACAAGACAGATTACAGAGGTCGTCTATGCCCTGTCTCTTATACACATCTGACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAATATCGTCTAGGCCATGTGTGACGCTCGGTCTGGGCTTCACGAACAGGGGGTCCGCCATGTACCGCGCGCTC
 +
 #>>3AAFFFBAAFAGGFFFFGFHHHGGGG0EFGFHHFGHBFFGFDGHFGEGFFEBEGFCBFGFGFF2F4B3EGFHHHEHEHGHHH3FGHFG3BEEFHHHGGEGHFFHHEFGHHFHFHHF1B?FFD/AD/FC/<@D-.FGBF1<<<<<<GH0=GE<C<AD.0--:-;::900000900---.000./..;/;/9;//9/;;//....--;..//....9//9.;--/..---..--.-9/////.------.

 
 */

class FastqRead {
public:
    FastqRead();
    FastqRead(Sequence, QualityScores);
    FastqRead(Sequence, QualityScores, string);
    FastqRead(string f); 
    FastqRead(string f, string n, string s, vector<int> sc); 
    FastqRead(ifstream&, bool&, string f);
    #ifdef USE_BOOST
    FastqRead(boost::iostreams::filtering_istream&, bool&, string f);
    #endif
    ~FastqRead() {}
    
    void setFormat(string f) { format = f; }
    string getFormat() { return format; }
    string getName() { return name; }
    void setName(string n) { name = n; }
    string getSeq() { return sequence; }
    void setSeq(string s) { sequence = s; }
    vector<int> getScores() { return scores; }
    void setScores(vector<int> s) { scores = s;  }
    void printFastq(ostream&);
    
    Sequence getSequence();
    QualityScores getQuality();
    

private:
    MothurOut* m;
    vector<int> scores;
    string name, comment;
    string sequence;
    string scoreString;
    string format;
    vector<char> convertTable;
    vector<int> convertBackTable;
    
    vector<int> convertQual(string qual);
    string convertQual(vector<int>);
    
};


#endif
