//
//  fastqread.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/26/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "fastqread.h"

/*******************************************************************************/
FastqRead::FastqRead() {
    try {
        m = MothurOut::getInstance();
        format = "illumina1.8+"; name = ""; sequence = ""; scores.clear();
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
/*******************************************************************************/
FastqRead::FastqRead(string f) {
    try {
        m = MothurOut::getInstance();
        format = f; name = ""; sequence = ""; scores.clear();
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
/*******************************************************************************/

FastqRead::FastqRead(string f, string n, string s, vector<int> sc) {
    try {
        m = MothurOut::getInstance();
        format = f; name = n; sequence = s; scores = sc;
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
/*******************************************************************************/

FastqRead::FastqRead(ifstream& in, bool& ignore, string f) {
    try {
        m = MothurOut::getInstance();
        
        ignore = false;
        format = f;
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
        }
        
        //read sequence name
        string line = m->getline(in); m->gobble(in);
        vector<string> pieces = m->splitWhiteSpace(line);
        name = "";  if (pieces.size() != 0) { name = pieces[0]; }
        if (name == "") {  m->mothurOut("[WARNING]: Blank fasta name, ignoring read."); m->mothurOutEndLine(); ignore=true;  }
        else if (name[0] != '@') { m->mothurOut("[WARNING]: reading " + name + " expected a name with @ as a leading character, ignoring read."); m->mothurOutEndLine(); ignore=true; }
        else { name = name.substr(1); }
        
        //read sequence
        sequence = m->getline(in); m->gobble(in);
        if (sequence == "") {  m->mothurOut("[WARNING]: missing sequence for " + name + ", ignoring."); ignore=true; }
        
        //read sequence name
        line = m->getline(in); m->gobble(in);
        pieces = m->splitWhiteSpace(line);
        string name2 = "";  if (pieces.size() != 0) { name2 = pieces[0]; }
        if (name2 == "") {  m->mothurOut("[WARNING]: expected a name with + as a leading character, ignoring."); ignore=true; }
        else if (name2[0] != '+') { m->mothurOut("[WARNING]: reading " + name2 + " expected a name with + as a leading character, ignoring."); ignore=true; }
        else { name2 = name2.substr(1); if (name2 == "") { name2 = name; } }
        
        //read quality scores
        string quality = m->getline(in); m->gobble(in);
        if (quality == "") {  m->mothurOut("[WARNING]: missing quality for " + name2 + ", ignoring."); ignore=true; }
        
        //sanity check sequence length and number of quality scores match
        if (name2 != "") { if (name != name2) { m->mothurOut("[WARNING]: names do not match. read " + name + " for fasta and " + name2 + " for quality, ignoring."); ignore=true; } }
        if (quality.length() != sequence.length()) { m->mothurOut("[WARNING]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores, ignoring read."); ignore=true; }
        
        scores = convertQual(quality);
        m->checkName(name);
        
        if (m->debug) { m->mothurOut("[DEBUG]: " + name + " " + sequence + " " + quality + "\n"); }
    
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<int> FastqRead::convertQual(string qual) {
    try {
        vector<int> qualScores;
        bool negativeScores = false;
        
        for (int i = 0; i < qual.length(); i++) {
            
            int temp = 0;
            temp = int(qual[i]);
            if (format == "illumina") {
                temp -= 64; //char '@'
            }else if (format == "illumina1.8+") {
                temp -= int('!'); //char '!'
            }else if (format == "solexa") {
                temp = int(convertTable[temp]); //convert to sanger
                temp -= int('!'); //char '!'
            }else {
                temp -= int('!'); //char '!'
            }
            
            if (temp < -5) { negativeScores = true; }
            qualScores.push_back(temp);
        }
        
        if (negativeScores) { m->mothurOut("[ERROR]: finding negative quality scores, do you have the right format selected? http://en.wikipedia.org/wiki/FASTQ_format#Encoding \n");  m->control_pressed = true;  }
        
        return qualScores;
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "convertQual");
        exit(1);
    }
}
//**********************************************************************************************************************
Sequence FastqRead::getSequence() {
    try {
        Sequence temp(name, sequence);
        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "getSequence");
        exit(1);
    }
}
//**********************************************************************************************************************
QualityScores FastqRead::getQuality() {
    try {
        QualityScores temp(name, scores);
        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "getQuality");
        exit(1);
    }
}
/*******************************************************************************/
