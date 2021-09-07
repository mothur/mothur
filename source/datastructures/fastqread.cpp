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
            convertBackTable.push_back(((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499)));
        }
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
/*******************************************************************************/
FastqRead::FastqRead(Sequence s, QualityScores q) {
    try {
        m = MothurOut::getInstance(); format = "illumina1.8+";
        
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
            convertBackTable.push_back(((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499)));
        }
        
        if (s.getName() != q.getName()) { m->mothurOut("[ERROR]: sequence name does not match quality score name. Found sequence named " + s.getName() + " quality scores named " + q.getName() + " Cannot construct fastq object.\n"); m->setControl_pressed(true); }
        else {
            name = s.getName();
            comment = s.getComment();
            sequence = s.getUnaligned();
            scores = q.getScores();
            scoreString = convertQual(scores);
        }
        
        
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
/*******************************************************************************/
FastqRead::FastqRead(Sequence s, QualityScores q, string f) {
    try {
        m = MothurOut::getInstance(); format = f;
        
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
            convertBackTable.push_back(((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499)));
        }
        
        if (s.getName() != q.getName()) { m->mothurOut("[ERROR]: sequence name does not match quality score name. Found sequence named " + s.getName() + " quality scores named " + q.getName() + " Cannot construct fastq object.\n"); m->setControl_pressed(true); }
        else {
            name = s.getName();
            comment = s.getComment();
            sequence = s.getUnaligned();
            scores = q.getScores();
            scoreString = convertQual(scores);
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
            convertBackTable.push_back(((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499)));
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
            convertBackTable.push_back(((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499)));
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
        string line = util.getline(in); util.gobble(in);
        vector<string> pieces = util.splitWhiteSpace(line);
        name = "";  if (pieces.size() != 0) { name = pieces[0]; }
        if (name == "") {  m->mothurOut("[WARNING]: Blank fasta name, ignoring read.\n");  ignore=true;  }
        else if (name[0] != '@') { m->mothurOut("[WARNING]: reading " + name + " expected a name with @ as a leading character, ignoring read.\n");  ignore=true; }
        else { name = name.substr(1); }
        if (pieces.size() > 1) { pieces.erase(pieces.begin()); comment = util.getStringFromVector(pieces, " "); }
        
        //read sequence
        sequence = util.getline(in); util.gobble(in);
        if (sequence == "") {  m->mothurOut("[WARNING]: missing sequence for " + name + ", ignoring."); ignore=true; }

        //read sequence name
        line = util.getline(in); util.gobble(in);
        pieces = util.splitWhiteSpace(line);
        string name2 = "";  if (pieces.size() != 0) { name2 = pieces[0]; }
        if (name2 == "") {  m->mothurOut("[WARNING]: expected a name with + as a leading character, ignoring."); ignore=true; }
        else if (name2[0] != '+') { m->mothurOut("[WARNING]: reading " + name2 + " expected a name with + as a leading character, ignoring."); ignore=true; }
        else { name2 = name2.substr(1); if (name2 == "") { name2 = name; } }
        
        //read quality scores
        string quality = util.getline(in); util.gobble(in);
        if (quality == "") {  m->mothurOut("[WARNING]: missing quality for " + name2 + ", ignoring."); ignore=true; }
        
        //sanity check sequence length and number of quality scores match
        if (name2 != "") { if (name != name2) { m->mothurOut("[WARNING]: names do not match. read " + name + " for fasta and " + name2 + " for quality, ignoring."); ignore=true; } }
        if (quality.length() != sequence.length()) { m->mothurOut("[WARNING]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores, ignoring read."); ignore=true; }
        
        scoreString = quality;
        scores = convertQual(quality);
        util.checkName(name);
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + " " + sequence + " " + quality + "\n"); }
    
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
//**********************************************************************************************************************
#ifdef USE_BOOST
FastqRead::FastqRead(boost::iostreams::filtering_istream& in, bool& ignore, string f) {
    try {
        m = MothurOut::getInstance();
        
        ignore = false;
        format = f;
        
        if (in.eof()) { ignore = true; }
        else {
            //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
            for (int i = -64; i < 65; i++) {
                char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
                convertTable.push_back(temp);
            }
            
            //read sequence name
            string line = util.getline(in); util.gobble(in);
            vector<string> pieces = util.splitWhiteSpace(line);
            name = "";  if (pieces.size() != 0) { name = pieces[0];  }
            if (name == "") {  m->mothurOut("[WARNING]: Blank fasta name, ignoring read.\n");  ignore=true;  }
            else if (name[0] != '@') { m->mothurOut("[WARNING]: reading " + name + " expected a name with @ as a leading character, ignoring read.\n");  ignore=true; }
            else { name = name.substr(1); }
            if (pieces.size() > 1) { pieces.erase(pieces.begin()); comment = util.getStringFromVector(pieces, " "); }
            
            //read sequence
            sequence = util.getline(in); util.gobble(in);
            if (sequence == "") {  m->mothurOut("[WARNING]: missing sequence for " + name + ", ignoring."); ignore=true; }
            
            //read sequence name
            line = util.getline(in); util.gobble(in);
            pieces = util.splitWhiteSpace(line);
            string name2 = "";  if (pieces.size() != 0) { name2 = pieces[0]; }
            if (name2 == "") {  m->mothurOut("[WARNING]: expected a name with + as a leading character, ignoring."); ignore=true; }
            else if (name2[0] != '+') { m->mothurOut("[WARNING]: reading " + name2 + " expected a name with + as a leading character, ignoring."); ignore=true; }
            else { name2 = name2.substr(1); if (name2 == "") { name2 = name; } }
            
            //read quality scores
            string quality = util.getline(in); util.gobble(in);
            if (quality == "") {  m->mothurOut("[WARNING]: missing quality for " + name2 + ", ignoring."); ignore=true; }
            
            //sanity check sequence length and number of quality scores match
            if (name2 != "") { if (name != name2) { m->mothurOut("[WARNING]: names do not match. read " + name + " for fasta and " + name2 + " for quality, ignoring."); ignore=true; } }
            if (quality.length() != sequence.length()) { m->mothurOut("[WARNING]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores, ignoring read."); ignore=true; }
            
            scoreString = quality;
            scores = convertQual(quality);
            util.checkName(name);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + " " + sequence + " " + quality + "\n"); }
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "FastqRead");
        exit(1);
    }
}
#endif
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
                temp -= int('!'); //char '!' //33
            }else if (format == "solexa") {
                temp = int(convertTable[temp]); //convert to sanger
                temp -= int('!'); //char '!' //33
            }else {
                temp -= int('!'); //char '!' //33
            }
            
            if (temp < 0) { negativeScores = true; temp = 0; }
            qualScores.push_back(temp);
        }
        
        if (negativeScores) { m->mothurOut("[ERROR]: finding negative quality scores, do you have the right format selected? http://en.wikipedia.org/wiki/FASTQ_format#Encoding \n");  m->setControl_pressed(true);  }
        
        return qualScores;
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "convertQual");
        exit(1);
    }
}
//**********************************************************************************************************************
string FastqRead::convertQual(vector<int> qual) {
    try {
        string scoreString = "";
        
        for (int i = 0; i < qual.size(); i++) {
            int controlChar = int('!');
            if (format == "illumina") {  controlChar = int('@');  }

            int temp = qual[i] + controlChar;
            
            if (format == "solexa") { temp = convertBackTable[temp];  }
            
            char qualChar = (char) temp;
            
            scoreString += qualChar;
        }
        
        return scoreString;
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "convertQual");
        exit(1);
    }
}
//**********************************************************************************************************************
void FastqRead::setScores(vector<int> qual) {
    try {
        scoreString = "";
        scores = qual;
        
        for (int i = 0; i < qual.size(); i++) {
            int controlChar = int('!');
            if (format == "illumina") {  controlChar = int('@');  }

            int temp = qual[i] + controlChar;
            
            if (format == "solexa") { temp = convertBackTable[temp];  }
            
            char qualChar = (char) temp;
            
            scoreString += qualChar;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "setScores");
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
void FastqRead::printFastq(ostream& out) {
    try {
        out << "@" << name << " " << comment << endl;
        out << sequence << endl;
        out << "+" << endl;
        out << scoreString << endl;
        
    }
    catch(exception& e) {
        m->errorOut(e, "FastqRead", "printFastq");
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
