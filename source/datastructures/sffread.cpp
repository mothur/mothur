//
//  sffread.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/9/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "sffread.hpp"

//***************************************************************************************
SffRead::~SffRead(){
    if (entireRead.size() != 0) {
        for (int i = 0; i < entireRead.size(); i++) { delete[] entireRead[i]; }
        entireRead.clear();
    }
}
//***************************************************************************************
SffRead::SffRead(int num){
    try {
        m = MothurOut::getInstance();
        numFlows = num;
        padSize1 = 0; padSize2 = 0;
        
        headerLength=0; nameLength=0; numBases=0; clipQualLeft=0; clipQualRight=0; clipAdapterLeft=0; clipAdapterRight=0; bases = ""; name = "";
        
        good = false;
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "SffRead");
        exit(1);
    }
}
//***************************************************************************************
SffRead::SffRead(ifstream& in, int num){
    try {
        m = MothurOut::getInstance();
        numFlows = num;
        padSize1 = 0; padSize2 = 0;
        
        headerLength=0; nameLength=0; numBases=0; clipQualLeft=0; clipQualRight=0; clipAdapterLeft=0; clipAdapterRight=0; bases = ""; name = "";
        
        good = readSff(in);
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "SffRead");
        exit(1);
    }
}
//***************************************************************************************
bool SffRead::readSff(ifstream& in) {
    try {
        bool goodRead = true;
        
        if (!in.eof()) {
            
            unsigned long long startSpotInFile = in.tellg();
            /*****************************************/
            //read header length
            char* readHeaderLength = new char[2];
            in.read(&(*readHeaderLength), 2);
            entireRead.push_back(readHeaderLength);
            headerLength = be_int2(*(unsigned short *)(readHeaderLength));
            
            //read name length
            char* readNameLength = new char [2];
            in.read(&(*readNameLength), 2);
            entireRead.push_back(readNameLength);
            nameLength = be_int2(*(unsigned short *)(readNameLength));
            
            //read num bases
            char* readNumBases = new char [4];
            in.read(&(*readNumBases), 4);
            entireRead.push_back(readNumBases);
            numBases =  be_int4(*(unsigned int *)(readNumBases));
            
            //read clip qual left
            char* rclipQualLeft = new char [2];
            in.read(&(*rclipQualLeft), 2);
            entireRead.push_back(rclipQualLeft);
            clipQualLeft = be_int2(*(unsigned short *)(rclipQualLeft));
            
            //read clip qual right
            char* rclipQualRight = new char [2];
            in.read(&(*rclipQualRight), 2);
            entireRead.push_back(rclipQualRight);
            clipQualRight = be_int2(*(unsigned short *)(rclipQualRight));
            
            //read clipAdapterLeft
            char* rclipAdapterLeft = new char  [2];
            in.read(&(*rclipAdapterLeft), 2);
            entireRead.push_back(rclipAdapterLeft);
            clipAdapterLeft = be_int2(*(unsigned short *)(rclipAdapterLeft));
            
            //read clipAdapterRight
            char* rclipAdapterRight = new char  [2];
            in.read(&(*rclipAdapterRight), 2);
            entireRead.push_back(rclipAdapterRight);
            clipAdapterRight = be_int2(*(unsigned short *)(rclipAdapterRight));
            
            //read name
            char* readName = new char[nameLength];
            in.read(&(*readName), nameLength);
            for (int i = 0; i < nameLength; i++) { name += readName[i]; }
            entireRead.push_back(readName);
            
            //extract info from name
            decodeName(timestamp, region, xy, name);
            
            /* Pad to 8 chars */
            unsigned long long spotInFile = in.tellg();
            unsigned long long spot = (spotInFile + 7)& ~7;
            char* padding = new char[spot-spotInFile];
            entireRead.push_back(padding);
            padSize1 = spot-spotInFile;
            in.seekg(spot);
        
            /*****************************************/
            //sequence read
            
            //read flowgram
            flowgram.resize(numFlows);
            char* flows = new char[numFlows*2];
            int count = 0;
            for (int i = 0; i < numFlows; i++) {
                char rflowgram [2];
                in.read(rflowgram, 2);
                flows[count] = rflowgram[0]; count++;
                flows[count] = rflowgram[1]; count++;
                flowgram[i] = be_int2(*(unsigned short *)(&rflowgram));
            }
            entireRead.push_back(flows);
            
            //read flowIndex
            flowIndex.resize(numBases);
            char* flowI = new char[numBases];
            count = 0;
            for (int i = 0; i < numBases; i++) {
                char flowINdex[1];
                in.read(flowINdex, 1);
                flowI[count] = flowINdex[0]; count++;
                flowIndex[i] = be_int1(*(unsigned char *)(&flowINdex));
            }
            entireRead.push_back(flowI);
            
            //read bases
            char* readBases = new char[numBases];
            in.read(&(*readBases), numBases);
            for (int i = 0; i < numBases; i++) { bases += readBases[i]; }
            entireRead.push_back(readBases);
            
            //read qual scores
            qualScores.resize(numBases, 0);
            char* scores = new char[numBases];
            count = 0;
            for (int i = 0; i < numBases; i++) {
                char score[1];
                in.read(score, 1);
                scores[count] = score[0]; count++;
                qualScores[i] = be_int1(*(unsigned char *)(&score));
            }
            entireRead.push_back(scores);
            
            /* Pad to 8 chars */
            spotInFile = in.tellg();
            spot = (spotInFile + 7)& ~7;
            
            size = spot - startSpotInFile;
            
            char* padding2 = new char[spot-spotInFile];
            entireRead.push_back(padding2);
            padSize2 = spot-spotInFile;

            goodRead = sanityCheck();
            
            //ensure good reset
            in.seekg(spot);
            
        }else { size = 0; goodRead = false;}
        
        good = goodRead;
        
        return goodRead;
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "read");
        exit(1);
    }
}
//***************************************************************************************

int SffRead::decodeName(string& timestamp, string& region, string& xy, string name) {
    try {
        Utils util;
        
        if (name.length() >= 6) {
            string time = name.substr(0, 6);
            unsigned int timeNum = util.fromBase36(time);
            
            int q1 = timeNum / 60;
            int sec = timeNum - 60 * q1;
            int q2 = q1 / 60;
            int minute = q1 - 60 * q2;
            int q3 = q2 / 24;
            int hr = q2 - 24 * q3;
            int q4 = q3 / 32;
            int day = q3 - 32 * q4;
            int q5 = q4 / 13;
            int mon = q4 - 13 * q5;
            int year = 2000 + q5;
            
            timestamp = toString(year) + "_" + toString(mon) + "_" + toString(day) + "_" + toString(hr) + "_" + toString(minute) + "_" + toString(sec);
        }
        
        if (name.length() >= 9) {
            region = name.substr(7, 2);
            
            string xyNum = name.substr(9);
            unsigned int myXy = util.fromBase36(xyNum);
            int x = myXy >> 12;
            int y = myXy & 4095;
            
            xy = toString(x) + "_" + toString(y);
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "decodeName");
        exit(1);
    }
}
//*********************************************************************************************
bool SffRead::sanityCheck() {
    try {
        bool okay = true;
        string message = "[WARNING]: Your sff file may be corrupted! Sequence: " + getName() + "\n";

        int readLength = getBases().length();
        int qualLength = getQualScores().size();
        unsigned short clipLeft = getClipQualLeft();
        unsigned short clipRight = getClipQualRight();
        
        if (clipLeft > readLength) {
            okay = false; message += "Clip Qual Left = " + toString(clipLeft) + ", but we only read " + toString(readLength) + " bases.\n";
        }
        if (clipRight > readLength) {
            okay = false; message += "Clip Qual Right = " + toString(clipRight) + ", but we only read " + toString(readLength) + " bases.\n";
        }
        if (clipLeft > qualLength) {
            okay = false; message += "Clip Qual Left = " + toString(clipLeft) + ", but we only read " + toString(qualLength) + " quality scores.\n";
        }
        if (clipRight > qualLength) {
            okay = false; message += "Clip Qual Right = " + toString(clipRight) + ", but we only read " + toString(qualLength) + " quality scores.\n";
        }
        
        if (!okay) { m->mothurOut(message+"\n");  }
        
        return okay;
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "sanityCheck");
        exit(1);
    }
}
//**********************************************************************************
void SffRead::printSff(ofstream& out) {
    try {
        if (entireRead.size() != 0) {
            out.write(entireRead[0], 2); //write header length
            out.write(entireRead[1], 2); //write name length
            out.write(entireRead[2], 4); //write num bases
            
            out.write(entireRead[3], 2); //write clip qual left
            out.write(entireRead[4], 2); //write clip qual right
            out.write(entireRead[5], 2); //write clipAdapterLeft
            out.write(entireRead[6], 2); //write clipAdapterRight
             
            out.write(entireRead[7], nameLength); //write name
            out.write(entireRead[8], padSize1); //write pad1
            
            out.write(entireRead[9], numFlows*2); //write flowgram
            out.write(entireRead[10], numBases); //write flowIndex
            out.write(entireRead[11], numBases); //write bases
            out.write(entireRead[12], numBases); //write qual scores
            
            out.write(entireRead[13], padSize2); //write pad2
        }
        else {
            m->mothurOut("[ERROR]: cannot print sff, did not read it, skipping.\n");
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "printSff");
        exit(1);
    }
}
//***************************************************************************************
void SffRead::printFasta(ofstream& out, bool trim) {
    try {
        string seq = bases;
        
        if (trim) {
            if(clipQualRight < clipQualLeft){
                if (clipQualRight == 0) { //don't trim right
                    seq = seq.substr(clipQualLeft-1);
                }else {
                    seq = "NNNN";
                }
            }
            else if((clipQualRight != 0) && ((clipQualRight-clipQualLeft) >= 0)){
                seq = seq.substr((clipQualLeft-1), (clipQualRight-clipQualLeft+1));
            }
            else {
                seq = seq.substr(clipQualLeft-1);
            }
        }else{
           
            int endValue = clipQualRight;
            //make the bases you want to clip lowercase and the bases you want to keep upper case
            if(endValue == 0){    endValue = seq.length();    }
            for (int i = 0; i < (clipQualLeft-1); i++) { seq[i] = tolower(seq[i]);  }
            for (int i = (clipQualLeft-1); i < (endValue-1); i++)  {   seq[i] = toupper(seq[i]);  }
            for (int i = (endValue-1); i < seq.length(); i++) {   seq[i] = tolower(seq[i]);  }
        }
        
        out << ">" << name  << " xy=" << xy << endl;
        out << seq << endl;

    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "printFasta");
        exit(1);
    }
}
//**********************************************************************************************************************
void SffRead::printQuality(ofstream& out, bool trim) {
    try {
        
        if (trim) {
            if(clipQualRight < clipQualLeft){
                if (clipQualRight == 0) { //don't trim right
                    out << ">" << name << " xy=" << xy << " length=" << (qualScores.size()-clipQualLeft) << endl;
                    for (int i = (clipQualLeft-1); i < qualScores.size(); i++) {   out << qualScores[i] << '\t';    }
                }else {
                    out << ">" << name << " xy=" << xy << endl;
                    out << "0\t0\t0\t0";
                }
            }
            else if((clipQualRight != 0) && ((clipQualRight-clipQualLeft) >= 0)){
                out << ">" << name << " xy=" << xy << " length=" << (clipQualRight-clipQualLeft+1) << endl;
                for (int i = (clipQualLeft-1); i < (clipQualRight); i++) {   out << qualScores[i] << '\t';    }
            }
            else{
                out << ">" << name << " xy=" << xy << " length=" << (clipQualRight-clipQualLeft) << endl;
                for (int i = (clipQualLeft-1); i < qualScores.size(); i++) {   out << qualScores[i] << '\t';    }
            }
        }else{
            out << ">" << name << " xy=" << xy << " length=" << qualScores.size() << endl;
            for (int i = 0; i < qualScores.size(); i++) {   out << qualScores[i] << '\t';  }
        }
        
        out << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "SffInfoCommand", "printQuality");
        exit(1);
    }
}

//**********************************************************************************************************************
void SffRead::printFlow(ofstream& out) {
    try {
        
        int endValue = clipQualRight;
        if (clipQualRight == 0) {
            endValue = flowIndex.size();
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + " has clipQualRight=0.\n"); }
        }
        if(endValue > clipQualLeft){
            
            int rightIndex = 0;
            for (int i = 0; i < endValue; i++) {  rightIndex +=  flowIndex[i];     }
            
            out << name << ' ' << rightIndex;
            for (int i = 0; i < flowgram.size(); i++) { out << setprecision(2) << ' ' << (flowgram[i]/(float)100);  }
            out << endl;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "printFlow");
        exit(1);
    }
}
//**********************************************************************************************************************
void SffRead::printSffTxt(ofstream& out) {
    try {
        printSffTxtHeader(out);
        
        out << "Flowgram: ";
        for (int i = 0; i < flowgram.size(); i++) { out << setprecision(2) << (flowgram[i]/(float)100) << '\t';  }
        
        out << endl <<  "Flow Indexes: ";
        int sum = 0;
        for (int i = 0; i < flowIndex.size(); i++) {  sum +=  flowIndex[i];  out << sum << '\t'; }
        
        //make the bases you want to clip lowercase and the bases you want to keep upper case
        int endValue = clipQualRight;
        if(endValue == 0){    endValue = bases.length();    }
        for (int i = 0; i < (clipQualLeft-1); i++) { bases[i] = tolower(bases[i]); }
        for (int i = (clipQualLeft-1); i < (endValue-1); i++) {   bases[i] = toupper(bases[i]);  }
        for (int i = (endValue-1); i < bases.length(); i++) {   bases[i] = tolower(bases[i]);  }
        
        out << endl <<  "Bases: " << bases << endl << "Quality Scores: ";
        for (int i = 0; i < qualScores.size(); i++) {   out << qualScores[i] << '\t';  }
        out << endl << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "printSffTxt");
        exit(1);
    }
}
//**********************************************************************************************************************
void SffRead::printSffTxtHeader(ofstream& out) {
    try {
        
        out << ">" << name << endl;
        out << "Run Prefix: " << timestamp << endl;
        out << "Region #:  " << region << endl;
        out << "XY Location: " << xy << endl << endl;
        
        out << "Run Name:  " << endl;
        out << "Analysis Name:  " << endl;
        out << "Full Path: " << endl << endl;
        
        out << "Read Len: " << headerLength << endl;
        out << "Name Length: " << nameLength << endl;
        out << "# of Bases: " << numBases << endl;
        out << "Clip Qual Left: " << clipQualLeft << endl;
        out << "Clip Qual Right: " << clipQualRight << endl;
        out << "Clip Adap Left: " << clipAdapterLeft << endl;
        out << "Clip Adap Right: " << clipAdapterRight << endl << endl;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SffRead", "printSffTxtHeader");
        exit(1);
    }
}
//****************************************************************************************
