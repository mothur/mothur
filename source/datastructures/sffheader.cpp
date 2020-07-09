//
//  sffheader.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/10/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "sffheader.hpp"

//***************************************************************************************
SffCommonHeader::~SffCommonHeader(){
    if (entireHeader.size() != 0) {
        for (int i = 0; i < entireHeader.size(); i++) { delete[] entireHeader[i]; }
        entireHeader.clear();
    }
}
//***************************************************************************************
SffCommonHeader::SffCommonHeader(){
    try {
        m = MothurOut::getInstance();
        padSize = 0;
        magicNumber=0; indexOffset=0; indexLength=0; numReads=0; headerLength=0; keyLength=0; numFlows=0; flogramFormatCode='s';
    }
    catch(exception& e) {
        m->errorOut(e, "SffCommonHeader", "SffCommonHeader");
        exit(1);
    }
}
//***************************************************************************************
SffCommonHeader::SffCommonHeader(ifstream& in){
    try {
        m = MothurOut::getInstance();
        read(in);
        padSize = 0;
        magicNumber=0; indexOffset=0; indexLength=0; numReads=0; headerLength=0; keyLength=0; numFlows=0; flogramFormatCode='s';
    }
    catch(exception& e) {
        m->errorOut(e, "SffCommonHeader", "SffCommonHeader");
        exit(1);
    }
}
//**********************************************************************************
bool SffCommonHeader::read(ifstream& in){
    try {
        
        bool goodHeader = true;
        
        if (!in.eof()) {
            
            unsigned long long startSpotInFile = in.tellg();
            
            //read magic number
            char* magic = new char[4];
            in.read(&(*magic), 4);
            magicNumber = be_int4(*(unsigned int *)(magic));
            entireHeader.push_back(magic);
            
            //read version
            char* cversion = new char[4];
            in.read(&(*cversion), 4);
            entireHeader.push_back(cversion);
            version = "";
            for (int i = 0; i < 4; i++) {  version += toString((int)(cversion[i]));  }
    
            //read offset - ignored in print
            char buffer2 [8];
            in.read(buffer2, 8);
            indexOffset =  be_int8(*(unsigned long long *)(&buffer2));
            
            //read index length - ignored in print
            char buffer3 [4];
            in.read(buffer3, 4);
            indexLength =  be_int4(*(unsigned int *)(&buffer3));
            
            //read num reads - ignored in print and set to samples numReads
            char rnumReads[4];
            in.read(rnumReads, 4);
            numReads = be_int4(*(unsigned int *)(&rnumReads));
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: numReads = " + toString(numReads) + "\n"); }
                
            //read header length
            char* hlength = new char  [2];
            in.read(&(*hlength), 2);
            entireHeader.push_back(hlength);
            headerLength =  be_int2(*(unsigned short *)(hlength));
                    
            //read key length
            char* klength = new char  [2];
            in.read(&(*klength), 2);
            entireHeader.push_back(klength);
            keyLength = be_int2(*(unsigned short *)(klength));
            
            //read number of flow reads
            char* nflows = new char  [2];
            in.read(&(*nflows), 2);
            entireHeader.push_back(nflows);
            numFlows =  be_int2(*(unsigned short *)(nflows));
                
            //read format code
            char* fcode = new char[1];
            in.read(&(*fcode), 1);
            entireHeader.push_back(fcode);
            flogramFormatCode = (int)(fcode[0]);
            
            //read flow chars
            char* tempBuffer = new char[numFlows];
            in.read(&(*tempBuffer), numFlows);
            flowChars = tempBuffer;
            if (flowChars.length() > numFlows) { flowChars = flowChars.substr(0, numFlows);  }
            entireHeader.push_back(tempBuffer);
            
            //read key
            char* tempBuffer2 = new char[keyLength];
            in.read(&(*tempBuffer2), keyLength);
            keySequence = tempBuffer2;
            if (keySequence.length() > keyLength) { keySequence = keySequence.substr(0, keyLength);  }
            entireHeader.push_back(tempBuffer2);
            
            /* Pad to 8 chars */
            unsigned long long spotInFile = in.tellg();
            unsigned long long spot = (spotInFile + 7)& ~7;  // ~ inverts
            
            char* padding = new char[spot-spotInFile];
            entireHeader.push_back(padding);
            
            padSize = spot-spotInFile;
            
            //ensure good reset
            in.seekg(spot);
            
            //check magic number and version
            if (magicNumber != 779314790) { m->mothurOut("[ERROR]: Magic Number is not correct, not a valid .sff file\n");  goodHeader = false; }
            if (version != "0001") { m->mothurOut("[ERROR]: Version is not supported, only support version 0001.\n");  goodHeader = false;  }
            
        }else{ m->mothurOut("Error reading sff common header.\n");  goodHeader = false; }
        
        return goodHeader;
    }
    catch(exception& e) {
        m->errorOut(e, "SffCommonHeader", "read");
        exit(1);
    }
}
//****************************************************************************************
void SffCommonHeader::printSampleCommonHeader(ofstream& out, int numReads){
    try {
        //magic number
        out.write(entireHeader[0], 4);
        
        //version
        out.write(entireHeader[1], 4);
        
        //offset - read and discard, we will set it to 0
        long long offset = 0;
        char offsetBuffer[8];
        offsetBuffer[0] = (offset >> 56) & 0xFF;
        offsetBuffer[1] = (offset >> 48) & 0xFF;
        offsetBuffer[2] = (offset >> 40) & 0xFF;
        offsetBuffer[3] = (offset >> 32) & 0xFF;
        offsetBuffer[4] = (offset >> 24) & 0xFF;
        offsetBuffer[5] = (offset >> 16) & 0xFF;
        offsetBuffer[6] = (offset >> 8) & 0xFF;
        offsetBuffer[7] = offset & 0xFF; //index = 15
        
        out.write(offsetBuffer, 8);
        
        offset = 0;
        char readIndexLength[4];
        readIndexLength[0] = (offset >> 24) & 0xFF;
        readIndexLength[1] = (offset >> 16) & 0xFF;
        readIndexLength[2] = (offset >> 8) & 0xFF;
        readIndexLength[3] = offset & 0xFF; //index = 19
                
        out.write(readIndexLength, 4);
        
        //change num reads
        char numSampleReads[4];
        numSampleReads[0] = (numReads >> 24) & 0xFF;
        numSampleReads[1] = (numReads >> 16) & 0xFF;
        numSampleReads[2] = (numReads >> 8) & 0xFF;
        numSampleReads[3] = numReads & 0xFF; //index = 23
        
        out.write(numSampleReads, 4);
        
        //read header length
        out.write(entireHeader[2], 2);
            
        //read key length
        out.write(entireHeader[3], 2);
            
        //read number of flow reads
        out.write(entireHeader[4], 2);
            
        //read format code
        out.write(entireHeader[5], 1);
            
        //read flow chars
        out.write(entireHeader[6], numFlows);
        
        //read key
        out.write(entireHeader[7], keyLength);
        
        /* Pad to 8 chars */
        out.write(entireHeader[8], padSize);
        
    }
    catch(exception& e) {
        m->errorOut(e, "SffInfoCommand", "printSampleCommonHeader");
        exit(1);
    }
}
//***********************************************************************************
void SffCommonHeader::printSFFTxt(ofstream& out) {
    try {
        out << "Common Header:\nMagic Number: " << magicNumber << endl;
        out << "Version: " << version << endl;
        out << "Index Offset: " << indexOffset << endl;
        out << "Index Length: " << indexLength << endl;
        out << "Number of Reads: " << numReads << endl;
        out << "Header Length: " << headerLength << endl;
        out << "Key Length: " << keyLength << endl;
        out << "Number of Flows: " << numFlows << endl;
        out << "Format Code: " << flogramFormatCode << endl;
        out << "Flow Chars: " << flowChars << endl;
        out << "Key Sequence: " << keySequence << endl << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "SffCommonHeader", "printSFFTxt");
        exit(1);
    }
}
//***********************************************************************************
