//
//  sffread.hpp
//  Mothur
//
//  Created by Sarah Westcott on 6/9/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef sffread_hpp
#define sffread_hpp

#include "mothurout.h"
#include "endiannessmacros.h"
#include "utils.hpp"


/* This class is a representation of a sff read.
 
 https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=format#sff
 
 readHeader:
 
 read_header_length         uint16_t

 name_length                uint16_t

 number_of_bases            uint32_t

 clip_qual_left             uint16_t

 clip_qual_right            uint16_t

 clip_adapter_left          uint16_t

 clip_adapter_right         uint16_t

 name                       char[name_length]

 eight_byte_padding         uint8_t[*]
 
 
 readInfo:
 
 flowgram_values            uint*_t[number_of_flows]

 flow_index_per_base        uint8_t[number_of_bases]

 bases                      char[number_of_bases]

 quality_scores             uint8_t[number_of_bases]

 eight_byte_padding         uint8_t[*]

 */

class SffRead {
    
    public:
    
    SffRead(ifstream&, int);
    SffRead(int num);
    ~SffRead();
    
    bool readSff(ifstream& in);
    bool isOkay() { return good; }
    
    void printFasta(ofstream& out, bool trim);
    void printQuality(ofstream& out, bool trim);
    void printFlow(ofstream& out);
    void printSff(ofstream& out);
    void printSffTxt(ofstream& out);
    
    //read header info
    string getName()        { return name;       }
    string getTimeStamp()   { return timestamp;  }
    string getRegion()      { return region;     }
    string getXY()          { return xy;         }
    unsigned short getHeaderLength()        { return headerLength;       }
    unsigned short getNameLength()          { return nameLength;         }
    unsigned short getClipQualLeft()        { return clipQualLeft;       }
    unsigned short getClipQualRight()       { return clipQualRight;      }
    unsigned short getClipAdapterLeft()     { return clipAdapterLeft;    }
    unsigned short getClipAdapterRight()    { return clipAdapterRight;   }
    unsigned int   getNumBases()            { return numBases;           }
    
    //read info
    vector<unsigned short> getFlowgrams()   { return flowgram;           }
    vector<unsigned int> getFlowIndex()     { return flowIndex;          }
    vector<unsigned int> getQualScores()    { return qualScores;         }
    string getBases()                       { return bases;              }
    
    void setName(string n)                      { name = n;       }
    void setTimeStamp(string n)                 { timestamp = n;  }
    void setRegion(string n)                    { region = n;     }
    void setXY(string n)                        { xy = n;         }
    void setHeaderLength(unsigned short n)      { headerLength = n;       }
    void setNameLength(unsigned short n)        { nameLength = n;         }
    void setClipQualLeft(unsigned short n)      { clipQualLeft = n;       }
    void setClipQualRight(unsigned short n)     { clipQualRight = n;      }
    void setClipAdapterLeft(unsigned short n)   { clipAdapterLeft = n;    }
    void setClipAdapterRight(unsigned short n)  { clipAdapterRight = n;   }
    void setNumBases(unsigned int n)            { numBases = n;           }
    
    //read info
    void setFlowgrams(vector<unsigned short> n)   { flowgram = n;           }
    void setFlowIndex(vector<unsigned int> n)     { flowIndex = n;          }
    void setQualScores(vector<unsigned int> n)    { qualScores = n;         }
    void setBases(string n)                       { bases = n;              }
    
private:
    MothurOut* m;
    vector<char*> entireRead;
    
    //header fields
    unsigned short headerLength;
    unsigned short nameLength;
    unsigned int numBases;
    unsigned short clipQualLeft;
    unsigned short clipQualRight;
    unsigned short clipAdapterLeft;
    unsigned short clipAdapterRight;
    string name; //length depends on nameLength
    string timestamp;
    string region;
    string xy;
    
    //readFields
    vector<unsigned short> flowgram;
    vector<unsigned int> flowIndex;
    string bases;
    vector<unsigned int> qualScores;
    
    int numFlows, padSize1, padSize2;
    unsigned long long size;
    bool good;
    
    void printSffTxtHeader(ofstream& out);
    int decodeName(string&, string&, string&, string);
    bool sanityCheck();

};


#endif /* sffread_hpp */
