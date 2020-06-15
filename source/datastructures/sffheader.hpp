//
//  sffheader.hpp
//  Mothur
//
//  Created by Sarah Westcott on 6/10/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef sffheader_hpp
#define sffheader_hpp

#include "mothurout.h"
#include "sequence.hpp"
#include "qualityscores.h"
#include "endiannessmacros.h"


/* This class is a representation of a sff common header.
 
 https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=format#sff
 
 /*
  magic_number               uint32_t
 
 version                    char[4]
 
 index_offset               uint64_t
 
 index_length               uint32_t
 
 number_of_reads            uint32_t
 
 header_length              uint16_t
 
 key_length                 uint16_t
 
 number_of_flows_per_read   uint16_t
 
 flowgram_format_code       uint8_t
 
 flow_chars                 char[number_of_flows_per_read]
 
 key_sequence               char[key_length]
 
 eight_byte_padding         uint8_t[*]
 
  
 The magic_number field value is 0x2E736666, the uint32_t encoding of the string ".sff"
 The version number corresponding to this proposal is 0001, or the byte array "\0\0\0\1".
 The index_offset and index_length fields are the offset and length of an optional index of the reads in the SFF file. If no index is included in the file, both fields must be 0.
 The number_of_reads field should be set to the number of reads stored in the file.
 The header_length field should be the total number of bytes required by this set of header fields, and should be equal to "31 + number_of_flows_per_read + key_length" rounded up to the next value divisible by 8.
 The key_length and key_sequence fields should be set to the length and nucleotide bases of the key sequence used for these reads.
 Note: The key_sequence field is not null-terminated.
 The number_of_flows_per_read should be set to the number of flows for each of the reads in the file.
 The flowgram_format_code should be set to the format used to encode each of the flowgram values for each read.
 Note: Currently, only one flowgram format has been adopted, so this value should be set to 1.
 The flowgram format code 1 stores each value as a uint16_t, where the floating point flowgram value is encoded as "(int) round(value * 100.0)", and decoded as "(storedvalue * 1.0 / 100.0)".
 The flow_chars should be set to the array of nucleotide bases ('A', 'C', 'G' or 'T') that correspond to the nucleotides used for each flow of each read. The length of the array should equal number_of_flows_per_read.
 Note: The flow_chars field is not null-terminated.
 If any eight_byte_padding bytes exist in the section, they should have a byte value of 0.
 If an index is included in the file, the index_offset and index_length values in the common header should point to the section of the file containing the index. To support different indexing methods, the index section should begin with the following two fields:
 
   */
/**********************************************************/

class SffCommonHeader {
    public:
    
    SffCommonHeader();
    SffCommonHeader(ifstream&);
    ~SffCommonHeader();
    
    bool read(ifstream& in);
    void printSFFTxt(ofstream&);
    void printSampleCommonHeader(ofstream& out, int numReads); 

    unsigned short getHeaderLength()        { return headerLength;       }
    unsigned short getKeyLength()           { return keyLength;          }
    unsigned short getNumFlows()            { return numFlows;           }
    unsigned int getMagicNumber()           { return magicNumber;        }
    unsigned long long getIndexLength()     { return indexLength;        }
    unsigned short getIndexOffset()         { return indexOffset;        }
    unsigned int getNumReads()              { return numReads;           }
    string getVersion()                     { return version;            }
    int getFlowgramFormat()                 { return flogramFormatCode;  }
    string getFlows()                       { return flowChars;          }
    string getKeySequence()                 { return keySequence;        }
    
    void print(ofstream&);
    
private:
    MothurOut* m;

    vector<char*> entireHeader;
    
    int padSize;
    unsigned int magicNumber;
    string version;
    unsigned long long indexOffset;
    unsigned int indexLength;
    unsigned int numReads;
    unsigned short headerLength;
    unsigned short keyLength;
    unsigned short numFlows;
    int flogramFormatCode;
    string flowChars; //length depends on number flow reads
    string keySequence; //length depends on key length
    
};

#endif /* sffheader_hpp */
