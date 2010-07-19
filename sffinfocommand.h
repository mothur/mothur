#ifndef SFFINFOCOMMAND_H
#define SFFINFOCOMMAND_H

/*
 *  sffinfocommand.h
 *  Mothur
 *
 *  Created by westcott on 7/7/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"


/**********************************************************/
struct CommonHeader {
	unsigned int magicNumber;
	string version;
	unsigned long int indexOffset;
	unsigned int indexLength;
	unsigned int numReads;
	unsigned short headerLength;
	unsigned short keyLength;
	unsigned short numFlowsPerRead;
	int flogramFormatCode;
	string flowChars; //length depends on number flow reads
	string keySequence; //length depends on key length
	
	CommonHeader(){ magicNumber=0; indexOffset=0; indexLength=0; numReads=0; headerLength=0; keyLength=0; numFlowsPerRead=0; flogramFormatCode='s'; } 
	~CommonHeader() { }
};
/**********************************************************/
struct Header {
	unsigned short headerLength;
	unsigned short nameLength;
	unsigned int numBases;
	unsigned short clipQualLeft;
	unsigned short clipQualRight;
	unsigned short clipAdapterLeft;
	unsigned short clipAdapterRight;
	string name; //length depends on nameLength
	
	Header() { headerLength=0; nameLength=0; numBases=0; clipQualLeft=0; clipQualRight=0; clipAdapterLeft=0; clipAdapterRight=0; }
	~Header() { } 
};
/**********************************************************/
struct seqRead {
	vector<unsigned short> flowgram;
	vector<unsigned int> flowIndex;
	string bases;
	vector<unsigned int> qualScores;
	
	seqRead() { } 
	~seqRead() { } 
};
/**********************************************************/

class SffInfoCommand : public Command {
	
public:
	SffInfoCommand(string);
	~SffInfoCommand();
	int execute();
	void help();
	
private:
	string sffFilename, outputDir, accnosName;
	vector<string> filenames, outputNames, accnosFileNames;
	bool abort, fasta, qual, trim, flow, sfftxt, hasAccnos;
	set<string> seqNames;
	
	int extractSffInfo(string, string);
	int readCommonHeader(ifstream&, CommonHeader&);
	int readHeader(ifstream&, Header&);
	int readSeqData(ifstream&, seqRead&, int, int);
	
	int printCommonHeader(ofstream&, CommonHeader&); 
	int printHeader(ofstream&, Header&);
	int printSffTxtSeqData(ofstream&, seqRead&);
	int printFlowSeqData(ofstream&, seqRead&, Header&);
	int printFastaSeqData(ofstream&, seqRead&, Header&);
	int printQualSeqData(ofstream&, seqRead&, Header&);
	int readAccnosFile(string);
		
};

/**********************************************************/
 
#endif


