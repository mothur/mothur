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

#define SFF_MAGIC   0x2e736666 /* ".sff" */
#define SFF_VERSION "\0\0\0\1"

/**********************************************************/
struct CommonHeader {
	unsigned int magicNumber;
	char* version;
	unsigned long int indexOffset;
	unsigned int indexLength;
	unsigned int numReads;
	unsigned short headerLength;
	unsigned short keyLength;
	unsigned short numFlowsPerRead;
	char flogramFormatCode;
	char* flowChars; //length depends on number flow reads
	char* keySequence; //length depends on key length
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
	char* name; //length depends on nameLength
};
/**********************************************************/
struct seqRead {
	vector<unsigned short> flowgram;
	vector<unsigned int> flowIndex;
	char* bases;
	vector<unsigned int> qualScores;
};
/**********************************************************/

class SffInfoCommand : public Command {
	
public:
	SffInfoCommand(string);
	~SffInfoCommand();
	int execute();
	void help();
	
private:
	string sffFilename, outputDir;
	vector<string> filenames, outputNames;
	bool abort;
	
	int extractSffInfo(string, string);
	int readCommonHeader(ifstream&, CommonHeader*&);
	int readHeader(ifstream&, Header*&);
	int readSeqData(ifstream&, seqRead*&, int, int);
	
	int printCommonHeader(ofstream&, CommonHeader*, bool); //bool is debug mode
	int printHeader(ofstream&, Header*, bool);
	int printSeqData(ofstream&, seqRead*, bool);
		
};

/**********************************************************/
 
#endif


