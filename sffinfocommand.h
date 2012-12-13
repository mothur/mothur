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
	unsigned long long indexOffset;
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
	string timestamp;
	string region;
	string xy;
	
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
	SffInfoCommand();
	~SffInfoCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sffinfo";					}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Sffinfo"; }
	string getDescription()		{ return "extract sequences reads from a .sff file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string sffFilename, sfftxtFilename, outputDir, accnosName, currentFileName, oligosfile, noMatchFile;
	vector<string> filenames, outputNames, accnosFileNames, oligosFileNames;
	bool abort, fasta, qual, trim, flow, sfftxt, hasAccnos, hasOligos;
	int mycount, split, numFPrimers, numLinkers, numSpacers, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs;
	set<string> seqNames;
    map<string, int> barcodes;
    map<string, int> primers;
    vector<string> linker, spacer, primerNameVector, barcodeNameVector, revPrimer;
    vector<vector<int> > numSplitReads;
    vector<vector<string> > filehandles, filehandlesHeaders;
    
	//extract sff file functions
	int extractSffInfo(string, string, string);
	int readCommonHeader(ifstream&, CommonHeader&);
	//int readHeader(ifstream&, Header&);
	int readSeqData(ifstream&, seqRead&, int, Header&);
	int decodeName(string&, string&, string&, string);
    bool readOligos(string oligosFile);
	
	int printCommonHeader(ofstream&, CommonHeader&); 
	int printHeader(ofstream&, Header&);
	int printSffTxtSeqData(ofstream&, seqRead&, Header&);
	int printFlowSeqData(ofstream&, seqRead&, Header&);
	int printFastaSeqData(ofstream&, seqRead&, Header&);
	int printQualSeqData(ofstream&, seqRead&, Header&);
	int readAccnosFile(string);
	int parseSffTxt();
	bool sanityCheck(Header&, seqRead&);
    int adjustCommonHeader(CommonHeader);
    int findGroup(Header header, seqRead read, int& barcode, int& primer);
    string reverseOligo(string oligo);
    
	//parsesfftxt file functions
	int parseHeaderLineToInt(ifstream&);
	vector<unsigned short> parseHeaderLineToFloatVector(ifstream&, int);
	vector<unsigned int> parseHeaderLineToIntVector(ifstream&, int);
	string parseHeaderLineToString(ifstream&);
};

/**********************************************************/
 
#endif


