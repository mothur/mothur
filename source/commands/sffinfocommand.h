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
#include "groupmap.h"
#include "oligos.h"
#include "trimoligos.h"
#include "sffread.hpp"
#include "sffheader.hpp"

/**********************************************************/

class SffInfoCommand : public Command {
	
public:
	SffInfoCommand(string);
	~SffInfoCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sff.info";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Sffinfo"; }
	string getDescription()		{ return "extract sequences reads from a .sff file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string sffFilename, sfftxtFilename,  accnosName, currentFileName, oligosfile, noMatchFile, groupfile;
	vector<string> outputNames;
	bool abort, fasta, qual, trim, flow, sfftxt, hasAccnos, hasOligos, hasGroup, reorient, pairedOligos;
	int mycount, split, numBarcodes, numFPrimers, numLinkers, numSpacers, numRPrimers, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, numNoMatch;
    unordered_set<string> seqNames;
    GroupMap* groupMap;
    map<string, int> GroupToFile;
    vector<vector<int> > numSplitReads;
    vector<vector<string> > filehandles;
    vector<vector<string> > filehandlesHeaders;
    Oligos* oligosObject;
    
	//extract sff file functions
	int extractSffInfo(string, string, string); //main function
	void assignToSample(SffRead*&, TrimOligos*&, TrimOligos*&);
    bool readOligos(string oligosFile);
    bool readGroup(string oligosFile);
    
    //assign read to sample when splitting
    int findGroup(SffRead*& read, int& barcode, int& primer, TrimOligos*&, TrimOligos*&);
    int findGroup(SffRead*&, int& barcode, int& primer, string);
    
    //common header functions
	int adjustCommonHeader(SffCommonHeader*&);
    
	//parsesfftxt file functions
    int parseSffTxt();
	int parseHeaderLineToInt(ifstream&);
    unsigned short parseHeaderLineToShort(ifstream& file);
	vector<unsigned short> parseHeaderLineToFloatVector(ifstream&, int);
	vector<unsigned int> parseHeaderLineToIntVector(ifstream&, int);
	string parseHeaderLineToString(ifstream&);
};

/**********************************************************/
 
#endif


