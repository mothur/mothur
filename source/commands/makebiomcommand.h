#ifndef Mothur_makebiomcommand_h
#define Mothur_makebiomcommand_h

//
//  makebiomcommand.h
//  Mothur
//
//  Created by Sarah Westcott on 4/16/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "inputdata.h"


class MakeBiomCommand : public Command {
	
public:
	MakeBiomCommand(string);
	MakeBiomCommand();	
	~MakeBiomCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "make.biom";	}
	string getCommandCategory()		{ return "General";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://biom-format.org/documentation/biom_format.html, http://www.mothur.org/wiki/Make.biom"; }
	string getDescription()		{ return "creates a biom file"; }
    
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    
	string sharedfile, relabundfile, contaxonomyfile, metadatafile, groups, outputDir, format, label, referenceTax, picrustOtuFile, inputFileName, fileFormat;
	vector<string> outputNames, Groups, sampleMetadata;
	set<string> labels;
    
	bool abort, allLines, picrust;
    
    int getBiom(SharedRAbundVectors*&);
    int getBiom(SharedRAbundFloatVectors*&);
    vector<string> getMetaData(SharedRAbundVectors*&);
    vector<string> getMetaData(SharedRAbundFloatVectors*&);
    vector<string> parseTax(string tax, vector<string>& scores);
    int getSampleMetaData(SharedRAbundVectors*&);
    int getSampleMetaData(SharedRAbundFloatVectors*&);
    //for picrust
    int getGreenGenesOTUIDs(SharedRAbundVectors*&, map<string, string>&);
    int getGreenGenesOTUIDs(SharedRAbundFloatVectors*&, map<string, string>&);
    map<string, string> readGGOtuMap();
};


#endif
