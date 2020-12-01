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
#include "picrust.hpp"
#include "biomsimple.hpp"
#include "biomhdf5.hpp"


class MakeBiomCommand : public Command {
	
public:
	MakeBiomCommand(string);
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
    
	string sharedfile, relabundfile, contaxonomyfile, metadatafile, groups,  format, label, referenceTax, picrustOtuFile, inputFileName, fileFormat, output;
	vector<string> outputNames, Groups;
	set<string> labels;
    
	bool abort, allLines, picrust;
    
    void getBiom(SharedRAbundVectors*&, Picrust*, vector<Taxonomy>, vector<string>);
    void getBiom(SharedRAbundFloatVectors*&, Picrust*, vector<Taxonomy>, vector<string>);
     
    vector<string> getSampleMetaData(SharedRAbundVectors*&);
    vector<string> getSampleMetaData(SharedRAbundFloatVectors*&);
    
};


#endif
