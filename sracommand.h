//
//  sracommand.h
//  Mothur
//
//  Created by SarahsWork on 10/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_sracommand_h
#define Mothur_sracommand_h

#include "command.hpp"
#include "trimoligos.h"

/**************************************************************************************************/

class SRACommand : public Command {
public:
    SRACommand(string);
    SRACommand();
    ~SRACommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "sra";			}
    string getCommandCategory()		{ return "Sequence Processing";		}
    
    string getOutputPattern(string);
    
	string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/sra"; }
    string getDescription()		{ return "create a Sequence Read Archive / SRA"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, isSFF, pairedOligos;
    int tdiffs, bdiffs, pdiffs, sdiffs, ldiffs;
    string sfffile, fastqfile, outputDir, groupfile, file, oligosfile, contactfile, inputfile;
    string libStrategy, libSource, libSelection, libLayout, platform, instrumentModel, fileType;
    string submissionName, lastName, firstName, email, centerName, centerType, description;
    vector<string> outputNames, Groups, revPrimer;
    vector<string> primerNameVector;
    vector<string> barcodeNameVector;
    map<string, string> Group2Barcode;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
	map<string, int> barcodes;
	map<string, int> primers;
    vector<string>  linker;
    vector<string>  spacer;

    bool checkCasesInstrumentModels(string&);
    bool checkCasesPlatforms(string&);
    bool checkCasesLibStrategy(string&);
    bool checkCasesLibSource(string&);
    bool checkCasesLibSelection(string&);
    int readFile(map<string, vector<string> >&);
    int readContactFile();
    int readOligos();
    int parseSffFile(map<string, vector<string> >&);
    int parseFastqFile(map<string, vector<string> >&);
    int checkGroups(map<string, vector<string> >&);
    int mapGroupToFile(map<string, vector<string> >&, vector<string>);
    string reverseOligo(string oligo);
    
};

/**************************************************************************************************/



#endif
