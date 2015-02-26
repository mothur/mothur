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
#include "oligos.h"

/**************************************************************************************************/

class SRACommand : public Command {
public:
    SRACommand(string);
    SRACommand();
    ~SRACommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.sra";			}
    string getCommandCategory()		{ return "Sequence Processing";		}
    
    string getOutputPattern(string);
    
	string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/Make.sra"; }
    string getDescription()		{ return "create a Sequence Read Archive / SRA"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    struct Grant {
        string grantId, grantTitle, grantAgency;
        Grant(string i, string a) : grantId(i), grantAgency(a), grantTitle("") {}
        Grant(string i, string a, string t) : grantId(i), grantAgency(a), grantTitle(t) {}
        Grant() : grantId(""), grantAgency(""), grantTitle("") {}
    };
    bool abort, isSFF, pairedOligos, setOligosParameter;
    int tdiffs, bdiffs, pdiffs, sdiffs, ldiffs, fileOption;
    
    string sfffile, fastqfile, outputDir, file, oligosfile, contactfile, inputfile, mimarksfile;
    string libStrategy, libSource, libSelection, libLayout, platform, instrumentModel, fileType, dataType, checkorient;
    string submissionName, lastName, firstName, email, centerName, centerType, ownership, description, website, orientation, packageType;
    string projectName, projectTitle, inputDir;
    vector<string> outputNames, Groups;
    vector<Grant> Grants;
    map<string, string> Group2Barcode;
    map<string, string> Group2Primer;
    vector<string> linkers;
    vector<string> spacers;
    map<string, string> Group2Organism;
    map<string, map<string, string> > mimarks;  //group -> <field -> valueForGroup> ex.  F003D001 -> <lat_lon -> 42.282026 -83.733850>
    set<string> uniqueNames;

    bool checkCasesInstrumentModels(string&);
    bool checkCasesPlatforms(string&);
    bool checkCasesLibStrategy(string&);
    bool checkCasesLibSource(string&);
    bool checkCasesLibSelection(string&);
    bool checkCasesDataType(string&);
    bool sanityCheckMiMarksGroups();
    int readFile(map<string, vector<string> >&);
    int readContactFile();
    int readMIMarksFile();
    int readOligos();
    int parseSffFile(map<string, vector<string> >&);
    int parseFastqFile(map<string, vector<string> >&);
    int checkGroups(map<string, vector<string> >&);
    int mapGroupToFile(map<string, vector<string> >&, vector<string>);
    int fixMap(map<string, vector<string> >&);
    int findFileOption();
    
};

/**************************************************************************************************/



#endif
