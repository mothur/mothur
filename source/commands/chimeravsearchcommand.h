//
//  chimeravsearchcommand.h
//  Mothur
//
//  Created by Sarah Westcott on 6/16/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__chimeravsearchcommand__
#define __Mothur__chimeravsearchcommand__

#include "command.hpp"
#include "sequenceparser.h"
#include "counttable.h"
#include "sequencecountparser.h"

/**************************************************************************************************/
struct vsearchVariables {
    bool dups, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useXa, hasCount, useMindiffs;
    string abskew, minh, mindiv, xn, dn, mindiffs;
    
    vsearchVariables() {}
    void setBooleans(bool dps, bool Abskew, bool calns, bool MinH, bool Mindiv, bool Xn, bool Dn, bool mindif, bool hc) {
        useAbskew = Abskew;
        chimealns = calns;
        useMinH = MinH;
        useMindiv = Mindiv;
        useMindiffs = mindif;
        useXn = Xn;
        useDn = Dn;
        hasCount = hc;
        dups = dps;
    }
    
    void setVariables(string abske, string min, string mindi, string x, string d, string mind) {
        abskew = abske;
        minh = min;
        mindiv = mindi;
        mindiffs = mind;
        xn = x;
        dn = d;
    }
    
};
/**************************************************************************************************/
struct vsearchData {
    string fastafile;
    string dupsfile;
    string outputFName;
    string accnos, alns, formattedFastaFilename, templatefile, vsearchLocation;
    string driverAccnos, driverAlns, driverOutputFName;
    map<string, vector<string> > parsedFiles;
    map<string, vector<string> > seqs2RemoveByGroup;
    
    int count, numChimeras;
    vector<string> groups;
    vsearchVariables* vars;
    MothurOut* m;
    Utils util;
    
    vsearchData(){}
    vsearchData(map<string, vector<string> > g2f, string o, string uloc, string t, string file, string f, string n, string ac,  string al, string nc, vector<string> gr, vsearchVariables* vs) {
        fastafile = f;
        dupsfile = n;
        formattedFastaFilename = file;
        outputFName = o;
        templatefile = t;
        accnos = ac;
        alns = al;
        m = MothurOut::getInstance();
        groups = gr;
        count = 0;
        numChimeras = 0;
        vsearchLocation = uloc;
        vars = vs;
        driverAccnos = ac;
        driverAlns = al;
        driverOutputFName = o;
        parsedFiles = g2f;
    }
    void setDriverNames(string o, string al, string ac) {
        driverAccnos = ac;
        driverAlns = al;
        driverOutputFName = o;
    }
    
};

/***********************************************************/

class ChimeraVsearchCommand : public Command {
public:
    ChimeraVsearchCommand(string);
    ~ChimeraVsearchCommand() {}
    
    vector<string> setParameters();
    string getCommandName()			{ return "chimera.vsearch";		}
    string getCommandCategory()		{ return "Sequence Processing"; }
    
    string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);
    string getCitation() { return "vsearch by https://github.com/torognes/vsearch.\nhttp://www.mothur.org/wiki/Chimera.vsearch\n"; }
    string getDescription()		{ return "detect chimeric sequences"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, ucl, useMindiffs, hasCount, dups;
    string fastafile, templatefile, countfile, abskew, minh, mindiv, xn, dn, mindiffs, vsearchLocation;
    int processors;
    vsearchVariables* vars;
    vector<string> outputNames;
    
    string getCountFile(string&);
    int readFasta(string, map<string, string>&);
    int deconvoluteResults(string, string, string, long long&);
    int prepFile(string filename, string);
    int createProcessesGroups(map<string, vector<string> >& groups2Files, string outputFName, string filename, string accnos, string alns, string newCountFile, vector<string> groups, map<string, vector<string> >&);
};


#endif
