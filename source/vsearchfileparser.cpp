//
//  vsearchfileparser.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/13/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "vsearchfileparser.h"


/***********************************************************************/
VsearchFileParser::VsearchFileParser(){
    try {
        m = MothurOut::getInstance();
        fastafile = "";
        namefile = "";
        countfile = "";
        format = "count";
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "VsearchFileParser");
        exit(1);
    }
}
/***********************************************************************/
VsearchFileParser::VsearchFileParser(string f){
    try {
        m = MothurOut::getInstance();
        fastafile = f;
        namefile = "";
        countfile = "";
        format = "count";
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "VsearchFileParser");
        exit(1);
    }
}
/***********************************************************************/
VsearchFileParser::VsearchFileParser(string f, string nameOrCount, string forma) {
    try {
        m = MothurOut::getInstance();
        fastafile = f;
        namefile = "";
        countfile = "";
        format = forma;
        
        if (format == "name") { namefile = nameOrCount; }
        else if (format == "count") { countfile = nameOrCount; }
        else {  m->mothurOut("[ERROR]: " + format + " is not a valid file format for the VsearchFileParser, quitting.\n"); m->setControl_pressed(true);  }
        
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "VsearchFileParser");
        exit(1);
    }			
}
/***********************************************************************/
string VsearchFileParser::getVsearchFile() {
    try {
        Utils util;
        if (fastafile == "") { m->mothurOut("[ERROR]: no fasta file given, cannot continue.\n"); m->setControl_pressed(true);  }
        
        //Run unique.seqs on the data if a name or count file is not given
        if ((namefile == "") && (countfile == ""))  {  getNamesFile(fastafile);                  }
        else if (namefile != "")                    {  counts = util.readNames(namefile);        }
        
        if (countfile != "") { CountTable countTable; countTable.readTable(countfile, false, false);  counts = countTable.getNameMap(); }
        
        if (m->getControl_pressed()) {  return 0; }
        
        //Remove gap characters from each sequence if needed
        //Append the number of sequences that each unique sequence represents to the end of the fasta file name
        //Sorts by abundance
        string vsearchFastafile = createVsearchFasta(fastafile);
        
        return vsearchFastafile;
        
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "getVsearchFile");
        exit(1);
    }
}
/**********************************************************************/

string VsearchFileParser::createVsearchFasta(string inputFile){
    try {
        Utils util;
        string vsearchFasta = util.getSimpleName(fastafile) + ".sorted.fasta.temp";
        
        vector<seqPriorityNode> seqs;
        map<string, int>::iterator it;
        
        ifstream in;
        util.openInputFile(inputFile, in);
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); return vsearchFasta; }
            
            Sequence seq(in); util.gobble(in);
            
            it = counts.find(seq.getName());
            if (it == counts.end()) {
                m->mothurOut("[ERROR]: " + seq.getName() + " is not in your name or countfile, quitting.\n"); m->setControl_pressed(true);
            }else {
                seqPriorityNode temp(it->second, seq.getUnaligned(), it->first);
                seqs.push_back(temp);
            }
            
        }
        in.close();
        
        util.printVsearchFile(seqs, vsearchFasta, ";size=", ";");
        
        return vsearchFasta;
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "createVsearchFasta");
        exit(1);
    }
}
/*************************************************************************/

string VsearchFileParser::getNamesFile(string& inputFile){
    try {
        
        m->mothurOut("\nNo namesfile given, running unique.seqs command to generate one.\n\n");
        
        //use unique.seqs to create new name and fastafile
        string inputString = "fasta=" + inputFile + ", format=count";
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: unique.seqs(" + inputString + ")\n");
        Command* uniqueCommand = new DeconvoluteCommand(inputString);
        uniqueCommand->execute();
        
        map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
        
        delete uniqueCommand;
        m->mothurOut("/******************************************/\n");
        
        countfile = filenames["count"][0];
        fastafile = filenames["fasta"][0];
        
        return countfile;
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "getNamesFile");
        exit(1);
    }
}
/*************************************************************************/
//S	1	275	*	*	*	*	*	GQY1XT001C44N8/ab=3677/	*
ListVector VsearchFileParser::createListFile(string inputFile, int numBins, string label, map<string, int>& ct){
    try {
        Utils util;
        map<string, string>::iterator itName;
        if (format == "name") { counts.clear(); util.readNames(namefile, nameMap); }
        
        ifstream in;
        util.openInputFile(inputFile, in);
        
        ListVector list(numBins); list.setLabel(label);
        
        int clusterNumber;
        string seqName, recordType, length, percentIdentity, strand, notUsed1, notUsed2, compressedAlignment, repSequence;
        
        while(!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            in >> recordType >> clusterNumber >> length >> percentIdentity >> strand >> notUsed1 >> notUsed2 >> compressedAlignment >> seqName >> repSequence; util.gobble(in);
            
            if (recordType != "S") {
                
                seqName = removeAbundances(seqName);
                
                if (format == "name") {
                    itName = nameMap.find(seqName);
                    if (itName == nameMap.end()) {  m->mothurOut("[ERROR]: " + seqName + " is not in your name file. Parsing error???\n"); m->setControl_pressed(true); }
                    else{  seqName = itName->second;  }
                }
                
                string bin = list.get(clusterNumber);
                if (bin == "")  {   bin = seqName;          }
                else            {   bin += ',' + seqName;   }
                list.set(clusterNumber, bin);
                
            }
            
        }
        in.close();
        ct = counts;
        
        return list;
        
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "createListFile");
        exit(1);
    }
}
/*************************************************************************/
//GQY1XT001C44N8/ab=3677/	*
string VsearchFileParser::removeAbundances(string seqName){
    try {
        
        int pos = seqName.find_last_of(";", seqName.length()-2); //don't look at the last /
        if (pos != string::npos) { seqName = seqName.substr(0, pos); }
        
        return seqName;
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "removeAbundances");
        exit(1);
    }
}
/*************************************************************************/
//GQY1XT001C44N8/ab=3677/	*
int VsearchFileParser::getNumBins(string logfile){
    try {
        
        int numBins = 0;
        
        ifstream in;
        Utils util; util.openInputFile(logfile, in);
        
        string line;
        while(!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            line = util.getline(in); util.gobble(in);
            
            int pos = line.find("Clusters:");
            if (pos != string::npos) {
                vector<string> pieces = util.splitWhiteSpace(line);
                if (pieces.size() > 1) { util.mothurConvert(pieces[1], numBins);  }
                break;
            }
        }
        in.close();
        
        return numBins;
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "getNumBins");
        exit(1);
    }
}
/***********************************************************************/

