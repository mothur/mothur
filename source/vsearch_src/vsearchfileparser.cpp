//
//  vsearchfileparser.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/13/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "vsearchfileparser.h"
#include "deconvolutecommand.h"
#include "sequence.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"

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
        else {  m->mothurOut("[ERROR]: " + format + " is not a valid file format for the VsearchFileParser, quitting.\n"); m->control_pressed = true;  }
        
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "VsearchFileParser");
        exit(1);
    }			
}
/***********************************************************************/
string VsearchFileParser::getVsearchFile() {
    try {
        //Run unique.seqs on the data if a name or count file is not given
        if ((namefile == "") && (countfile == ""))  {  countfile = getNamesFile(fastafile);     }
        else if (namefile != "")                    {  counts = m->readNames(namefile);        }
        
        if (countfile != "") { CountTable countTable; countTable.readTable(countfile, false, false);  counts = countTable.getNameMap(); }
        
        if (m->control_pressed) {  return 0; }
        
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
        string vsearchFasta = m->getSimpleName(fastafile) + ".sorted.fasta.temp";
        
        vector<seqPriorityNode> seqs;
        map<string, int>::iterator it;
        
        ifstream in;
        m->openInputFile(inputFile, in);
        
        while (!in.eof()) {
            
            if (m->control_pressed) { in.close(); return vsearchFasta; }
            
            Sequence seq(in); m->gobble(in);
            
            it = counts.find(seq.getName());
            if (it == counts.end()) {
                m->mothurOut("[ERROR]: " + seq.getName() + " is not in your name or countfile, quitting.\n"); m->control_pressed = true;
            }else {
                seqPriorityNode temp(it->second, seq.getUnaligned(), it->first);
                seqs.push_back(temp);
            }
            
        }
        in.close();
        
        m->printVsearchFile(seqs, vsearchFasta);
        
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
        string nameFile = "";
        
        m->mothurOutEndLine(); m->mothurOut("No namesfile given, running unique.seqs command to generate one."); m->mothurOutEndLine(); m->mothurOutEndLine();
        
        //use unique.seqs to create new name and fastafile
        string inputString = "fasta=" + inputFile + ", format=count";
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine();
        m->mothurCalling = true;
        
        Command* uniqueCommand = new DeconvoluteCommand(inputString);
        uniqueCommand->execute();
        
        map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
        
        delete uniqueCommand;
        m->mothurCalling = false;
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        countfile = filenames["count"][0];
        fastafile = filenames["fasta"][0];
        
        return nameFile;
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "getNamesFile");
        exit(1);
    }
}
/*************************************************************************/
//S	1	275	*	*	*	*	*	GQY1XT001C44N8/ab=3677/	*
int VsearchFileParser::createListFile(string inputFile, string listFile, string sabundFile, string rabundFile, int numBins, string label){
    try {
        map<string, string>::iterator itName;
        if (format == "name") { counts.clear(); m->readNames(namefile, nameMap); }
        
        ifstream in;
        m->openInputFile(inputFile, in);
        
        ListVector list(numBins); list.setLabel(label);
        
        int clusterNumber;
        string seqName, recordType, length, percentIdentity, strand, notUsed1, notUsed2, compressedAlignment, repSequence;
        
        while(!in.eof()) {
            if (m->control_pressed) { break; }
            
            in >> recordType >> clusterNumber >> length >> percentIdentity >> strand >> notUsed1 >> notUsed2 >> compressedAlignment >> seqName >> repSequence; m->gobble(in);
            
            seqName = removeAbundances(seqName);
            
            if (format == "name") {
                itName = nameMap.find(seqName);
                if (itName == nameMap.end()) {  m->mothurOut("[ERROR]: " + seqName + " is not in your name file. Parsing error???\n"); m->control_pressed = true; }
                else{  seqName = itName->second;  }
            }
            
            string bin = list.get(clusterNumber);
            if (bin == "")  {   bin = seqName;          }
            else            {   bin += ',' + seqName;   }
            list.set(clusterNumber, bin);
            
        }
        in.close();
        
        ofstream out;
        m->openOutputFile(listFile,	out);
        list.printHeaders(out);
        
        if (countfile != "") {
            list.print(out, counts);
        }else {
            list.print(out);
            
            //print sabund and rabund
            ofstream sabund, rabund;
            m->openOutputFile(sabundFile, sabund);
            m->openOutputFile(rabundFile, rabund);
            
            list.getRAbundVector().print(rabund);  rabund.close();
            list.getSAbundVector().print(sabund);  sabund.close();
        }
        out.close();
        
        return 0;
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
        
        int pos = seqName.find_last_of("/", seqName.length()-2); //don't look at the last /
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
        m->openInputFile(logfile, in);
        
        string line;
        while(!in.eof()) {
            if (m->control_pressed) { break; }
            
            line = m->getline(in); m->gobble(in);
            
            int pos = line.find("Clusters:");
            if (pos != string::npos) {
                vector<string> pieces = m->splitWhiteSpace(line);
                if (pieces.size() > 1) { m->mothurConvert(pieces[1], numBins);  }
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

