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

/***********************************************************************/
VsearchFileParser::VsearchFileParser(string f){
    try {
        m = MothurOut::getInstance();
        fastafile = f;
        namefile = "";
        countfile = "";
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "VsearchFileParser");
        exit(1);
    }
}
/***********************************************************************/
VsearchFileParser::VsearchFileParser(string f, string nameOrCount, string format) {
    try {
        m = MothurOut::getInstance();
        fastafile = f;
        namefile = "";
        countfile = "";
        
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
        map<string, int> nameMap;
        if ((namefile == "") && (countfile == ""))  {  countfile = getNamesFile(fastafile);     }
        else if (namefile != "")                    {  nameMap = m->readNames(namefile);        }
        
        if (countfile != "") { CountTable countTable; countTable.readTable(countfile, false, false);  nameMap = countTable.getNameMap(); }
        
        if (m->control_pressed) {  return 0; }
        
        //Remove gap characters from each sequence if needed
        //Append the number of sequences that each unique sequence represents to the end of the fasta file name
        //Sorts by abundance
        string vsearchFastafile = createVsearchFasta(fastafile, nameMap);  nameMap.clear();
        
        return vsearchFastafile;
        
    }
    catch(exception& e) {
        m->errorOut(e, "VsearchFileParser", "getVsearchFile");
        exit(1);
    }
}
/**********************************************************************/

string VsearchFileParser::createVsearchFasta(string inputFile, map<string, int>& nameMap){
    try {
        string vsearchFasta = m->getSimpleName(fastafile) + ".sorted.fasta.temp";
        
        vector<seqPriorityNode> seqs;
        map<string, int>::iterator it;
        
        ifstream in;
        m->openInputFile(inputFile, in);
        
        while (!in.eof()) {
            
            if (m->control_pressed) { in.close(); return vsearchFasta; }
            
            Sequence seq(in); m->gobble(in);
            
            it = nameMap.find(seq.getName());
            if (it == nameMap.end()) {
                m->mothurOut("[ERROR]: " + seq.getName() + " is not in your name or countfile, quitting.\n"); m->control_pressed = true;
            }else {
                seqPriorityNode temp(it->second, seq.getUnaligned(), it->first);
                seqs.push_back(temp);
                nameMap.erase(it);
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
/***********************************************************************/

