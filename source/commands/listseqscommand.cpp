/*
 *  listseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "listseqscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"
#include "counttable.h"
#include "fastqread.h"

//**********************************************************************************************************************
vector<string> ListSeqsCommand::setParameters(){	
	try {
        CommandParameter pfastq("fastq", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pfastq);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pfasta);
        CommandParameter pqfile("qfile", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pqfile);
		CommandParameter pname("name", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false); parameters.push_back(palignreport);
        CommandParameter pcontigsreport("contigsreport", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false); parameters.push_back(pcontigsreport);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["accnos"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The list.seqs command reads a fasta, name, group, count, list, taxonomy, fastq, qfile, alignreport or contigsreport file and outputs a .accnos file containing sequence names.\n";
		helpString += "The list.seqs command parameters are fasta, name, group, count, list, taxonomy, fastq, contigsreport and alignreport.  You must provide one of these parameters.\n";
		helpString += "The list.seqs command should be in the following format: list.seqs(fasta=yourFasta).\n";
		helpString += "Example list.seqs(fasta=amazon.fasta).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "accnos") {  pattern = "[filename],accnos"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ListSeqsCommand::ListSeqsCommand(string option) : Command()  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastafiles = validParameter.validFiles(parameters, "fasta");
            if (fastafiles.size() != 0) {
                if (fastafiles[0] == "not open") { abort = true; }
                else { current->setFastaFile(fastafiles[0]); }
            }
            
			namefiles = validParameter.validFiles(parameters, "name");
            if (namefiles.size() != 0) {
                if (namefiles[0] == "not open") { abort = true; }
                else { current->setNameFile(namefiles[0]); }
            }
			
			groupfiles = validParameter.validFiles(parameters, "group");
            if (groupfiles.size() != 0) {
                if (groupfiles[0] == "not open") { abort = true; }
                else { current->setGroupFile(groupfiles[0]); }
            }
			
			alignfiles = validParameter.validFiles(parameters, "alignreport");
            if (alignfiles.size() != 0) {
                if (alignfiles[0] == "not open") { abort = true; }
            }
			
            contigsreportfiles = validParameter.validFiles(parameters, "contigsreport");
            if (contigsreportfiles.size() != 0) {
                if (contigsreportfiles[0] == "not open") { abort = true; }
                else { current->setContigsReportFile(contigsreportfiles[0]); }
            }
			
			listfiles = validParameter.validFiles(parameters, "list");
            if (listfiles.size() != 0) {
                if (listfiles[0] == "not open") { abort = true; }
                else { current->setListFile(listfiles[0]); }
            }
			
			taxfiles = validParameter.validFiles(parameters, "taxonomy");
            if (taxfiles.size() != 0) {
                if (taxfiles[0] == "not open") { abort = true; }
                else { current->setTaxonomyFile(taxfiles[0]); }
            }
            
            countfiles = validParameter.validFiles(parameters, "count");
            if (countfiles.size() != 0) {
                if (countfiles[0] == "not open") { abort = true; }
                else { current->setCountFile(countfiles[0]); }
            }
            
            fastqfiles = validParameter.validFiles(parameters, "fastq");
            if (fastqfiles.size() != 0) {
                if (fastqfiles[0] == "not open") { abort = true; }
            }
			
            qualityfiles = validParameter.validFiles(parameters, "qfile");
            if (qualityfiles.size() != 0) {
                if (qualityfiles[0] == "not open") { abort = true; }
                else { current->setQualFile(qualityfiles[0]); }
            }
            
			if ((qualityfiles.size() == 0) && (fastqfiles.size() == 0) && (countfiles.size() == 0) && (fastafiles.size() == 0) && (namefiles.size() == 0) && (listfiles.size() == 0) && (groupfiles.size() == 0) && (alignfiles.size() == 0) && (taxfiles.size() == 0) && (contigsreportfiles.size() == 0))  { m->mothurOut("You must provide a file.\n"); abort = true; }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "ListSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void addName(bool empty, string name, set<string>& names, set<string>& newNames) {
    if (empty) { newNames.insert(name); } //for first file or single file
    else {
        if (names.count(name) != 0) { newNames.insert(name); } //present in files so far so add to newNames
    }
}
//**********************************************************************************************************************
void readFastq(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }
        Utils util;
        
        while(!in.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            bool ignore;
            FastqRead fread(in, ignore, "illumina1.8+"); util.gobble(in);
            
            if (!ignore) { addName(empty, fread.getName(), names, newNames); }
        }
        
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readFastq");
        exit(1);
    }
}
//**********************************************************************************************************************
void readQual(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }
        Utils util;
        
        while(!in.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            QualityScores currSeq(in); util.gobble(in);
            
            if (currSeq.getName() != "") { addName(empty, currSeq.getName(), names, newNames); }
        }
        
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readQual");
        exit(1);
    }
}
//**********************************************************************************************************************
void readFasta(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }
        Utils util;
        
        while(!in.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            Sequence currSeq(in); util.gobble(in);
            
            if (currSeq.getName() != "") { addName(empty, currSeq.getName(), names, newNames); }
        }
        
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readFasta");
        exit(1);
    }
}
//**********************************************************************************************************************
void readList(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }
        
        Utils util; string tag = "Otu"; string readHeaders = ""; //Tells mothur to try and read headers from the file
        
        if(!in.eof()){
            ListVector list(in, readHeaders, tag); //read in list vector
            
            //for each bin
            for (int i = 0; i < list.getNumBins(); i++) {
                if (m->getControl_pressed()) { break; }
                
                string bin = list.get(i);
                vector<string> binnames; util.splitAtComma(bin, binnames);
                
                for (int j = 0; j < binnames.size(); j++) { addName(empty, binnames[j], names, newNames); }
            }
        }
        
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readList");
        exit(1);
    }
}
//**********************************************************************************************************************
void readNameTaxGroup(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }

        Utils util; string name;
        
        while(!in.eof()){
        
            if (m->getControl_pressed()) { break; }

            in >> name; util.getline(in); util.gobble(in);
            
            addName(empty, name, names, newNames);
        }
        
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readNameTaxGroup");
        exit(1);
    }
}
//**********************************************************************************************************************
void readCount(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }
        
        CountTable ct; ct.readTable(in, false, false);
        
        if (m->getControl_pressed()) { return; }
        
        vector<string> cnames = ct.getNamesOfSeqs();
        
        for (int j = 0; j < cnames.size(); j++) { addName(empty, cnames[j], names, newNames); }
        
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readCount");
        exit(1);
    }
}
//**********************************************************************************************************************
void readAlignContigs(set<string>& names, ifstream& in, MothurOut*& m){
    try {
        set<string> newNames;
        bool empty = true;
        if (names.size() != 0) { empty=false; }
        string name;
        
        Utils util; util.getline(in);  util.gobble(in);
        
        while(!in.eof()){
            if (m->getControl_pressed()) { break; }

            in >> name; util.getline(in); util.gobble(in);
            
            addName(empty, name, names, newNames);
        }
 
        names = newNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "readAlignContigs");
        exit(1);
    }
}
//**********************************************************************************************************************

int ListSeqsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        set<string> names;
        
		//read functions fill names vector
		if (fastafiles.size() != 0)		    {   process(fastafiles, names, &readFasta);	                }
        if (qualityfiles.size() != 0)       {   process(qualityfiles, names, &readQual);                }
        if (fastqfiles.size() != 0)	        {	process(fastqfiles, names, "fastq", &readFastq);	    }
		if (namefiles.size() != 0)	        {	process(namefiles, names, &readNameTaxGroup);           }
		if (groupfiles.size() != 0)	        {	process(groupfiles, names, &readNameTaxGroup);          }
        if (taxfiles.size() != 0)           {   process(taxfiles, names, &readNameTaxGroup);            }
		if (alignfiles.size() != 0)	        {   process(alignfiles, names, &readAlignContigs);          }
        if (contigsreportfiles.size() != 0) {   process(contigsreportfiles, names, &readAlignContigs);  }
		if (listfiles.size() != 0)	        {	process(listfiles, names, &readList);                   }
        if (countfiles.size() != 0)	        {	process(countfiles, names, &readCount);                 }
        
		if (m->getControl_pressed()) { outputTypes.clear();  return 0; }
		
		if (outputdir == "") {  outputdir += util.hasPath(inputFileName);  }
		
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
		string outputFileName = getOutputFileName("accnos", variables);

        util.printAccnos(outputFileName, names);
        
		outputNames.push_back(outputFileName); outputTypes["accnos"].push_back(outputFileName);
		
		if (m->getControl_pressed()) { outputTypes.clear();  util.mothurRemove(outputFileName); return 0; }
		
		current->setAccnosFile(outputFileName);
		
		m->mothurOut("\nOutput File Names: \n" + outputFileName + "\n\n");
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		return 0;		
	}catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void ListSeqsCommand::process(vector<string> files, set<string>& names, string isFastq, void f(set<string>&, ifstream&, MothurOut*&)){
    try {
        Utils util;
        
        //determine if the files are compressed. If so,
        
        for (int i = 0; i < files.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            inputFileName = files[i];
            
            ifstream in; util.openInputFile(inputFileName, in);
            
            f(names, in, m);
            
            in.close();
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************
void ListSeqsCommand::process(vector<string> files, set<string>& names, void f(set<string>&, ifstream&, MothurOut*&)){
    try {
        Utils util;
        for (int i = 0; i < files.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            inputFileName = files[i];
            
            ifstream in; util.openInputFile(inputFileName, in);
            
            f(names, in, m);
            
            in.close();
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************
