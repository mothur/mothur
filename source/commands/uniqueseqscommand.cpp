/*
 *  deconvolute.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "uniqueseqscommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> UniqueSeqsCommand::setParameters(){
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","count",false,false,true); parameters.push_back(pcount);
        CommandParameter pformat("format", "Multiple", "count-name", "count", "", "", "","",false,false, true); parameters.push_back(pformat);
        CommandParameter poutput("output", "Multiple", "count-name", "count", "", "", "","",false,false, true); parameters.push_back(poutput);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
       
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UniqueSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string UniqueSeqsCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The unique.seqs command reads a fastafile and creates a name or count file.\n";
		helpString += "The unique.seqs command parameters are fasta, name, count and format.  fasta is required, unless there is a valid current fasta file.\n";
        helpString += "The name parameter is used to provide an existing name file associated with the fasta file. \n";
        helpString += "The count parameter is used to provide an existing count file associated with the fasta file. \n";
        helpString += "The format parameter is used to indicate what type of file you want outputted.  Choices are name and count, default=count unless name file used then default=name.\n";
		helpString += "The unique.seqs command should be in the following format: \n";
		helpString += "unique.seqs(fasta=yourFastaFile) \n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "UniqueSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string UniqueSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],unique,[extension]"; } 
        else if (type == "name") {  pattern = "[filename],names-[filename],[tag],names"; } 
        else if (type == "count") {  pattern = "[filename],count_table-[filename],[tag],count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "UniqueSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
/**************************************************************************************/
UniqueSeqsCommand::UniqueSeqsCommand(string option) : Command()  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {
				fastafile = current->getFastaFile();
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required.\n");  abort = true; }
			}else { current->setFastaFile(fastafile); }
			
            if (outputdir == ""){ outputdir += util.hasPath(fastafile);  }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){	namefile = "";	}
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a unique.seqs command you must enter ONLY ONE of the following: count or name.\n");  abort = true; }
			
            //allow format parameter to have two names - format or output
            format = validParameter.valid(parameters, "format");
            if(format == "not found"){
                format = validParameter.valid(parameters, "output");
                if(format == "not found"){
                    if (namefile != "") { format = "name";    }
                    else { format = "count";                  }
                }
            }
            
            if ((format != "name") && (format != "count")) {
                m->mothurOut(format + " is not a valid format option. Options are count or name.");
                if (countfile == "") { m->mothurOut("I will use count.\n"); format = "count"; }
                else {  m->mothurOut("I will use count.\n"); format = "count"; }
            }
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "UniqueSeqsCommand", "UniqueSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
int UniqueSeqsCommand::execute() {
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (countfile != "")        {  processCount(countfile); }
        else if (namefile != "")    {  processName(namefile);   }
        else { //no existing duplicate data
            if (format == "count")      {  processCount("");        }
            else                        {  processName("");         }
        }
		if (m->getControl_pressed()) { return 0; }
		
		m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "UniqueSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
string UniqueSeqsCommand::processCount(string countfile) { //countfile can be blank, indicating no countfile provided
    try {
        if (format == "name") {
            if (countfile == "") {  return (processName("")); }
            else {
                CountTable ct; ct.readTable(countfile, true, false);
                map<string, int> nameMap = ct.getNameMap();
                
                string newNamesFile = createNewNameFile(countfile, nameMap);
                
                return (processName(newNamesFile));
            }
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        string outCountFile = getOutputFileName("count", variables);
        variables["[extension]"] = util.getExtension(fastafile);
        string outFastaFile = getOutputFileName("fasta", variables);
        
        CountTable ct; CountTable newCt;
        if (countfile != "")  {
            ct.readTable(countfile, true, false);
            if (countfile == outCountFile){
                //prepare filenames and open files
                map<string, string> mvariables;
                mvariables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
                mvariables["[tag]"] = "unique";
                outCountFile = getOutputFileName("count", mvariables);
            }
            newCt.copy(&ct);
        }
        
        if (m->getControl_pressed()) { return 0; }
        
        ifstream in;  util.openInputFile(fastafile, in);
        ofstream outFasta; util.openOutputFile(outFastaFile, outFasta);
        outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
        
        map<string, string> sequenceStrings; //sequenceString -> list of names.  "atgc...." -> seq1,seq2,seq3.
        map<string, string>::iterator itStrings;
        set<string> nameInFastaFile; //for sanity checking
        set<string>::iterator itname;
        vector<string> nameFileOrder;
        
        int count = 0;
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); gobble(in);
            
            if (seq.getName() != "") { //not end of file
                
                //sanity checks
                itname = nameInFastaFile.find(seq.getName());
                if (itname == nameInFastaFile.end()) { nameInFastaFile.insert(seq.getName());  }
                else { m->mothurOut("[ERROR]: You already have a sequence named " + seq.getName() + " in your fasta file, sequence names must be unique, please correct.\n");  }
                
                itStrings = sequenceStrings.find(seq.getAligned());
                string seqName = seq.getName();
                int numCurrentReps = 0;
                
                if (itStrings == sequenceStrings.end()) { //this is a new unique sequence
                    
                    seq.printSequence(outFasta); //output to unique fasta file
                    sequenceStrings[seq.getAligned()] = seqName;    nameFileOrder.push_back(seq.getAligned());
                    
                    if (countfile != "") { numCurrentReps = ct.getNumSeqs(seqName); } //checks to make sure seq is in table
                    else { newCt.push_back(seqName); }
                    
                }else { //this is a dup
                    
                    if (countfile != "") {
                        
                        numCurrentReps = newCt.getNumSeqs(seq.getName()); //checks to make sure seq is in table
                        
                        if (numCurrentReps != 0) { //its in the table
                            newCt.mergeCounts(itStrings->second, seq.getName()); //merges counts and saves in uniques name
                        }
                    }else {
                        numCurrentReps = newCt.getNumSeqs(itStrings->second);
                        newCt.setNumSeqs(itStrings->second, numCurrentReps+1);
                    }
                }
                count++;
            }
            
            
            if(count % 1000 == 0)    { m->mothurOutJustToScreen(toString(count) + "\t" + toString(sequenceStrings.size()) + "\n");    }
        }
        
        if(count % 1000 != 0)    { m->mothurOut(toString(count) + "\t" + toString(sequenceStrings.size())); m->mothurOutEndLine();    }
        
        in.close(); outFasta.close();
        if (m->getControl_pressed()) {  util.mothurRemove(outFastaFile); }
        
        //print new names file
        ofstream outCount; util.openOutputFile(outCountFile, outCount); outputTypes["count"].push_back(outCountFile); outputNames.push_back(outCountFile);
        
        newCt.printCompressedHeaders(outCount);
        
        for (int i = 0; i < nameFileOrder.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            itStrings = sequenceStrings.find(nameFileOrder[i]);
            
            if (itStrings != sequenceStrings.end()) {
                newCt.printCompressedSeq(outCount, itStrings->second);
            }else{ m->mothurOut("[ERROR]: mismatch in namefile print.\n");  m->setControl_pressed(true); }
        }
        outCount.close();
        
        if (m->getControl_pressed()) {  util.mothurRemove(outFastaFile);  util.mothurRemove(outCountFile);  outputTypes.clear(); outputNames.clear(); }
        
        return outFastaFile;
        
    }
    catch(exception& e) {
        m->errorOut(e, "UniqueSeqsCommand", "processCount");
        exit(1);
    }
}
/**************************************************************************************/
string UniqueSeqsCommand::processName(string namefile) { //namefile can be blank, indicating no namefile provided
    try {
        
        if (format == "count") { //if user want to convert the count file to a names file, do that first
            if (namefile == "") {  return (processCount("")); }
            else {
                //convert name to count
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(namefile));
                string newCountFile = getOutputFileName("count", variables);
                
                CountTable ct; ct.createTable(namefile, "", nullVector);
                ct.printCompressedTable(newCountFile);
                
                return (processCount(newCountFile));
            }
        }
        
        //assumes output is name and namefile is given
        //prepare filenames
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        string outNameFile = getOutputFileName("name", variables);
        variables["[extension]"] = util.getExtension(fastafile);
        string outFastaFile = getOutputFileName("fasta", variables);
        
        //check for name overwrite
        if (namefile == outNameFile){
            map<string, string> mvariables;
            mvariables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
            mvariables["[tag]"] = "unique";
            outNameFile = getOutputFileName("name", mvariables);
        }
        
        map<string, string> nameMap; map<string, string>::iterator itNames;
        if (namefile != "") {  util.readNames(namefile, nameMap); } //add existing duplicates
        
        //bail if error reading namefile
        if (m->getControl_pressed()) { return ""; }
        
        //open files
        ifstream in;  util.openInputFile(fastafile, in);
        ofstream outFasta; util.openOutputFile(outFastaFile, outFasta);
        outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
        
        map<string, string> sequenceStrings; //sequenceString -> list of names.  "atgc...." -> seq1,seq2,seq3.
        map<string, string>::iterator itStrings;
        set<string> nameInFastaFile; //for sanity checking
        set<string>::iterator itname;
        vector<string> nameFileOrder;
        
        int count = 0;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); gobble(in);
            
            if (seq.getName() != "") { //not end of file
                
                //sanity checks
                itname = nameInFastaFile.find(seq.getName());
                if (itname == nameInFastaFile.end()) { nameInFastaFile.insert(seq.getName());  }
                else { m->mothurOut("[ERROR]: You already have a sequence named " + seq.getName() + " in your fasta file, sequence names must be unique, please correct.\n");  }
                
                string key = seq.getName();
                
                if (namefile != "") {
                    itNames = nameMap.find(seq.getName());
                    
                    if (itNames == nameMap.end()) { //namefile and fastafile do not match
                        m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file, and not in your namefile, please correct.\n");
                    }else { key = itNames->second;  }
                }

                itStrings = sequenceStrings.find(seq.getAligned());
                
                if (itStrings == sequenceStrings.end()) { //this is a new unique sequence
                    
                    seq.printSequence(outFasta); //output to unique fasta file
                    
                    //add new unique sequence to seqStrings
                    sequenceStrings[seq.getAligned()] = key;    nameFileOrder.push_back(seq.getAligned());
                    
                }else { //this is a dup
                     sequenceStrings[seq.getAligned()] += "," + key;
                }
                count++;
            }
            
            
            if(count % 1000 == 0)    { m->mothurOutJustToScreen(toString(count) + "\t" + toString(sequenceStrings.size()) + "\n");    }
        }
        
        if(count % 1000 != 0)    { m->mothurOut(toString(count) + "\t" + toString(sequenceStrings.size())); m->mothurOutEndLine();    }
        
        in.close(); outFasta.close();
        if (m->getControl_pressed()) {  return outFastaFile;  }
        
        ofstream outNames; util.openOutputFile(outNameFile, outNames);
        outputNames.push_back(outNameFile); outputTypes["name"].push_back(outNameFile);
         
        for (int i = 0; i < nameFileOrder.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            itStrings = sequenceStrings.find(nameFileOrder[i]);
            
            if (itStrings != sequenceStrings.end()) {
                //get rep name
                int pos = (itStrings->second).find_first_of(',');
                
                if (pos == string::npos) { // only reps itself
                    outNames << itStrings->second << '\t' << itStrings->second << endl;
                }else {
                    outNames << (itStrings->second).substr(0, pos) << '\t' << itStrings->second << endl;
                }
                
            }else{ m->mothurOut("[ERROR]: mismatch in namefile print.\n");  m->setControl_pressed(true); }
        }
        outNames.close();
        
        if (m->getControl_pressed()) {  util.mothurRemove(outFastaFile);  util.mothurRemove(outNameFile);  outputTypes.clear(); outputNames.clear(); }
        
        return outFastaFile;
    }
    catch(exception& e) {
        m->errorOut(e, "UniqueSeqsCommand", "processName");
        exit(1);
    }
}
/**************************************************************************************/
string UniqueSeqsCommand::createNewNameFile(string countfile, map<string, int> nameMap) { //namefile can be blank, indicating no namefile provided
    try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(countfile));
        string outNameFile = getOutputFileName("name", variables);
        
        ofstream out; util.openOutputFile(outNameFile, out);
         
        for (map<string, int>::iterator it = nameMap.begin(); it != nameMap.end(); it++) {
            string seqName = it->first;
            string expandedName = seqName;
            int numSeqs = it->second;
            
            for (int i = 1; i < numSeqs; i++) {  expandedName += "," + seqName + "_" + toString(i);  }
            
            out << seqName << '\t' << expandedName << endl;
        }

        out.close();
        
        return outNameFile;
    }
    catch(exception& e) {
        m->errorOut(e, "UniqueSeqsCommand", "createNewNameFile");
        exit(1);
    }
}
/**************************************************************************************/
