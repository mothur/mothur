/*
 *  removeseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "removeseqscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"
#include "counttable.h"
#include "fastqread.h"
#include "inputdata.h"
#include "contigsreport.hpp"
#include "alignreport.hpp"

//**********************************************************************************************************************
vector<string> RemoveSeqsCommand::setParameters(){	
	try {
        CommandParameter pfastq("fastq", "InputTypes", "", "", "none", "FNGLT", "none","fastq",false,false,true); parameters.push_back(pfastq);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "FNGLT", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "FNGLT", "none","group",false,false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none","list",false,false,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none","taxonomy",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "none", "FNGLT", "none","alignreport",false,false); parameters.push_back(palignreport);
        CommandParameter pcontigsreport("contigsreport", "InputTypes", "", "", "FNGLT", "FNGLT", "none","contigsreport",false,false); parameters.push_back(pcontigsreport);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "FNGLT", "none","qfile",false,false); parameters.push_back(pqfile);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(paccnos);
		CommandParameter pdups("dups", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pdups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["fastq"] = tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["alignreport"] = tempOutNames;
        outputTypes["contigsreport"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.seqs command reads an .accnos file and at least one of the following file types: fasta, name, group, count, list, taxonomy, quality, fastq, contigsreport or alignreport file.\n";
		helpString += "It outputs a file containing the sequences NOT in the .accnos file.\n";
		helpString += "The remove.seqs command parameters are accnos, fasta, name, group, count, list, taxonomy, qfile, alignreport, contigsreport, fastq and dups.  You must provide accnos and at least one of the file parameters.\n";
        helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
		helpString += "The dups parameter allows you to remove the entire line from a name file if you remove any name from the line. default=true. \n";
		helpString += "The remove.seqs command should be in the following format: remove.seqs(accnos=yourAccnos, fasta=yourFasta).\n";
		helpString += "Example remove.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],pick,[extension]";    }
        else if (type == "fastq")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],pick,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],pick,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "list")        {   pattern = "[filename],[distance],pick,[extension]";    }
        else if (type == "qfile")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "alignreport")      {   pattern = "[filename],pick.[extension]";    }
        else if (type == "contigsreport")    {   pattern = "[filename],pick.[extension]";  }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

RemoveSeqsCommand::RemoveSeqsCommand(string accnos, string dupsFile, string dupsFileType, string output) {
    try {
        names = util.readAccnos(accnos);
    
        outputdir = output; dups = true; abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
        if (dupsFile != "") {
            if (dupsFileType == "count")    { readCount(dupsFile); }
            else                            {  readName(dupsFile); }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveSeqsCommand", "RemoveSeqsCommand - mothurRun");
        exit(1);
    }
}
//**********************************************************************************************************************

RemoveSeqsCommand::RemoveSeqsCommand(unordered_set<string> n, string dupsFile, string dupsFileType, string output) {
    try {
        names = n;
    
        outputdir = output; dups = true; abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
        if (dupsFile != "") {
            if (dupsFileType == "count")    { readCount(dupsFile); }
            else                            {  readName(dupsFile); }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveSeqsCommand", "RemoveSeqsCommand - mothurRun");
        exit(1);
    }
}
//**********************************************************************************************************************
RemoveSeqsCommand::RemoveSeqsCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  
				accnosfile = current->getAccnosFile();
				if (accnosfile != "") {  m->mothurOut("Using " + accnosfile + " as input file for the accnos parameter.\n");  }
				else { 
					m->mothurOut("[ERROR]: You have no valid accnos file and accnos is required.\n"); abort = true;
				}  
			}else { current->setAccnosFile(accnosfile); }	
			
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
                
            string usedDups = "true";
			string temp = validParameter.valid(parameters, "dups");
			if (temp == "not found") { 
				if (namefiles.size() != 0)  {  temp = "true";					}
				else				        {  temp = "false"; usedDups = "";	}
			}
			dups = util.isTrue(temp);
            
            format = validParameter.valid(parameters, "format");		if (format == "not found"){	format = "illumina1.8+";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  {
                m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting.\n" );  abort=true;
            }
		}
        
        if (!abort) { names = util.readAccnos(accnosfile); }
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "RemoveSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveSeqsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (m->getControl_pressed()) { return 0; }
		
		//read through the correct file and output lines you want to keep
        if (namefiles.size() != 0)          { for (int i = 0; i < namefiles.size(); i++)    { readName(namefiles[i]);       }  }
        if (fastafiles.size() != 0)         { for (int i = 0; i < fastafiles.size(); i++)   { readFasta(fastafiles[i]);     }  }
        if (qualityfiles.size() != 0)       { for (int i = 0; i < qualityfiles.size(); i++) { readQual(qualityfiles[i]);    }  }
        if (groupfiles.size() != 0)         { for (int i = 0; i < groupfiles.size(); i++)   { readGroup(groupfiles[i]);     }  }
        if (taxfiles.size() != 0)           { for (int i = 0; i < taxfiles.size(); i++)     { readTax(taxfiles[i]);         }  }
        if (listfiles.size() != 0)          { for (int i = 0; i < listfiles.size(); i++)     { readList(listfiles[i]);         }  }
        if (alignfiles.size() != 0)         { for (int i = 0; i < alignfiles.size(); i++)   { readAlign(alignfiles[i]);     }  }
        if (countfiles.size() != 0)         { for (int i = 0; i < countfiles.size(); i++)   { readCount(countfiles[i]);     }  }
        if (fastqfiles.size() != 0)         { for (int i = 0; i < fastqfiles.size(); i++)   { readFastq(fastqfiles[i]);     }  }
        if (contigsreportfiles.size() != 0) { for (int i = 0; i < contigsreportfiles.size(); i++)   { readContigs(contigsreportfiles[i]);    }  }


		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
	
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names:\n");
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]+"\n"); 	}
			m->mothurOutEndLine();
			
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
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
			}
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
			}
			
			itTypes = outputTypes.find("qfile");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
			}	
            
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
			}
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveSeqsCommand::readFasta(string fastafile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string outputFileName = getOutputFileName("fasta", variables);
		
		ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(fastafile, in);
		
        string name; bool wroteSomething = false; int removedCount = 0;
        set<string> uniqueNames;
        
		while(!in.eof()){
            if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName); return; }
			
			Sequence currSeq(in);
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(currSeq.getName());
                if (it != uniqueMap.end()) { currSeq.setName(it->second); }
            }

			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) == 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        uniqueNames.insert(name);
                        wroteSomething = true;
					
                        currSeq.printSequence(out);
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fasta file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }
				}else {  removedCount++;  }
			}
			util.gobble(in);
		}
		in.close();	out.close();
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + fastafile + " contains only sequences from the .accnos file.\n");  }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + fastafile + ".\n");
		
		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readGZFastq(string fastqfile){
    try {
        
#ifdef USE_BOOST
        
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastqfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastqfile));
        variables["[extension]"] = ".fastq" + util.getExtension(fastqfile);
        string outputFileName = getOutputFileName("fastq", variables);
        
        ifstream in; boost::iostreams::filtering_istream inBoost;
        util.openInputFileBinary(fastqfile, in, inBoost);
        
        ofstream file; ostream* out; boost::iostreams::filtering_streambuf<boost::iostreams::output> outBoost;
        util.openOutputFileBinary(outputFileName, file, out, outBoost);
        
        bool wroteSomething = false; int removedCount = 0; set<string> uniqueNames;
        
        while(!inBoost.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            //read sequence name
            bool ignore;
            FastqRead fread(inBoost, ignore, format);  util.gobble(inBoost);
            
            if (!ignore) {
                string name = fread.getName();
                
                if (names.count(name) == 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        wroteSomething = true;
                        fread.printFastq(*out);
                        uniqueNames.insert(name);
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fastq file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }
                }else { removedCount++; }
            }
            util.gobble(inBoost);
        }
        in.close(); inBoost.pop();
        boost::iostreams::close(outBoost);
        file.close(); delete out;
        
        if (m->getControl_pressed()) { util.mothurRemove(outputFileName); return; }
        
        if (wroteSomething == false) { m->mothurOut("[WARNING]: " + fastqfile + " contains only sequences from the .accnos file.\n");   }
        outputNames.push_back(outputFileName);  outputTypes["fastq"].push_back(outputFileName);
        
        m->mothurOut("Removed " + toString(removedCount) + " sequences from " + fastqfile + ".\n");
#else
        m->mothurOut("[ERROR]: mothur requires the boost libraries to read and write compressed files. Please decompress your files and rerun.\n");

#endif
        
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveSeqsCommand", "readFastq");
        exit(1);
    }
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readFastq(string fastqfile){
	try {
        
        bool gz = util.isGZ(fastqfile)[1];
        
        if (gz) { readGZFastq(fastqfile); return; }
        
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(fastqfile);  }
		map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastqfile));
        variables["[extension]"] = util.getExtension(fastqfile);
		string outputFileName = getOutputFileName("fastq", variables);
		
        ofstream out; util.openOutputFile(outputFileName, out);
        ifstream in; util.openInputFile(fastqfile, in);
        
        bool wroteSomething = false; int removedCount = 0;
		set<string> uniqueNames;
        
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName); return; }
			
            //read sequence name
            bool ignore; FastqRead fread(in, ignore, format); util.gobble(in);
            
            if (!ignore) {
                string name = fread.getName();
                
                if (names.count(name) == 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        wroteSomething = true;
                        fread.printFastq(out);
                        uniqueNames.insert(name);
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fastq file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }

                }else { removedCount++; }
            }
            
			util.gobble(in);
		}
		in.close(); out.close();
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + fastqfile + " contains only sequences from the .accnos file.\n");   }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + fastqfile + ".\n");

		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readFastq");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readQual(string qualfile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(qualfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(qualfile));
        variables["[extension]"] = util.getExtension(qualfile);
		string outputFileName = getOutputFileName("qfile", variables);
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(qualfile, in);
		
        string name; bool wroteSomething = false; int removedCount = 0;
		set<string> uniqueNames;
        
		while(!in.eof()){	
			
            QualityScores qual(in); util.gobble(in);
			
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(qual.getName());
                if (it != uniqueMap.end()) { qual.setName(it->second); }
            }
            
            string name = qual.getName();
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
				
                    qual.printQScores(out);
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your qfile more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {  removedCount++;  }
			
			util.gobble(in);
		}
		in.close(); out.close();
		
		if (wroteSomething == false) { m->mothurOut("[WARNING]: " + qualfile + " contains only sequences from the .accnos file.\n");   }
		outputNames.push_back(outputFileName);  outputTypes["qfile"].push_back(outputFileName); 
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + qualfile + ".\n");
		
		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readQual");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readCount(string countfile){
	try {
        
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(countfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[extension]"] = util.getExtension(countfile);
		string outputFileName = getOutputFileName("count", variables);
		
        CountTable ct; ct.readTable(countfile, true, false); int originalCount = ct.getNumSeqs();
        
        for (auto it = names.begin(); it != names.end(); it++) {
            ct.zeroOutSeq(*it);
            if (m->getControl_pressed()) {  return; }
        }
        
        ct.printTable(outputFileName);
        
        int removedCount = originalCount - ct.getNumSeqs();
        
		outputTypes["count"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + countfile + ".\n");
        
        if (ct.getNumSeqs() == 0) {  m->mothurOut("[WARNING]: " + countfile + " contains only sequences from the .accnos file.\n");  }
        
		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readCount");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readList(string listfile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(listfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        InputData input(listfile, "list", nullVector);
        ListVector* list = input.getListVector();
		
		bool wroteSomething = false; int removedCount = 0;
        
		while(list != nullptr) {
			
			removedCount = 0;
            set<string> uniqueNames;
			
            //make a new list vector
			ListVector newList;
			newList.setLabel(list->getLabel());
            
			variables["[distance]"] = list->getLabel();
            string outputFileName = getOutputFileName("list", variables);
			
			ofstream out;
			util.openOutputFile(outputFileName, out);
			outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            vector<string> binLabels = list->getLabels();
            vector<string> newBinLabels;
            
            if (m->getControl_pressed()) { out.close();  return; }

			//for each bin
			for (int i = 0; i < list->getNumBins(); i++) {
				if (m->getControl_pressed()) {  out.close();  util.mothurRemove(outputFileName);  return; }
			
				//parse out names that are in accnos file
				string bin = list->get(i);
                vector<string> bnames;
                util.splitAtComma(bin, bnames);
				
				string newNames = "";
                for (int j = 0; j < bnames.size(); j++) {
					string name = bnames[j];
                    //if that name is in the .accnos file, add it
					if (names.count(name) == 0) {
                        if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                            uniqueNames.insert(name);
                            newNames += name + ",";
                        }else {
                            m->mothurOut("[WARNING]: " + name + " is in your list file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                        }
                    }
					else {  removedCount++;  }
                }

				//if there are names in this bin add to new list
				if (newNames != "") {  
					newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
					newList.push_back(newNames);
                    newBinLabels.push_back(binLabels[i]);
				}
			}
				
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.setLabels(newBinLabels);
				newList.print(out, false);
			}

            out.close();
            
            delete list;
            list = input.getListVector();
		}
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + listfile + " contains only sequences from the .accnos file.\n");  }
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + listfile + ".\n");
		
		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readName(string namefile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(namefile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
        variables["[extension]"] = util.getExtension(namefile);
		string outputFileName = getOutputFileName("name", variables);
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(namefile, in);
        
		string name, firstCol, secondCol; bool wroteSomething = false; int removedCount = 0;
        set<string> uniqueNames;
        
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
			in >> firstCol;		util.gobble(in);		
			in >> secondCol;	util.gobble(in);
			
			vector<string> parsedNames;
			util.splitAtComma(secondCol, parsedNames);
			
            vector<string> validSecond;  validSecond.clear(); vector<string> parsedNames2;
            bool parsedError = false;
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
                    if (uniqueNames.count(parsedNames[i]) == 0) { //this name hasn't been seen yet
                        uniqueNames.insert(parsedNames[i]);
                        validSecond.push_back(parsedNames[i]);
                        parsedNames2.push_back(parsedNames[i]);
                    }else {
                        m->mothurOut("[WARNING]: " + parsedNames[i] + " is in your name file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                        parsedError = true;
                    }
				}
			}
            
            if (parsedError) {  parsedNames = parsedNames2; }
			
			if ((dups) && (validSecond.size() != parsedNames.size())) {  //if dups is true and we want to get rid of anyone, get rid of everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]);  }
				removedCount += parsedNames.size();
			}else {
                if (validSecond.size() != 0) {
                    removedCount += parsedNames.size()-validSecond.size();
                    //if the name in the first column is in the set then print it and any other names in second column also in set
                    if (names.count(firstCol) == 0) {
                        
                        wroteSomething = true;
                        
                        out << firstCol << '\t';
                        
                        //you know you have at least one valid second since first column is valid
                        for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
                        out << validSecond[validSecond.size()-1] << endl;
                        
                        //make first name in set you come to first column and then add the remaining names to second column
                    }else {
                        
                        //you want part of this row
                        if (validSecond.size() != 0) {
                            
                            wroteSomething = true;
                            
                            out << validSecond[0] << '\t';
                            //we are changing the unique name in the fasta file
                            uniqueMap[firstCol] = validSecond[0];
                            
                            //you know you have at least one valid second since first column is valid
                            for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
                            out << validSecond[validSecond.size()-1] << endl;
                        }
                    }
                }
			}
		}
		in.close(); out.close();
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + namefile + " contains only sequences from the .accnos file.\n"); }
		outputTypes["name"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + namefile + ".\n");
		
		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveSeqsCommand::readGroup(string groupfile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(groupfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[extension]"] = util.getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);	
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(groupfile, in);
        
		string name, group; bool wroteSomething = false; int removedCount = 0;
        set<string> uniqueNames;
        
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
			in >> name;			util.gobble(in);		//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                    out << name << '\t' << group << endl;
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your group file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {  removedCount++;  }
					
			util.gobble(in);
		}
		in.close(); out.close();
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + groupfile + " contains only sequences from the .accnos file.\n");   }
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + groupfile + ".\n");

		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readTax(string taxfile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(taxfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
        variables["[extension]"] = util.getExtension(taxfile);
		string outputFileName = getOutputFileName("taxonomy", variables);
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(taxfile, in);
		
        string name, tax; bool wroteSomething = false; int removedCount = 0;
        set<string> uniqueNames;
        
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
            in >> name; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
            
			//if this name is in the accnos file
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
            
                    out << name << '\t' << tax << endl;
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your taxonomy file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {  removedCount++;  }
		}
		in.close(); out.close();
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + taxfile + " contains only sequences from the .accnos file.\n");   }
		outputTypes["taxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + taxfile + ".\n");
		
		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
void RemoveSeqsCommand::readAlign(string alignfile){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(alignfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(alignfile));
        variables["[extension]"] = util.getExtension(alignfile);
		string outputFileName = getOutputFileName("alignreport", variables);
		
		ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(alignfile, in);
		
		bool wroteSomething = false; int removedCount = 0;
        set<string> uniqueNames;
		
        AlignReport report; report.readHeaders(in); util.gobble(in); report.printHeaders(out);
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
			report.read(in); util.gobble(in);
            string name = report.getQueryName();
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                    
                    report.print(out);
                    
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your alignreport file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else { removedCount++;   }
		}
		in.close(); out.close();
		
		if (wroteSomething == false) {  m->mothurOut("[WARNING]: " + alignfile + " contains only sequences from the .accnos file.\n");  ofstream out1; util.openOutputFile(outputFileName, out1); out1.close(); }
		outputTypes["alignreport"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from " + alignfile + ".\n");

		return;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************
//contigsreport file has a column header line then all other lines contain 8 columns.  we just want the first column since that contains the name
void RemoveSeqsCommand::readContigs(string contigsreportfile){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(contigsreportfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(contigsreportfile));
        variables["[extension]"] = util.getExtension(contigsreportfile);
        string outputFileName = getOutputFileName("contigsreport", variables);
        
        ofstream out; util.openOutputFile(outputFileName, out);
        ifstream in; util.openInputFile(contigsreportfile, in);
        
        bool wroteSomething = false; int removedCount = 0;
        set<string> uniqueNames;

        ContigsReport report; report.readHeaders(in); util.gobble(in); report.printHeaders(out);
        
        while(!in.eof()){
        
            if (m->getControl_pressed()) { break; }

            report.read(in); util.gobble(in);
            string name = report.getName();
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
            
            if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                
                    report.print(out);
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your contigsreport file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
            }else { removedCount++; }
        }
        in.close(); out.close();
        
        if (wroteSomething == false) { m->mothurOut("[WARNING]: " + contigsreportfile + " only contains sequences from the .accnos file, everything removed.\n");  ofstream out1; util.openOutputFile(outputFileName, out1); out1.close(); } //reopening file clears header line
        outputNames.push_back(outputFileName);  outputTypes["contigsreport"].push_back(outputFileName);
        
        m->mothurOut("Removed " + toString(removedCount) + " sequences from " + contigsreportfile + ".\n");
        
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveSeqsCommand", "readContigs");
        exit(1);
    }
}
//**********************************************************************************************************************


