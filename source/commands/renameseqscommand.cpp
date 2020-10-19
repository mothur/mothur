//
//  renameseqscommand.cpp
//  Mothur
//
//  Created by SarahsWork on 5/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "renameseqscommand.h"
#include "sequence.hpp"
#include "groupmap.h"
#include "counttable.h"
#include "qualityscores.h"
#include "contigsreport.hpp"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> RenameSeqsCommand::setParameters(){
	try {
        CommandParameter pfile("file", "InputTypes", "", "", "fileFasta-file", "fileFasta", "none","fasta",false,false,true); parameters.push_back(pfile);
        CommandParameter pmap("map", "InputTypes", "", "", "none", "none", "none","fasta",false,false,true); parameters.push_back(pmap);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "fileFasta-file", "fileFasta", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter plist("list", "InputTypes", "", "", "fileFasta-file", "fileFasta", "none","fasta",false,false,true); parameters.push_back(plist);
        CommandParameter pqfile("qfile", "InputTypes", "", "", "file", "none", "none","qfile",false,false,true); parameters.push_back(pqfile);
        CommandParameter pcontigsreport("contigsreport", "InputTypes", "", "", "file", "none", "none","contigsreport",false,false,true); parameters.push_back(pcontigsreport);
        CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "file", "none","taxonomy",false,false,true); parameters.push_back(ptaxonomy);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount-file", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup-file", "none", "none","count",false,false,true); parameters.push_back(pcount);
        CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup-file", "none", "none","group",false,false,true); parameters.push_back(pgroup);
        CommandParameter pdelim("delim", "String", "", "_", "", "", "","",false,false); parameters.push_back(pdelim);
        CommandParameter pplacement("placement", "Multiple", "front-back", "back", "", "", "","",false,false); parameters.push_back(pplacement);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        outputTypes["map"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;
        outputTypes["contigsreport"] = tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RenameSeqsCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The rename.seqs command renames sequences in the input files. By default, mothur will generate new names based on your inputs. Alternatively, you can provide a map file.\n";
        helpString += "The rename.seqs command parameters are " + getCommandParameters() + ".\n";
        helpString += "The list parameter allows you to provide an associated list file.\n";
        helpString += "The fasta parameter allows you to provide an associated fasta file.\n";
        helpString += "The qfile parameter allows you to provide an associated quality file.\n";
        helpString += "The taxonomy parameter allows you to provide an associated taxonomy file.\n";
        helpString += "The contigsreport allows you to provide an associated contigsreport file.\n";
        helpString += "The file option allows you to provide a 2 or 3 column file. The first column contains the file type: fasta or qfile. The second column is the filename, and the optional third column can be a group name. If there is a third column, all sequences in the file will be assigned to that group.  This can be helpful when renaming data separated into samples. \n";
        helpString += "The placement parameter allows you to indicate whether you would like the group name appended to the front or back of the sequence number.  Options are front or back. Default=back.\n";
        helpString += "The delim parameter allow you to enter the character or characters you would like to separate the sequence number from the group name. Default='_'.\n";
        helpString += "The rename.seqs command should be in the following format: \n";
        helpString += "The rename.seqs command should be in the following format: \n";
		helpString += "rename.seqs(fasta=yourFastaFile, group=yourGroupFile) \n";
		helpString += "Example rename.seqs(fasta=abrecovery.unique.fasta, group=abrecovery.group).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RenameSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")                    {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "name")                {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "group")               {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "count")               {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "taxonomy")            {   pattern = "[filename],renamed,[extension]";    }
        else if (type == "qfile")               {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "contigsreport")       {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "list")                {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "map")                 {  pattern = "[filename],renamed_map"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
/**************************************************************************************/
RenameSeqsCommand::RenameSeqsCommand(string option)  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta");
			if (fastaFile == "not open") { abort = true; }
			else if (fastaFile == "not found") { fastaFile = ""; }
			else { current->setFastaFile(fastaFile); }
            
            fileFile = validParameter.validFile(parameters, "file");
            if (fileFile == "not open") { abort = true; }
            else if (fileFile == "not found") { fileFile = ""; }
            else { current->setFileFile(fileFile); }
			
            groupfile = validParameter.validFile(parameters, "group");
            if (groupfile == "not open") { abort = true; }
            else if (groupfile == "not found") {  groupfile = ""; }
			else { current->setGroupFile(groupfile); }

            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = ""; }
            else {  current->setCountFile(countfile);  }
            
            nameFile = validParameter.validFile(parameters, "name");
            if (nameFile == "not open") { abort = true; }
            else if (nameFile == "not found"){ nameFile =""; }
            else { current->setNameFile(nameFile); }
            
            qualfile = validParameter.validFile(parameters, "qfile");
            if (qualfile == "not open") { abort = true; }
            else if (qualfile == "not found"){ qualfile =""; }
            else { current->setQualFile(qualfile); }
            
            listfile = validParameter.validFile(parameters, "list");
            if (listfile == "not open") { abort = true; }
            else if (listfile == "not found"){ listfile =""; }
            else { current->setListFile(listfile); }
            
            mapFile = validParameter.validFile(parameters, "map");
            if (mapFile == "not open") { abort = true; }
            else if (mapFile == "not found"){ mapFile = ""; }
            
            contigsfile = validParameter.validFile(parameters, "contigsreport");
            if (contigsfile == "not open") { abort = true; }
            else if (contigsfile == "not found"){ contigsfile = ""; }
            
            taxfile = validParameter.validFile(parameters, "taxonomy");
            if (taxfile == "not open") { taxfile = ""; abort = true; }
            else if (taxfile == "not found") {  taxfile = "";  }
            else { current->setTaxonomyFile(taxfile); }
            
            if ((countfile != "") && (nameFile != "")) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name.\n"); abort = true; }
            
            if ((fileFile != "") && (fastaFile != "")) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: file or fasta.\n"); abort = true; }
            
            if ((countfile != "") && (groupfile != "")) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group.\n");  abort = true; }
            
            if ((fileFile != "") && ((listfile != "") || (nameFile != "") || (groupfile != "") || (qualfile != "") || (contigsfile != "") || (countfile != "") || (fastaFile != "")) ) {
                m->mothurOut("[ERROR]: The file option cannot be used with any other files except the map file.\n");  abort = true;
            }else if ((fileFile == "") && (listfile == "") && (nameFile == "") && (groupfile == "") && (qualfile == "") && (contigsfile == "") && (countfile == "") && (fastaFile == "")) {
                m->mothurOut("[ERROR]: No input files provided, please correct.\n");  abort = true;
            }
            
            placement = validParameter.valid(parameters, "placement");		if (placement == "not found") { placement = "back"; }
            if ((placement == "front") || (placement == "back")) { }
            else { m->mothurOut("[ERROR]: " + placement + " is not a valid placement option.  Valid placement options are front or back.\n"); abort = true; }
            
            delim = validParameter.valid(parameters, "delim");			if (delim == "not found") { delim = "_"; }
            
            if (countfile == "") {
                if (nameFile == "")  {
                    vector<string> files; files.push_back(fastaFile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "RenameSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
int RenameSeqsCommand::execute() {
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        ignoreNew = false;
        
        map<string, string> renameMap; bool printMap = false;
        if (mapFile != "") { readMapFile(renameMap); ignoreNew = true; }
        else { printMap = true; }
        
        if (fileFile != "") {  processFile(renameMap); }
        else {
            
            map<string, string> old2NewNameMap;
            if ((nameFile != "") || (countfile != "") || (groupfile != "")) {
                processNameGroupCountFiles(renameMap, old2NewNameMap);
            }else if (mapFile != "") {
                old2NewNameMap = renameMap;
            }
            renameMap.clear();
            
            if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]);  } return 0; }
            
            if (listfile != "")         {   readList(old2NewNameMap);     }
            if (fastaFile != "")         { readFasta(old2NewNameMap);     }
            if (qualfile != "")         { readQual(old2NewNameMap);       }
            if (contigsfile != "")      { readContigs(old2NewNameMap);     }
            if (taxfile != "")          { readTax(old2NewNameMap);        }
            
            if (printMap) {
                string thisOutputDir = outputdir;
                map<string, string> variables;
                string outMapFile = thisOutputDir + util.getRootName(util.getSimpleName(fastaFile));
                variables["[filename]"] = outMapFile;
                outMapFile = getOutputFileName("map", variables);
                outputNames.push_back(outMapFile); outputTypes["map"].push_back(outMapFile);
                ofstream outMap; util.openOutputFile(outMapFile, outMap);
            
                //print map
                for(map<string, string>::iterator it = old2NewNameMap.begin(); it != old2NewNameMap.end(); it++) {
                    outMap << it->second << '\t' << it->first << endl;
                }
                outMap.close();
            }
        }
        
        if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]);  } return 0; }

        m->mothurOut("\nOutput File Names:\n ");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
        
        //set fasta file as new current fastafile
        string currentName = "";
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
        }
        
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }
        
        itTypes = outputTypes.find("name");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
        }
        
        itTypes = outputTypes.find("group");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
        }
        
        itTypes = outputTypes.find("count");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
        }
        
        itTypes = outputTypes.find("qfile");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
        }
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void RenameSeqsCommand::processNameGroupCountFiles(map<string, string>& oldMap, map<string, string>& old2NewNameMap){
    try {
        bool oldMapEmpty = true;
        if (oldMap.size() != 0) { oldMapEmpty = false; }
        
        GroupMap* groupMap = NULL;
        CountTable* countTable = NULL;
        
        bool hasGroups = false;
        vector<string> Groups;
        if (groupfile != "") {
            groupMap = new GroupMap(groupfile);
            int groupError = groupMap->readMap();
            if (groupError == 1) { delete groupMap; return; }
            Groups = groupMap->getNamesOfGroups();
            hasGroups = true;
        }else if (countfile != "") {
            countTable = new CountTable();
            countTable->readTable(countfile, true, false);
            hasGroups = countTable->hasGroupInfo();
            if (hasGroups) {     Groups = countTable->getNamesOfGroups(); Groups.push_back("Multi"); }
        }
        
        //set up for reads
        map<string, int> counts; for (int i = 0; i < Groups.size(); i++) {  counts[Groups[i]] = 1; }
        string thisOutputDir = outputdir;
        
        if (nameFile != "") {
            map<string, vector<string> > nameMap;
            
            thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir += util.hasPath(nameFile);  }
            string outNameFile = thisOutputDir + util.getRootName(util.getSimpleName(nameFile));
            map<string, string> variables;
            variables["[filename]"] = outNameFile;
            variables["[extension]"] = util.getExtension(nameFile);
            outNameFile = getOutputFileName("name", variables);
            outputNames.push_back(outNameFile); outputTypes["name"].push_back(outNameFile);
            
            ofstream outName; util.openOutputFile(outNameFile, outName);
            util.readNames(nameFile, nameMap);
            
            for (map<string, vector<string> >::iterator itNames = nameMap.begin(); itNames != nameMap.end(); itNames++) {
                vector<string> dups = itNames->second;
                
                if (m->getControl_pressed()) { break; }
                
                if (oldMapEmpty) {
                    for (int i = 0; i < dups.size(); i++) {
                        string group = "";
                    
                        if (groupfile != "") {
                            group = groupMap->getGroup(dups[i]);
                        }else if (countfile != "") {
                            if (hasGroups) {
                                vector<string> thisReadsGroups = countTable->getGroups(dups[i]);
                                if (thisReadsGroups.size() == 0)        {   group = "not found";            }
                                else if (thisReadsGroups.size() == 1)   {   group = thisReadsGroups[0];     }
                                else                                    {   group = "Multi";                }
                            }
                        }
                    
                        if (group == "not found") {  m->mothurOut("[ERROR]: " + dups[i] + " is not in your file, please correct.\n"); m->setControl_pressed(true); }
                        else {
                            
                            string newName = toString(counts[group]); counts[group]++;
                            if ((placement == "back") && (group != "")) { newName += delim + group; }
                            else if (group != "") { newName = group + delim + newName; }
                            
                            if (i == 0)  { outName << newName << '\t' << newName;  }
                            else { outName << "," << newName;  }
                            
                            oldMap[newName] = dups[i];
                            old2NewNameMap[dups[i]] = newName;
                        }
                    }
                    outName << endl;
                }else {
                    for (int i = 0; i < dups.size(); i++) {
                        //get new name
                        string newName = "";
                        
                        map<string, string>::iterator itMap = oldMap.find(dups[i]);
                        if (itMap == oldMap.end()) { m->mothurOut("[ERROR]: " + dups[i] + " is not in your map file, please correct.\n"); m->setControl_pressed(true);}
                        else { newName = itMap->second; }
                        
                        if (i == 0)  { outName << newName << '\t' << newName;  }
                        else { outName << "," << newName;  }
                        
                        old2NewNameMap[dups[i]] = newName;
                    }
                    outName << endl;
                }
            }
            outName.close();
        }
        
        if (groupfile != "") {
            thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir += util.hasPath(groupfile);  }
            string outGroupFile = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
            map<string, string> variables;
            variables["[filename]"] = outGroupFile;
            variables["[extension]"] = util.getExtension(groupfile);
            outGroupFile = getOutputFileName("group", variables);
            outputNames.push_back(outGroupFile); outputTypes["group"].push_back(outGroupFile);
            
            ofstream outGroup; util.openOutputFile(outGroupFile, outGroup);
            
            vector<string> namesOfSeqs = groupMap->getNamesSeqs();
            
            for (int i = 0; i < namesOfSeqs.size(); i++) {
                
                string group = groupMap->getGroup(namesOfSeqs[i]);
                string newName = "";
                
                if (group == "not found") {  m->mothurOut("[ERROR]: " + namesOfSeqs[i] + " is not in your file, please correct.\n"); m->setControl_pressed(true); }
                
                if (m->getControl_pressed()) { break; }
                
                if (oldMapEmpty) {
                    map<string, string>::iterator itMap = old2NewNameMap.find(namesOfSeqs[i]);
                    if (itMap == old2NewNameMap.end()) {
                        newName = toString(counts[group]); counts[group]++;
                        if ((placement == "back") && (group != "")) { newName += delim + group; }
                        else if (group != "") { newName = group + delim + newName; }
                    }else {
                        newName = itMap->second; //newName = name from namefile
                    }
                    oldMap[newName] = namesOfSeqs[i];
                }else {
                    map<string, string>::iterator itMap = oldMap.find(namesOfSeqs[i]);
                    if (itMap == oldMap.end()) { m->mothurOut("[ERROR]: " + namesOfSeqs[i] + " is not in your map file, please correct.\n"); m->setControl_pressed(true);}
                    else { newName = itMap->second; }
                }
                outGroup << newName << '\t' << group << endl;
            
                old2NewNameMap[namesOfSeqs[i]] = newName;
            }
            outGroup.close();
        }
        
        if (countfile != "") {
            thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir += util.hasPath(countfile);  }
            string outCountFile = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
            map<string, string> variables;
            variables["[filename]"] = outCountFile;
            variables["[extension]"] = util.getExtension(countfile);
            outCountFile = getOutputFileName("count", variables);
            outputNames.push_back(outCountFile); outputTypes["count"].push_back(outCountFile);
            
            vector<string> namesOfSeqs = countTable->getNamesOfSeqs();
            
            for (int i = 0; i < namesOfSeqs.size(); i++) {
                
                if (m->getControl_pressed()) { break; }
                
                string newName = "";
                if (oldMapEmpty) {
                    map<string, string>::iterator itMap = old2NewNameMap.find(namesOfSeqs[i]);
                    if (itMap == old2NewNameMap.end()) {
                        string group = "";
                        if (hasGroups) {
                            vector<string> thisReadsGroups = countTable->getGroups(namesOfSeqs[i]);
                            if (thisReadsGroups.size() == 0)        {   group = "not found";            }
                            else if (thisReadsGroups.size() == 1)   {   group = thisReadsGroups[0];     }
                            else                                    {   group = "Multi";                }
                        }
                        
                        if (group == "not found") {  m->mothurOut("[ERROR]: " + namesOfSeqs[i] + " is not in your file, please correct.\n"); m->setControl_pressed(true); }
                        
                        newName = toString(counts[group]); counts[group]++;
                        if ((placement == "back") && (group != "")) { newName += delim + group; }
                        else if (group != "") { newName = group + delim + newName; }
                    }else { newName = itMap->second; }
                    oldMap[newName] = namesOfSeqs[i];
                }else {
                    
                    map<string, string>::iterator itMap = oldMap.find(namesOfSeqs[i]);
                    if (itMap == oldMap.end()) { m->mothurOut("[ERROR]: " + namesOfSeqs[i] + " is not in your map file, please correct.\n"); m->setControl_pressed(true);}
                    else { newName = itMap->second; }
                }
                countTable->renameSeq(namesOfSeqs[i], newName);
        
                old2NewNameMap[namesOfSeqs[i]] = newName;
            }
            
            countTable->printTable(outCountFile);
        }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } }
        
        if (groupMap != NULL) { delete groupMap; }
        if (countTable != NULL) { delete countTable; }
            
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
void RenameSeqsCommand::readFasta(map<string, string>& oldMap){
    try {
        //prepare filenames and open files
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastaFile);  }
        string outFastaFile = thisOutputDir + util.getRootName(util.getSimpleName(fastaFile));
        map<string, string> variables;
        variables["[filename]"] = outFastaFile;
        variables["[extension]"] = util.getExtension(fastaFile);
        outFastaFile = getOutputFileName("fasta", variables);
        outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);

        ifstream in; util.openInputFile(fastaFile, in);
        ofstream out; util.openOutputFile(outFastaFile, out);
        
        map<string, string>::iterator it;
        int count = 0;
        while(!in.eof()){
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); util.gobble(in);
            
            it = oldMap.find(seq.getName());
            if (it == oldMap.end()) { //not in other files, create name
                if (!ignoreNew) {
                    seq.setName(toString(count));
                    oldMap[seq.getName()] = toString(count); count++;
                }
            }else {
                seq.setName(it->second);
            }
            
            seq.printSequence(out);
        }
        in.close(); out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readFasta");
        exit(1);
    }
}
//**********************************************************************************************************************
void RenameSeqsCommand::readQual(map<string, string>& oldMap){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(qualfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(qualfile));
        variables["[extension]"] = util.getExtension(qualfile);
        string outputFileName = getOutputFileName("qfile", variables);
        ofstream out;
        util.openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName); outputTypes["qfile"].push_back(outputFileName);

        ifstream in;
        util.openInputFile(qualfile, in);
        
        map<string, string>::iterator it;
        int count = 0;
        while(!in.eof()){
            if (m->getControl_pressed()) { break; }
            
            QualityScores qual(in); util.gobble(in);
            
            it = oldMap.find(qual.getName());
            if (it == oldMap.end()) {
                if (!ignoreNew) {
                    qual.setName(toString(count));
                    oldMap[qual.getName()] = toString(count); count++;
                }
            }else {
                qual.setName(it->second);
            }
            
            qual.printQScores(out);
        }
        in.close(); out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readQual");
        exit(1);
    }
}
//**********************************************************************************************************************
void RenameSeqsCommand::readTax(map<string, string>& oldMap){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(taxfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
        variables["[extension]"] = util.getExtension(taxfile);
        string outputFileName = getOutputFileName("taxonomy", variables);
        outputNames.push_back(outputFileName);  outputTypes["taxonomy"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out);
        
        ifstream in; util.openInputFile(taxfile, in);
        string name, tax; int count = 0;
        
        map<string, string>::iterator it;
        while(!in.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            in >> name; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
            
            it = oldMap.find(name);
            if (it == oldMap.end()) {
                if (!ignoreNew) {
                    name = toString(count);
                    oldMap[name] = toString(count); count++;
                }
            }else {
                name = it->second;
            }
            
            out << name << '\t' << tax << endl;
        }
        in.close(); out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readTax");
        exit(1);
    }
}
//**********************************************************************************************************************
void RenameSeqsCommand::readContigs(map<string, string>& oldMap){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(contigsfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(contigsfile));
        variables["[extension]"] = util.getExtension(contigsfile);
        string outputFileName = getOutputFileName("contigsreport", variables);
        ofstream out;
        util.openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName); outputTypes["contigsreport"].push_back(outputFileName);

        ifstream in; util.openInputFile(contigsfile, in);
        ContigsReport report;
        report.readHeaders(in); util.gobble(in);
        report.printHeaders(out);
        
        map<string, string>::iterator it;
        int count = 0;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            report.read(in); util.gobble(in);
            
            it = oldMap.find(report.getName());
            
            if (it != oldMap.end()) { report.setName(it->second); }
            else {
                if (!ignoreNew) {
                    oldMap[report.getName()] = toString(count);
                    report.setName(toString(count)); count++;
                }
            }
            report.print(out);
        }
        in.close(); out.close();

    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readContigs");
        exit(1);
    }
}
//**********************************************************************************************************************
void RenameSeqsCommand::readList(map<string, string>& oldMap){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(listfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        string outputFileName = getOutputFileName("list", variables);
        ofstream out; util.openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName); outputTypes["list"].push_back(outputFileName);
        
        InputData input(listfile, "list", nullVector);
        set<string> processedLabels;
        set<string> userLabels;
        string lastLabel = "";
        bool printHeaders = true;
        
        ListVector* list = util.getNextList(input, true, userLabels, processedLabels, lastLabel);
        
        while (list != NULL) {
            
            if (m->getControl_pressed()) { delete list; break; }
            
            list->setPrintedLabels(printHeaders);
            
            //process list
            int count = 0;
            for (int i = 0; i < list->getNumBins(); i++) {
                string bin = list->get(i);
                vector<string> names; util.splitAtComma(bin, names);
                
                for (int j = 0; j < names.size(); j++) {
                    map<string, string>::iterator it = oldMap.find(names[j]);
                    if (it == oldMap.end()) {
                        if (!ignoreNew) {
                            string newName = toString(count); count++;
                            oldMap[names[j]] = newName;
                            names[j] = newName;
                        }
                    }else {
                        names[j] = it->second;
                    }
                }
                bin = util.getStringFromVector(names, ",");
                list->set(i, bin);
            }
            
            //print list
            list->print(out); printHeaders = false;
            
            delete list;
            
            list = util.getNextList(input, true, userLabels, processedLabels, lastLabel);
        }
        out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readList");
        exit(1);
    }
}
//**********************************************************************************************************************
int RenameSeqsCommand::readMapFile(map<string, string>& readMap){
    try {
        ifstream in;
        util.openInputFile(mapFile, in);
        
        map<string, string>::iterator it;
        string oldname, newname;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            in >> oldname; util.gobble(in);
            in >> newname; util.gobble(in);
            
            it = readMap.find(oldname);
            if (it != readMap.end()) {
                m->mothurOut("[ERROR]: " + oldname + " is already in your map file. Sequence names must be unique, quitting.\n"); m->setControl_pressed(true);
            }else {
                readMap[oldname] = newname;
            }
            
        }
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readMapFile");
        exit(1);
    }
}
//**********************************************************************************************************************
int RenameSeqsCommand::processFile(map<string, string>& readMap){
    try {
        
        vector<map<string, string> > files = readFiles();
        
        for (int i = 0; i < files.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            string thisFile = "";
            string thisFileType = "";
            string group = (files[i].begin())->second;
            string temp = (files[i].begin())->first;
            
            int pos = temp.find_first_of('-');
            if (pos == string::npos) {  m->mothurOut("[Opps]: I should never get here...\n");  }
            else {
                thisFileType = temp.substr(0, pos);
                thisFile = temp.substr(pos+1);
            }
            
            string thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir += util.hasPath(thisFile);  }
            string outMapFile = thisOutputDir + util.getSimpleName(thisFile);
            map<string, string> variables;
            variables["[filename]"] = outMapFile + ".";
            outMapFile = getOutputFileName("map", variables);
            outputNames.push_back(outMapFile); outputTypes["map"].push_back(outMapFile);
            ofstream outMap; util.openOutputFile(outMapFile, outMap);
            
            //prepare filenames and open files
            string outFile = thisOutputDir + util.getRootName(util.getSimpleName(thisFile));
            variables["[filename]"] = outFile;
            variables["[extension]"] = util.getExtension(thisFile);
            outFile = getOutputFileName(thisFileType, variables);
            outputNames.push_back(outFile); outputTypes[thisFileType].push_back(outFile);
            
            ofstream out;
            util.openOutputFile(outFile, out);
            
            ifstream in;
            util.openInputFile(thisFile, in);
            
            int count = 1;
            while (!in.eof()) {
                if (m->getControl_pressed()) { break; }
                
                Sequence* seq;  QualityScores* qual;
                string name = "";
                if (thisFileType == "fasta")        {  seq = new Sequence(in);   name = seq->getName();        }
                else if (thisFileType == "qfile")   {  qual = new QualityScores(in);  name = qual->getName();  }
                else if (thisFileType == "group")   {  in >> name; util.gobble(in); in >> group;               }
                util.gobble(in);
 
                //get new name
                string newName = "";
                if (mapFile != "") {
                    map<string, string>::iterator itMap = readMap.find(name);
                    if (itMap == readMap.end()) { m->mothurOut("[ERROR]: " + name + " is not in your map file, please correct.\n"); m->setControl_pressed(true);}
                    else { newName = itMap->second; }
                }else {
                    newName = toString(count); count++;
                    if ((placement == "back") && (group != "")) { newName += delim + group; }
                    else if (group != "") { newName = group + delim + newName; }
                }
                
                
                if (thisFileType == "fasta")        {  seq->setName(newName);  seq->printSequence(out); delete seq;     }
                if (thisFileType == "qfile")        {  qual->setName(newName);  qual->printQScores(out); delete qual;   }
                if (thisFileType == "group")        {  out << newName << '\t' << group << endl;                         }
                
                outMap << newName << '\t' << name << endl;
            }
            in.close();
            outMap.close();
            out.close();
        }
        
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "processFile");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<map<string, string> > RenameSeqsCommand::readFiles(){
    try {
        
        FileFile dataFile(fileFile, "renameSeqs");
        vector< vector<string> > dataFiles = dataFile.getFiles();
        int dataFileFormat = dataFile.getFileFormat();
        map<int, string> file2Group = dataFile.getFile2Group();
        
        if (dataFileFormat == 3) { //4 column
            m->mothurOut("[ERROR]: Your file contains " + toString(4) + " columns. The rename.seqs file option allows you to provide a 2 or 3 column file. The first column contains the file type: fasta or qfile. The second column is the filename, and the optional third column can be a group name. If there is a third column, all sequences in the file will be assigned to that group.  This can be helpful when renaming data separated into samples.\n"); m->setControl_pressed(true);
        }

        vector<map<string, string> > files;
        for (int i = 0; i < dataFiles.size(); i++) {
            string group = "";
            string thisFileName; thisFileName = "";
            string thisFileName1, thisFileName2;
            thisFileName1 = dataFiles[i][0]; thisFileName2 = dataFiles[i][1];
            
            if (dataFileFormat == 1) { //2 column
                thisFileName = thisFileName1+"-"+thisFileName2;
            }else if (dataFileFormat == 2) { //3 column
                thisFileName = thisFileName1+"-"+thisFileName2;
                group = file2Group[i];
            }
            
            map<string, string> temp; temp[thisFileName] = group;
            files.push_back(temp);
        }

        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "readFiles");
        exit(1);
    }
}
/**************************************************************************************/

