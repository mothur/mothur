/*
 *  countseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "countseqscommand.h"
#include "counttable.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> CountSeqsCommand::setParameters(){	
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "NameSHared-sharedGroup", "NameSHared", "none","count",false,false,true); parameters.push_back(pshared);
        CommandParameter pcount("count", "InputTypes", "", "", "NameSHared-sharedGroup", "NameSHared", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pname("name", "InputTypes", "", "", "NameSHared", "NameSHared", "none","count",false,false,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "sharedGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
        CommandParameter pcompress("compress", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pcompress);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
       
        vector<string> tempOutNames;
        outputTypes["count"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The count.seqs aka. make.table command reads a name or shared file and outputs a .count_table file.  You may also provide a group with the names file to get the counts broken down by group.\n";
        helpString += "You can also inflate or deflate an existing count table using the count and compress parameters. ie. count.seqs(count=current, compress=t)\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include in the counts, by default all groups in your groupfile are used.\n";
        helpString += "The compress parameter allows you to indicate you want the count table printed in compressed format. Default=t.\n";
		helpString += "When you use the groups parameter and a sequence does not represent any sequences from the groups you specify it is not included in the .count.summary file.\n";
		helpString += "The count.seqs command should be in the following format: count.seqs(name=yourNameFile).\n";
		helpString += "Example count.seqs(name=amazon.names) or make.table(name=amazon.names).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        if (type == "count") {  pattern = "[filename],count_table-[filename],[distance],count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "CountSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

CountSeqsCommand::CountSeqsCommand(string option) : Command()  {
	try {
        allLines = true;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){	namefile = ""; }
            else { current->setNameFile(namefile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found"){	sharedfile = ""; }
            else { current->setSharedFile(sharedfile); }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found"){	countfile = ""; }
            else { current->setCountFile(countfile); }
            
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { current->setGroupFile(groupfile); }
            
            if ((namefile == "") && (sharedfile == "") && (countfile == "")) {
                namefile = current->getNameFile();
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter.\n");  }
				else {
                    sharedfile = current->getSharedFile();
                    if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
                    else {
                        m->mothurOut("You have no current namefile or sharedfile and the name or shared parameter is required, unless inflating or deflating an existing count file.\n");  abort = true;
                    }
                }
			}

			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "all"; }
			util.splitAtDash(groups, Groups);
            if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
            
            string temp = validParameter.valid(parameters, "compress");			if (temp == "not found") { temp = "t"; }
            compress = util.isTrue(temp);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "CountSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountSeqsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        map<string, string> variables;

        if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, true, false, Groups);
            
            if (outputdir == "") { outputdir = util.hasPath(countfile); }
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(countfile));
            
            if (compress) {
                variables["[distance]"] = "sparse";
                string outputFileName = getOutputFileName("count", variables);
                outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
                
                ct.printCompressedTable(outputFileName);
            }else {
                variables["[distance]"] = "full";
                string outputFileName = getOutputFileName("count", variables);
                outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
                
                ct.printTable(outputFileName, false);
            }
        }else if (namefile != "") {
            if (outputdir == "") { outputdir = util.hasPath(namefile); }
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(namefile));
            string outputFileName = getOutputFileName("count", variables);
            
            long start = time(nullptr);
            unsigned long long total = process(outputFileName);
            
            if (m->getControl_pressed()) { util.mothurRemove(outputFileName); return 0; }
            
            m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " secs to create a table for " + toString(total) + " sequences.\n\n");
            m->mothurOut("Total number of sequences: " + toString(total) + "\n");
        }else {
            if (outputdir == "") { outputdir = util.hasPath(sharedfile); }
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
            
            InputData input(sharedfile, "sharedfile", Groups);
            set<string> processedLabels;
            set<string> userLabels = labels;
            string lastLabel = "";
            
            SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            vector<string> currentLabels = lookup->getOTUNames();
            Groups = lookup->getNamesGroups();
            
            while (lookup != nullptr) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                processShared(data, variables, currentLabels);
                for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear(); delete lookup;
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
        }
        
        //set rabund file as new current rabundfile
		itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { string currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
        
        m->mothurOut("\nOutput File Names: \n"); 
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
		m->mothurOutEndLine();
        
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

unsigned long long CountSeqsCommand::processShared(vector<SharedRAbundVector*>& lookup, map<string, string> variables, vector<string> currentLabels){
    try {
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("count", variables);
        outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
        
        CountTable ct;
        
        for (int i = 0; i < lookup.size(); i++) { ct.addGroup(lookup[i]->getGroup()); }
        
        for (int j = 0; j < lookup[0]->getNumBins(); j++) {
            if (m->getControl_pressed()) { break; }
            
            vector<int> outputs;
            for (int i = 0; i < lookup.size(); i++) {
                outputs.push_back(lookup[i]->get(j));
            }
            ct.push_back(currentLabels[j], outputs);
        }
        
        if (compress) {
            ct.printCompressedTable(outputFileName);
        }else {
            ct.printTable(outputFileName);
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountSeqsCommand", "processShared");
        exit(1);
    }
}
//**********************************************************************************************************************

unsigned long long CountSeqsCommand::process(string outputFileName){
	try {
        CountTable ct; ct.createTable(namefile, groupfile, Groups);

        if (compress) {
            ct.printCompressedTable(outputFileName);
        }else {
            ct.printTable(outputFileName, false);
        }
        
        outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
        
        return ct.getNumSeqs();
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
map<int, string> CountSeqsCommand::processNameFile(string name) {
	try {
        map<int, string> indexToNames;
        
        ofstream out; util.openOutputFile(name, out);
		ifstream in; util.openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        int count = 0;
        
		while (!in.eof()) {
			if (m->getControl_pressed()) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = util.splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    util.checkName(firstCol);
                    util.checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    util.splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  out << theseNames[i] << '\t' << count << endl;  }
                    indexToNames[count] = firstCol;
                    pairDone = false; 
                    count++;
                }
            }
		}
		in.close();
       
		
        if (rest != "") {
            vector<string> pieces = util.splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    util.checkName(firstCol);
                    util.checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    util.splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  out << theseNames[i] << '\t' << count << endl;  }
                    indexToNames[count] = firstCol;
                    pairDone = false; 
                    count++;
                }
            }

        }
        out.close();
        
        return indexToNames;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "processNameFile");
		exit(1);
	}
}
/**************************************************************************************************/
map<int, string> CountSeqsCommand::getGroupNames(string filename, set<string>& namesOfGroups) {
	try {
        map<int, string> indexToGroups;
        map<string, int> groupIndex;
        map<string, int>::iterator it;
        
        ofstream out; util.openOutputFile(filename, out);
		ifstream in; util.openInputFile(groupfile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        int count = 0;
        
		while (!in.eof()) {
			if (m->getControl_pressed()) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = util.splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    util.checkName(firstCol);
                    it = groupIndex.find(secondCol);
                    if (it == groupIndex.end()) { //add group, assigning the group and number so we can use vectors above
                        groupIndex[secondCol] = count;
                        count++;
                    }
                    out << firstCol << '\t' << groupIndex[secondCol] << endl; 
                    namesOfGroups.insert(secondCol);
                    pairDone = false; 
                }
            }
		}
		in.close();
        
        
        if (rest != "") {
            vector<string> pieces = util.splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    util.checkName(firstCol);
                    it = groupIndex.find(secondCol);
                    if (it == groupIndex.end()) { //add group, assigning the group and number so we can use vectors above
                        groupIndex[secondCol] = count;
                        count++;
                    }
                    out << firstCol << '\t' << groupIndex[secondCol] << endl; 
                    namesOfGroups.insert(secondCol);
                    pairDone = false; 
                }
            }
        }
        out.close();
		
        for (it = groupIndex.begin(); it != groupIndex.end(); it++) {  indexToGroups[it->second] = it->first;  }
        
        return indexToGroups;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getGroupNames");
		exit(1);
	}
}
//**********************************************************************************************************************



