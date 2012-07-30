/*
 *  countseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "countseqscommand.h"
#include "groupmap.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> CountSeqsCommand::setParameters(){	
	try {
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
        CommandParameter plarge("large", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(plarge);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The count.seqs aka. make.table command reads a name file and outputs a .count.table file.  You may also provide a group file to get the counts broken down by group.\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include in the counts, by default all groups in your groupfile are used.\n";
        helpString += "The large parameter indicates the name and group files are too large to fit in RAM.\n";
		helpString += "When you use the groups parameter and a sequence does not represent any sequences from the groups you specify it is not included in the .count.summary file.\n";
		helpString += "The count.seqs command should be in the following format: count.seqs(name=yourNameFile).\n";
		helpString += "Example count.seqs(name=amazon.names) or make.table(name=amazon.names).\n";
		helpString += "Note: No spaces between parameter labels (i.e. name), '=' and parameters (i.e.yourNameFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountSeqsCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "counttable") {  outputFileName =  "count.table"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
CountSeqsCommand::CountSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["counttable"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "CountSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

CountSeqsCommand::CountSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["counttable"] = tempOutNames;
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){					
				namefile = m->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current namefile and the name parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { m->setGroupFile(groupfile); }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "all"; }
			m->splitAtDash(groups, Groups);
            
            string temp = validParameter.validFile(parameters, "large", false);		if (temp == "not found") {	temp = "F";	}
			large = m->isTrue(temp);
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(namefile);		}

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
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(namefile)) + getOutputFileNameTag("counttable");
		
        int total = 0;
        if (!large) { total = processSmall(outputFileName); }
        else { total = processLarge(outputFileName);  }
				
		if (m->control_pressed) { m->mothurRemove(outputFileName); return 0; }
		
        //set rabund file as new current rabundfile
		itTypes = outputTypes.find("counttable");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { string current = (itTypes->second)[0]; m->setCountTableFile(current); }
		}
        
        m->mothurOutEndLine();
		m->mothurOut("Total number of sequences: " + toString(total)); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountSeqsCommand::processSmall(string outputFileName){
	try {
        ofstream out;
        m->openOutputFile(outputFileName, out); outputTypes["counttable"].push_back(outputFileName);
        outputNames.push_back(outputFileName); outputTypes["counttable"].push_back(outputFileName);
		out << "Representative_Sequence\ttotal\t";
        
        GroupMap* groupMap;
		if (groupfile != "") { 
			groupMap = new GroupMap(groupfile); groupMap->readMap(); 
			
			//make sure groups are valid. takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			vector<string> nameGroups = groupMap->getNamesOfGroups();
			util->setGroups(Groups, nameGroups);
			delete util;
			
			//sort groupNames so that the group title match the counts below, this is needed because the map object automatically sorts
			sort(Groups.begin(), Groups.end());
			
			//print groupNames
			for (int i = 0; i < Groups.size(); i++) {
				out << Groups[i] << '\t';
			}
		}
		out << endl;
		
		//open input file
		ifstream in;
		m->openInputFile(namefile, in);
        
		int total = 0;
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol; m->gobble(in); in >> secondCol; m->gobble(in);
			
			vector<string> names;
			m->splitAtChar(secondCol, names, ',');
			
			if (groupfile != "") {
				//set to 0
				map<string, int> groupCounts;
				int total = 0;
				for (int i = 0; i < Groups.size(); i++) { groupCounts[Groups[i]] = 0; }
				
				//get counts for each of the users groups
				for (int i = 0; i < names.size(); i++) {
					string group = groupMap->getGroup(names[i]);
					
					if (group == "not found") { m->mothurOut("[ERROR]: " + names[i] + " is not in your groupfile, please correct."); m->mothurOutEndLine(); }
					else {
						map<string, int>::iterator it = groupCounts.find(group);
						
						//if not found, then this sequence is not from a group we care about
						if (it != groupCounts.end()) {
							it->second++;
							total++;
						}
					}
				}
				
				if (total != 0) {
					out << firstCol << '\t' << total << '\t';
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << it->second << '\t';
					}
					out << endl;
				}
			}else {
				out << firstCol << '\t' << names.size() << endl;
			}
			
			total += names.size();
		}
		in.close();
        out.close();
		
		if (groupfile != "") { delete groupMap; }

        return total;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "processSmall");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountSeqsCommand::processLarge(string outputFileName){
	try {
        set<string> namesOfGroups;
        map<string, int> initial;
        for (set<string>::iterator it = namesOfGroups.begin(); it != namesOfGroups.end(); it++) { initial[(*it)] = 0;  }
        ofstream out;
        m->openOutputFile(outputFileName, out); 
        outputNames.push_back(outputFileName); outputTypes["counttable"].push_back(outputFileName);
		out << "Representative_Sequence\ttotal\t";
        if (groupfile == "") { out << endl; }
        
        map<string, unsigned long long> namesToIndex;
        string outfile = m->getRootName(groupfile) + "sorted.groups.temp";
        string outName = m->getRootName(namefile) + "sorted.name.temp";
        map<int, string> indexToName;
        map<int, string> indexToGroup;
        if (groupfile != "") { 
            time_t estart = time(NULL);
            //convert name file to redundant -> unique.  set unique name equal to index so we can use vectors, save name for later.
            string newNameFile = m->getRootName(namefile) + ".name.temp";
            string newGroupFile = m->getRootName(groupfile) + ".group.temp";
            indexToName = processNameFile(newNameFile);
            indexToGroup = getGroupNames(newGroupFile, namesOfGroups);
            
            //sort file by first column so the names of sequences will be easier to find
            //use the unix sort 
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                string command = "sort -n " + newGroupFile + " -o " + outfile;
                system(command.c_str());
                command = "sort -n " + newNameFile + " -o " + outName;
                system(command.c_str());
            #else //sort using windows sort
                string command = "sort " + newGroupFile + " /O " + outfile;
                system(command.c_str());
                command = "sort " + newNameFile + " /O " + outName;
                system(command.c_str());
            #endif
            m->mothurRemove(newNameFile);
            m->mothurRemove(newGroupFile);
            
            m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to sort and index the group and name files. "); m->mothurOutEndLine();
        }else { outName = namefile; }
         
        time_t estart = time(NULL);
        //open input file
		ifstream in;
		m->openInputFile(outName, in);
        
        //open input file
		ifstream in2;
		
		int total = 0;
        vector< vector<int> > nameMapCount;
        if (groupfile != "") {
            m->openInputFile(outfile, in2);
            nameMapCount.resize(indexToName.size());
            for (int i = 0; i < nameMapCount.size(); i++) {
                nameMapCount[i].resize(indexToGroup.size(), 0);
            }
        }
        
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol;
			in >> firstCol;  m->gobble(in);
			
			if (groupfile != "") {
                int uniqueIndex;
                in >> uniqueIndex; m->gobble(in);
                
                string name; int groupIndex;
                in2 >> name >> groupIndex; m->gobble(in2);
                
                if (name != firstCol) { m->mothurOut("[ERROR]: found " + name + " in your groupfile, but " + firstCol + " was in your namefile, please correct.\n"); m->control_pressed = true; }
                
                nameMapCount[uniqueIndex][groupIndex]++;
                total++;
            }else { 
                string secondCol;
                in >> secondCol; m->gobble(in);
                int num = m->getNumNames(secondCol);
                out << firstCol << '\t' << num << endl;
                total += num;
            }
		}
		in.close();
        
        if (groupfile != "") {
            m->mothurRemove(outfile);
            m->mothurRemove(outName);
            in2.close();
            for (map<int, string>::iterator it = indexToGroup.begin(); it != indexToGroup.end(); it++) { out << it->second << '\t';  }
            out << endl;
            for (int i = 0; i < nameMapCount.size(); i++) {
                string totalsLine = "";
                int seqTotal = 0;
                for (int j = 0; j < nameMapCount[i].size(); j++) {
                    seqTotal += nameMapCount[i][j];
                    totalsLine += toString(nameMapCount[i][j]) + '\t';
                }
                out << indexToName[i] << '\t' << seqTotal << '\t' << totalsLine << endl;
            }
        }
        
        out.close();
        
        m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to create the count table file. "); m->mothurOutEndLine();
        
        return total;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "processLarge");
		exit(1);
	}
}
/**************************************************************************************************/
map<int, string> CountSeqsCommand::processNameFile(string name) {
	try {
        map<int, string> indexToNames;
        
        ofstream out;
        m->openOutputFile(name, out);
        
        //open input file
		ifstream in;
		m->openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        int count = 0;
        
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    //parse names into vector
                    vector<string> theseNames;
                    m->splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  out << theseNames[i] << '\t' << count << endl;  }
                    indexToNames[count] = firstCol;
                    pairDone = false; 
                    count++;
                }
            }
		}
		in.close();
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
        
        ofstream out;
        m->openOutputFile(filename, out);
        
        //open input file
		ifstream in;
		m->openInputFile(groupfile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        int count = 0;
        
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
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



