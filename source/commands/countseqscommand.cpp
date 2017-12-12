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
		CommandParameter pname("name", "InputTypes", "", "", "NameSHared", "NameSHared", "none","count",false,false,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "sharedGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The groups parameter allows you to indicate which groups you want to include in the counts, by default all groups in your groupfile are used.\n";
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
CountSeqsCommand::CountSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["count"] = tempOutNames;
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
        allLines = 1;
		
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
			outputTypes["count"] = tempOutNames;
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){	namefile = ""; }
            else { current->setNameFile(namefile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found"){	sharedfile = ""; }
            else { current->setSharedFile(sharedfile); }
            
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { current->setGroupFile(groupfile); }
            
            if ((namefile == "") && (sharedfile == "")) {
                namefile = current->getNameFile();
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else {
                    sharedfile = current->getSharedFile();
                    if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
                    else {
                        m->mothurOut("You have no current namefile or sharedfile and the name or shared parameter is required."); m->mothurOutEndLine(); abort = true;
                    }
                }
			}

			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "all"; }
			util.splitAtDash(groups, Groups);
            if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
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

        if (namefile != "") {
            if (outputDir == "") { outputDir = util.hasPath(namefile); }
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(namefile));
            string outputFileName = getOutputFileName("count", variables);
            
            long start = time(NULL);
            unsigned long long total = process(outputFileName);
            
            if (m->getControl_pressed()) { util.mothurRemove(outputFileName); return 0; }
            
            m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to create a table for " + toString(total) + " sequences.\n\n");
            m->mothurOut("Total number of sequences: " + toString(total) + "\n");
        }else {
            if (outputDir == "") { outputDir = util.hasPath(sharedfile); }
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
            
            InputData input(sharedfile, "sharedfile", Groups);
            SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
            Groups = lookup->getNamesGroups();
            string lastLabel = lookup->getLabel();
            vector<string> currentLabels = lookup->getOTUNames();
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->getControl_pressed()) { delete lookup; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
                
                if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
                    
                    m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                    processShared(data, variables, currentLabels);
                    for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear();
                    
                    processedLabels.insert(lookup->getLabel());
                    userLabels.erase(lookup->getLabel());
                }
                
                if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = lookup->getLabel();
                    
                    delete lookup;
                    lookup = input.getSharedRAbundVectors(lastLabel);
                    m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    
                    vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                    processShared(data, variables, currentLabels);
                    for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear();
                    
                    processedLabels.insert(lookup->getLabel());
                    userLabels.erase(lookup->getLabel());
                    
                    //restore real lastlabel to save below
                    lookup->setLabels(saveLabel);
                }
                
                lastLabel = lookup->getLabel();
                //prevent memory leak
                delete lookup;
                
                if (m->getControl_pressed()) { return 0; }
                
                //get next line to process
                lookup = input.getSharedRAbundVectors();
            }
            
            if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  return 0; }
            
            //output error messages about any remaining user labels
            set<string>::iterator it;
            bool needToRun = false;
            for (it = userLabels.begin(); it != userLabels.end(); it++) {
                m->mothurOut("Your file does not include the label " + *it);
                if (processedLabels.count(lastLabel) != 1) {
                    m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
                    needToRun = true;
                }else {
                    m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
                }
            }
            
            //run last label if you need to
            if (needToRun )  {
                delete lookup;
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                
                vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                processShared(data, variables, currentLabels);
                for(int i = 0; i < data.size(); i++) {  delete data[i]; } data.clear();
                
               delete lookup;
            }
            
        }
        
        //set rabund file as new current rabundfile
		itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { string currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
        
        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
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
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        out << "OTU_Label\ttotal";
        for (int i = 0; i < lookup.size(); i++) { out << '\t' << lookup[i]->getGroup(); } out << endl;
        
        for (int j = 0; j < lookup[0]->getNumBins(); j++) {
            if (m->getControl_pressed()) { break; }
            
            int total = 0;
            string output = "";
            for (int i = 0; i < lookup.size(); i++) {
                total += lookup[i]->get(j);
                output += '\t' + toString(lookup[i]->get(j));
            }
            out << currentLabels[j] << '\t' << total << output << endl;
        }
        out.close();
        
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
        ofstream out;
        util.openOutputFile(outputFileName, out); outputTypes["count"].push_back(outputFileName);
        outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
		out << "Representative_Sequence\ttotal";
        
        GroupMap* groupMap = NULL;
		if (groupfile != "") { 
			groupMap = new GroupMap(groupfile); groupMap->readMap(); 
			
			vector<string> nameGroups = groupMap->getNamesOfGroups();
            if (Groups.size() == 0) { Groups = nameGroups; }
            
			//sort groupNames so that the group title match the counts below, this is needed because the map object automatically sorts
			sort(Groups.begin(), Groups.end());
			
			//print groupNames
			for (int i = 0; i < Groups.size(); i++) { out << '\t' << Groups[i]; }
		}
		out << endl;
        out.close();
        
        unsigned long long total = driver(outputFileName, groupMap);
        
        if (groupfile != "") { delete groupMap; }
        
        return total;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "processSmall");
		exit(1);
	}
}
/**************************************************************************************************/
unsigned long long CountSeqsCommand::driver(string outputFileName, GroupMap*& groupMap) {
	try {
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        ifstream in;
		util.openInputFile(namefile, in);
		
        unsigned long long total = 0;
		while (!in.eof()) {
			if (m->getControl_pressed()) { break; }
			
			string firstCol, secondCol;
			in >> firstCol; util.gobble(in); in >> secondCol; util.gobble(in);
            
            util.checkName(firstCol);
            util.checkName(secondCol);
            
			vector<string> names;
			util.splitAtChar(secondCol, names, ',');
			
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
					out << firstCol << '\t' << total;
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << '\t' << it->second;
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
        
        return total;
    }
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
map<int, string> CountSeqsCommand::processNameFile(string name) {
	try {
        map<int, string> indexToNames;
        
        ofstream out;
        util.openOutputFile(name, out);
        
        //open input file
		ifstream in;
		util.openInputFile(namefile, in);
        
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
        
        ofstream out;
        util.openOutputFile(filename, out);
        
        //open input file
		ifstream in;
		util.openInputFile(groupfile, in);
        
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



