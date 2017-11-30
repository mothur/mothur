//
//  mergecountcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/3/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//

#include "mergecountcommand.hpp"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> MergeCountCommand::setParameters(){
    try {
        CommandParameter pcount("count", "InputTypes", "", "", "", "", "","count",false,false,true); parameters.push_back(pcount);
        CommandParameter poutput("output", "String", "", "", "", "", "","",false,true,true); parameters.push_back(poutput);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeCountCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string MergeCountCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The merge.count command takes a list of count files separated by dashes and merges them into one file.";
        helpString += "The merge.count command parameters are count and output.";
        helpString += "Example merge.count(count=final.count_table-new.count_table, output=complete.count_table).";
        helpString += "Note: No spaces between parameter labels (i.e. output), '=' and parameters (i.e.yourOutputFileName).\n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeCountCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
MergeCountCommand::MergeCountCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["count"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeCountCommand", "MergeCountCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

MergeCountCommand::MergeCountCommand(string option)  {
    try {
        abort = false; calledHelp = false;
        
        if(option == "help") {
            help();
            abort = true; calledHelp = true;
        }else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            
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
            
            string fileList = validParameter.valid(parameters, "count");
            if(fileList == "not found") { m->mothurOut("[ERROR]: you must enter two or more count file names"); m->mothurOutEndLine();  abort=true;  }
            else{ 	util.splitAtDash(fileList, fileNames);	}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            string outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found")	{	outputDir = "";		}
            
            
            numInputFiles = fileNames.size();
            ifstream testFile;
            if(numInputFiles == 0){
                m->mothurOut("you must enter two or more file names and you entered " + toString(fileNames.size()) +  " file names"); m->mothurOutEndLine();
                abort=true;
            }
            else{
                for(int i=0;i<numInputFiles;i++){
                    if (inputDir != "") {
                        string path = util.hasPath(fileNames[i]);
                        //if the user has not given a path then, add inputdir. else leave path alone.
                        if (path == "") {	fileNames[i] = inputDir + fileNames[i];		}
                    }
                    
                    map<string, string> file; file["file"] = fileNames[i];
                    fileNames[i] = validParameter.validFile(file, "file");
                    if(fileNames[i] == "not found"){ 	abort = true;	}
                }
            }
            
            outputFileName = validParameter.valid(parameters, "output");
            if (outputFileName == "not found") { m->mothurOut("you must enter an output file name"); m->mothurOutEndLine();  abort=true;  }
            else if (outputDir != "") { outputFileName = outputDir + util.getSimpleName(outputFileName);  }
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeCountCommand", "MergeCountCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeCountCommand::execute(){
    try {
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        util.mothurRemove(outputFileName);
        
        //read headers from each file to confirm all contain groupinfo or all do not
        //Also collect all group names
        bool allContainGroups = true; bool allNoGroups = true;
        set<string> allGroups;
        for(int i = 0; i < numInputFiles; i++) {
            
            if (m->getControl_pressed()) { return 0; }
            
            vector<string> thisTablesGroups;
            CountTable table;
            bool hasGroups = table.testGroups(fileNames[i], thisTablesGroups);
            
            if (hasGroups) {
                allNoGroups = false;
                for (int j = 0; j < thisTablesGroups.size(); j++) { allGroups.insert(thisTablesGroups[j]);  }
            }else {  allContainGroups = false; }
        }
        int numGroups = allGroups.size();
        
        //check to make sure all files are one type - quit if not
        if (!allContainGroups && !allNoGroups) { m->mothurOut("[ERROR]: your have countfiles that contains group information and count files that do not. These cannot be combined without loss of information, please correct.\n"); m->setControl_pressed(true); return 0; }
        
        if (m->getControl_pressed()) { return 0; }
        
        //Create Blank Table - (set<string>&, map<string, string>&, set<string>&); //seqNames, seqName->group, groupNames
        set<string> seqNames; map<string, string> seqGroup; set<string> g;
        CountTable completeTable;
        completeTable.createTable(seqNames, seqGroup, g);
        
        //append first one to get headers
        map<string, int> groupIndex;
        if (allNoGroups) { util.appendBinaryFiles(fileNames[0], outputFileName);  }
        else { //create groupMap to save time setting abundance vector
            int count = 0;
            for (set<string>::iterator it = allGroups.begin(); it != allGroups.end(); it++) {
                completeTable.addGroup(*it);
                groupIndex[*it] = count; count++;
            }
        }
        
        //for each file
        for(int i = 0; i < numInputFiles; i++) {
            
            if (m->getControl_pressed()) { break; }

            if (allContainGroups) {
                
                CountTable table; table.readTable(fileNames[i], true, false);
                vector<string> groups = table.getNamesOfGroups();
                
                vector<string> seqs = table.getNamesOfSeqs();
                for (int j = 0; j < seqs.size(); j++) {
                    if (m->getControl_pressed()) { break; }
                    vector<int> abunds = table.getGroupCounts(seqs[j]);
                    vector<int> newAbunds; newAbunds.resize(numGroups, 0);
                    for (int k = 0; k < abunds.size(); k++) {
                        if (abunds[k] != 0) { //we need to set abundance in vector with all groups
                            //groups and abunds are in matching order. we know all groups are in groupIndex from above.
                            int newIndex = groupIndex[groups[k]];
                            newAbunds[newIndex] = abunds[k];
                        }
                    }
                    completeTable.push_back(seqs[j], newAbunds);
                }
            }
            else { util.appendFilesWithoutHeaders(fileNames[i], outputFileName); } //No group info so simple append
        }
        
        if (m->getControl_pressed()) {  util.mothurRemove(outputFileName); return 0;  }
        
        //print new table
        if (allContainGroups) {  completeTable.printTable(outputFileName); }
        
        if (m->getControl_pressed()) {  util.mothurRemove(outputFileName); return 0;  }
        
        //update current count file
        current->setCountFile(outputFileName);
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: "); m->mothurOutEndLine();
        m->mothurOut(outputFileName); m->mothurOutEndLine();	outputNames.push_back(outputFileName); outputTypes["merge"].push_back(outputFileName);
        m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeCountCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
