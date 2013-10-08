/*
 *  sharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedcommand.h"
#include "sharedutilities.h"
#include "counttable.h"

//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSharedRabunds(SharedRAbundVector* left, SharedRAbundVector* right){
	return (left->getGroup() < right->getGroup());	
}
//**********************************************************************************************************************
vector<string> SharedCommand::setParameters(){	
	try {
        CommandParameter pbiom("biom", "InputTypes", "", "", "BiomListGroup", "BiomListGroup", "none","shared",false,false); parameters.push_back(pbiom);
		CommandParameter plist("list", "InputTypes", "", "", "BiomListGroup", "BiomListGroup", "ListGroup","shared",false,false,true); parameters.push_back(plist);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "GroupCount", "none","",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "GroupCount", "ListGroup","",false,false,true); parameters.push_back(pgroup);
		//CommandParameter pordergroup("ordergroup", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pordergroup);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","group",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SharedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.shared command reads a list and group file or a biom file and creates a shared file. If a list and group are provided a rabund file is created for each group.\n";
		helpString += "The make.shared command parameters are list, group, biom, groups, count and label. list and group or count are required unless a current file is available or you provide a biom file.\n";
        helpString += "The count parameter allows you to provide a count file containing the group info for the list file.\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include, group names should be separated by dashes. ex. groups=A-B-C. Default is all groups in your groupfile.\n";
		helpString += "The label parameter is only valid with the list and group option and allows you to indicate which labels you want to include, label names should be separated by dashes. Default is all labels in your list file.\n";
		//helpString += "The ordergroup parameter allows you to indicate the order of the groups in the sharedfile, by default the groups are listed alphabetically.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared") {  pattern = "[filename],shared"; } 
        else if (type == "rabund") {  pattern = "[filename],[group],rabund"; } 
        else if (type == "group") {  pattern = "[filename],[group],groups"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SharedCommand::SharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedCommand::SharedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			
			 vector<string> myArray = setParameters();
			 
			 OptionParser parser(option);
			 map<string, string> parameters = parser.getParameters();
			 
			 ValidParameters validParameter;
			 map<string, string>::iterator it;
			 
			 //check to make sure all parameters are valid for command
			 for (it = parameters.begin(); it != parameters.end(); it++) { 
				 if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			 }
			 
			 //if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			 if (inputDir == "not found"){	inputDir = "";		}
			 else {
				 string path;
				 it = parameters.find("list");
				 //user has given a template file
				 if(it != parameters.end()){ 
					 path = m->hasPath(it->second);
					 //if the user has not given a path then, add inputdir. else leave path alone.
					 if (path == "") {	parameters["list"] = inputDir + it->second;		}
				 }
			 
				 it = parameters.find("group");
				 //user has given a template file
				 if(it != parameters.end()){ 
					 path = m->hasPath(it->second);
					 //if the user has not given a path then, add inputdir. else leave path alone.
					 if (path == "") {	parameters["group"] = inputDir + it->second;		}
				 }
			 
				 it = parameters.find("count");
				 //user has given a template file
				 if(it != parameters.end()){ 
					 path = m->hasPath(it->second);
					 //if the user has not given a path then, add inputdir. else leave path alone.
					 if (path == "") {	parameters["count"] = inputDir + it->second;		}
				 }
                 
                 it = parameters.find("biom");
				 //user has given a template file
				 if(it != parameters.end()){ 
					 path = m->hasPath(it->second);
					 //if the user has not given a path then, add inputdir. else leave path alone.
					 if (path == "") {	parameters["biom"] = inputDir + it->second;		}
				 }
			 }
			 
             vector<string> tempOutNames;
             outputTypes["rabund"] = tempOutNames;
             outputTypes["shared"] = tempOutNames;
             outputTypes["group"] = tempOutNames;
			 
			 //if the user changes the output directory command factory will send this info to us in the output parameter 
			 outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			 
			 //check for required parameters
			 listfile = validParameter.validFile(parameters, "list", true);
			 if (listfile == "not open") { listfile = ""; abort = true; }
			 else if (listfile == "not found") { listfile = "";  }
			 else { m->setListFile(listfile); }	
            
             biomfile = validParameter.validFile(parameters, "biom", true);
             if (biomfile == "not open") { biomfile = ""; abort = true; }
             else if (biomfile == "not found") { biomfile = "";  }
             else { m->setBiomFile(biomfile); }		
							
			 ordergroupfile = validParameter.validFile(parameters, "ordergroup", true);
			 if (ordergroupfile == "not open") { abort = true; }	
			 else if (ordergroupfile == "not found") { ordergroupfile = ""; }
			 			 
			 groupfile = validParameter.validFile(parameters, "group", true);
			 if (groupfile == "not open") { groupfile = ""; abort = true; }	
			 else if (groupfile == "not found") { groupfile = ""; }
			 else {  m->setGroupFile(groupfile); }
            
             countfile = validParameter.validFile(parameters, "count", true);
             if (countfile == "not open") { countfile = ""; abort = true; }	
             else if (countfile == "not found") { countfile = ""; }
             else {  
                 m->setCountTableFile(countfile); 
                 CountTable temp;
                 if (!temp.testGroups(countfile)) { m->mothurOut("[ERROR]: Your count file does not have group info, aborting."); m->mothurOutEndLine(); abort=true; }
             }
			 
            if ((biomfile == "") && (listfile == "")) { 
				//is there are current file available for either of these?
				//give priority to list, then biom
				listfile = m->getListFile(); 
				if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					biomfile = m->getBiomFile(); 
					if (biomfile != "") {  m->mothurOut("Using " + biomfile + " as input file for the biom parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list or biom file before you can use the make.shared command."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			else if ((biomfile != "") && (listfile != "")) { m->mothurOut("When executing a make.shared command you must enter ONLY ONE of the following: list or biom."); m->mothurOutEndLine(); abort = true; }
			
			if (listfile != "") {
				if ((groupfile == "") && (countfile == "")) { 
					groupfile = m->getGroupFile(); 
					if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = m->getCountTableFile(); 
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a groupfile or countfile if you are going to use the list format."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
					}	
				}
			}

                        
			 string groups = validParameter.validFile(parameters, "groups", false);			
			 if (groups == "not found") { groups = ""; }
			 else { 
				 m->splitAtDash(groups, Groups);
				 m->setGroups(Groups);
			 }
			 
			 //check for optional parameter and set defaults
			 // ...at some point should added some additional type checking...
			 string label = validParameter.validFile(parameters, "label", false);			
			 if (label == "not found") { label = ""; }
			 else { 
				 if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				 else { allLines = 1;  }
			 }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SharedCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        		
		//getting output filename
        string filename = "";
		if (listfile != "") { filename = listfile; }
        else { filename = biomfile; }
		
		if (outputDir == "") { outputDir += m->hasPath(filename); }
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(filename));
		filename = getOutputFileName("shared",variables);
		outputNames.push_back(filename); outputTypes["shared"].push_back(filename);
		
        if (listfile != "") {  createSharedFromListGroup(filename);  }
        else {   createSharedFromBiom(filename);  }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  }
        
		//set rabund file as new current rabundfile
		string current = "";
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}	
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedCommand::createSharedFromBiom(string filename) {
	try {
        ofstream out;
        m->openOutputFile(filename, out);
        
        /*{
            "id":"/Users/SarahsWork/Desktop/release/temp.job2.shared-unique",
            "format": "Biological Observation Matrix 0.9.1",
            "format_url": "http://biom-format.org",
            "type": "OTU table",
            "generated_by": "mothur1.24.0",
            "date": "Tue Apr 17 13:12:07 2012", */
        
        ifstream in;
        m->openInputFile(biomfile, in);
         
        string matrixFormat = "";
        int numRows = 0;
        int numCols = 0;
        int shapeNumRows = 0;
        int shapeNumCols = 0;
        vector<string> otuNames;
        vector<string> groupNames;
        map<string, string> fileLines;
        vector<string> names;
        int countOpenBrace = 0;
        int countClosedBrace = 0;
        int openParen = -1; //account for opening brace
        int closeParen = 0;
        bool ignoreCommas = false;
        bool atComma = false;
        string line = "";
        string matrixElementType = "";
        
        while (!in.eof()) { //split file by tags, so each "line" will have something like "id":"/Users/SarahsWork/Desktop/release/final.tx.1.subsample.1.pick.shared-1"
            if (m->control_pressed) { break; }
            
            char c = in.get(); m->gobble(in);
            
            if (c == '[')               { countOpenBrace++;     }
            else if (c == ']')          { countClosedBrace++;   }
            else if (c == '{')          { openParen++;          }
            else if (c == '}')          { closeParen++;         }
            else if ((!ignoreCommas) && (c == ','))          { atComma = true;       }  
            
            if ((countOpenBrace != countClosedBrace) && (countOpenBrace != countClosedBrace)) { ignoreCommas = true;  }
            else if ((countOpenBrace == countClosedBrace) && (countOpenBrace == countClosedBrace)) { ignoreCommas = false;  }
            if (atComma && !ignoreCommas) { 
                if (fileLines.size() == 0) { //clip first {
                    line = line.substr(1);
                }
                string tag = getTag(line);
                fileLines[tag] = line;
                line = "";
                atComma = false;
                ignoreCommas = false;
                
            }else {  line += c;  }
            
        }
        if (line != "") {
            line = line.substr(0, line.length()-1);
            string tag = getTag(line);
            fileLines[tag] = line;
        }
        in.close();
        
        map<string, string>::iterator it;
        it = fileLines.find("type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a type provided.\n"); }
        else {
            string thisLine = it->second;
            string type = getTag(thisLine);
            if ((type != "OTU table") && (type != "OTUtable")) { m->mothurOut("[ERROR]: " + type + " is not a valid biom type for mothur. Only type allowed is OTU table.\n"); m->control_pressed = true;  }
        }
        
        if (m->control_pressed) { out.close(); m->mothurRemove(filename); return 0; }
        
        it = fileLines.find("matrix_type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a matrix_type provided.\n"); }
        else {
            string thisLine = it->second;
            matrixFormat = getTag(thisLine);
            if ((matrixFormat != "sparse") && (matrixFormat != "dense")) { m->mothurOut("[ERROR]: " + matrixFormat + " is not a valid biom matrix_type for mothur. Types allowed are sparse and dense.\n"); m->control_pressed = true; }
        }
        
        if (m->control_pressed) { out.close(); m->mothurRemove(filename); return 0; }
        
        it = fileLines.find("matrix_element_type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a matrix_element_type provided.\n"); }
        else {
            string thisLine = it->second;
            matrixElementType = getTag(thisLine);
            if ((matrixElementType != "int") && (matrixElementType != "float")) { m->mothurOut("[ERROR]: " + matrixElementType + " is not a valid biom matrix_element_type for mothur. Types allowed are int and float.\n"); m->control_pressed = true; }
            if (matrixElementType == "float") { m->mothurOut("[WARNING]: the shared file only uses integers, any float values will be rounded down to the nearest integer.\n"); }
        }
        
        if (m->control_pressed) { out.close(); m->mothurRemove(filename); return 0; }
        
        it = fileLines.find("rows");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a rows provided.\n"); }
        else {
            string thisLine = it->second;
            otuNames = readRows(thisLine, numRows);  
        }
        
        if (m->control_pressed) { out.close(); m->mothurRemove(filename); return 0; }
        
        it = fileLines.find("columns");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a columns provided.\n"); }
        else {
            string thisLine = it->second;
            //read sample names
            groupNames = readRows(thisLine, numCols); 
            
            //if users selected groups, then remove the groups not wanted.
            SharedUtil util;
            vector<string> Groups = m->getGroups();
            vector<string> allGroups = groupNames;
            util.setGroups(Groups, allGroups);
            m->setGroups(Groups);
            
            //fill filehandles with neccessary ofstreams
            int i;
            ofstream* temp;
            for (i=0; i<Groups.size(); i++) {
                temp = new ofstream;
                filehandles[Groups[i]] = temp;
            }
            
            //set fileroot
            fileroot = outputDir + m->getRootName(m->getSimpleName(biomfile));
            
            //clears file before we start to write to it below
            for (int i=0; i<Groups.size(); i++) {
                m->mothurRemove((fileroot + Groups[i] + ".rabund"));
                outputNames.push_back((fileroot + Groups[i] + ".rabund"));
                outputTypes["rabund"].push_back((fileroot + Groups[i] + ".rabund"));
            }
        }
        
        if (m->control_pressed) { for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { delete it3->second; } out.close(); m->mothurRemove(filename); return 0; }

        it = fileLines.find("shape");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a shape provided.\n"); }
        else {
            string thisLine = it->second;
            getDims(thisLine, shapeNumRows, shapeNumCols);
            
            //check shape
            if (shapeNumCols != numCols) { m->mothurOut("[ERROR]: shape indicates " + toString(shapeNumCols) + " columns, but I only read " + toString(numCols) + " columns.\n"); m->control_pressed = true; }
            
            if (shapeNumRows != numRows) { m->mothurOut("[ERROR]: shape indicates " + toString(shapeNumRows) + " rows, but I only read " + toString(numRows) + " rows.\n"); m->control_pressed = true; }
        }
        
        if (m->control_pressed) { for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { delete it3->second; } out.close(); m->mothurRemove(filename); return 0; }
        
        it = fileLines.find("data");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a data provided.\n"); }
        else {
            string thisLine = it->second;
            m->currentBinLabels = otuNames;
            
            //read data
            vector<SharedRAbundVector*> lookup = readData(matrixFormat, thisLine, matrixElementType, groupNames, otuNames.size());
            
            m->mothurOutEndLine(); m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
            lookup[0]->printHeaders(out); 
            printSharedData(lookup, out);
        }
        
        for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { delete it3->second; } 
        out.close();
        
        if (m->control_pressed) {  m->mothurRemove(filename); return 0; }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "createSharedFromBiom");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<SharedRAbundVector*> SharedCommand::readData(string matrixFormat, string line, string matrixElementType, vector<string>& groupNames, int numOTUs) {
	try {
        
        vector<SharedRAbundVector*> lookup; 
        
        //creates new sharedRAbunds
        for (int i = 0; i < groupNames.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(numOTUs); //sets all abunds to 0
            temp->setLabel("dummy");
            temp->setGroup(groupNames[i]);
            lookup.push_back(temp);
        }
        
        bool dataStart = false;
        bool inBrackets = false;
        string num = "";
        vector<int> nums;
        int otuCount = 0;
        for (int i = 0; i < line.length(); i++) {
            
            if (m->control_pressed) { return lookup; }
            
            //look for opening [ to indicate data is starting
            if ((line[i] == '[') && (!dataStart)) { dataStart = true; i++;  if (!(i < line.length())) { break; } }
            else if ((line[i] == ']') && dataStart && (!inBrackets)) { break; } //we are done reading data
                
            if (dataStart) {
                if ((line[i] == '[') && (!inBrackets)) { inBrackets = true; i++;  if (!(i < line.length())) { break; } }
                else if ((line[i] == ']') && (inBrackets)) { 
                    inBrackets = false; 
                    int temp;
                    float temp2;
                    if (matrixElementType == "float") { m->mothurConvert(num, temp2); temp = (int)temp2; }
                    else { m->mothurConvert(num, temp); }
                    nums.push_back(temp);
                    num = "";
                    
                    //save info to vectors
                    if (matrixFormat == "dense") {
                        
                        //sanity check
                        if (nums.size() != lookup.size()) { m->mothurOut("[ERROR]: trouble parsing OTU data.  OTU " + toString(otuCount) + " causing errors.\n"); m->control_pressed = true; }
                        
                        //set abundances for this otu
                        //nums contains [abundSample0, abundSample1, abundSample2, ...] for current OTU
                        for (int j = 0; j < lookup.size(); j++) { lookup[j]->set(otuCount, nums[j], groupNames[j]); }
                        
                        otuCount++;
                    }else {
                        //sanity check
                        if (nums.size() != 3) { m->mothurOut("[ERROR]: trouble parsing OTU data.\n"); m->control_pressed = true; }
                        
                        //nums contains [otuNum, sampleNum, abundance]
                        lookup[nums[1]]->set(nums[0], nums[2], groupNames[nums[1]]);
                    }
                    nums.clear();
                }
                
                if (inBrackets) {
                    if (line[i] == ',') {
                        int temp;
                        m->mothurConvert(num, temp);
                        nums.push_back(temp);
                        num = "";
                    }else { if (!isspace(line[i])) { num += line[i]; }  }
                }
            }
        }
        
        SharedUtil util;
        
		bool remove = false;
		for (int i = 0; i < lookup.size(); i++) {
			//if this sharedrabund is not from a group the user wants then delete it.
			if (util.isValidGroup(lookup[i]->getGroup(), m->getGroups()) == false) { 
				remove = true;
				delete lookup[i]; lookup[i] = NULL;
				lookup.erase(lookup.begin()+i); 
				i--; 
			}
		}
		
		if (remove) { eliminateZeroOTUS(lookup); }

        
        return lookup;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "readData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedCommand::eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup) {
    try {
        
        vector<SharedRAbundVector*> newLookup;
        for (int i = 0; i < thislookup.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
            temp->setLabel(thislookup[i]->getLabel());
            temp->setGroup(thislookup[i]->getGroup());
            newLookup.push_back(temp);
        }
        
        //for each bin
        vector<string> newBinLabels;
        string snumBins = toString(thislookup[0]->getNumBins());
        for (int i = 0; i < thislookup[0]->getNumBins(); i++) {
            if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
            
            //look at each sharedRabund and make sure they are not all zero
            bool allZero = true;
            for (int j = 0; j < thislookup.size(); j++) {
                if (thislookup[j]->getAbundance(i) != 0) { allZero = false;  break;  }
            }
            
            //if they are not all zero add this bin
            if (!allZero) {
                for (int j = 0; j < thislookup.size(); j++) {
                    newLookup[j]->push_back(thislookup[j]->getAbundance(i), thislookup[j]->getGroup());
                }
                
                //if there is a bin label use it otherwise make one
                string binLabel = "Otu";
                string sbinNumber = toString(i+1);
                if (sbinNumber.length() < snumBins.length()) { 
                    int diff = snumBins.length() - sbinNumber.length();
                    for (int h = 0; h < diff; h++) { binLabel += "0"; }
                }
                binLabel += sbinNumber; 
                if (i < m->currentBinLabels.size()) {  binLabel = m->currentBinLabels[i]; }
                
                newBinLabels.push_back(binLabel);
            }
        }
        
        for (int j = 0; j < thislookup.size(); j++) {  delete thislookup[j];  }
        
        thislookup = newLookup;
        m->currentBinLabels = newBinLabels;
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedCommand", "eliminateZeroOTUS");
        exit(1);
    }
}
//**********************************************************************************************************************
int SharedCommand::getDims(string line, int& shapeNumRows, int& shapeNumCols) {
	try {
        //get shape
        bool inBar = false;
        string num = "";
        
        for (int i = 0; i < line.length(); i++) {
            
            //you want to ignore any ; until you reach the next '
            if ((line[i] == '[') && (!inBar)) {  inBar = true; i++;  if (!(i < line.length())) { break; } } 
            else if ((line[i] == ']') && (inBar)) {  
                inBar= false;  
                m->mothurConvert(num, shapeNumCols);
                break;
            } 
            
            if (inBar) {  
                if (line[i] == ',') {
                    m->mothurConvert(num, shapeNumRows);
                    num = "";
                }else { if (!isspace(line[i])) { num += line[i]; }  }
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getDims");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SharedCommand::readRows(string line, int& numRows) {
	try {
        /*"rows":[
         {"id":"Otu01", "metadata":{"taxonomy":["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Porphyromonadaceae", "unclassified"], "bootstrap":[100, 100, 100, 100, 100, 100]}},
         {"id":"Otu02", "metadata":{"taxonomy":["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Rikenellaceae", "Alistipes"], "bootstrap":[100, 100, 100, 100, 100, 100]}},
         ...
         ],*/
        vector<string> names;
        int countOpenBrace = 0;
        int countClosedBrace = 0;
        int openParen = 0;
        int closeParen = 0;
        string nextRow = "";
        bool end = false;
        
        for (int i = 0; i < line.length(); i++) {
            
            if (m->control_pressed) { return names; }
            
            if (line[i] == '[')         { countOpenBrace++;     }
            else if (line[i] == ']')    { countClosedBrace++;   }
            else if (line[i] == '{')    { openParen++;          }
            else if (line[i] == '}')    { closeParen++;         }
            else if (openParen != 0)    { nextRow += line[i];   }  //you are reading the row info
            
            //you have reached the end of the rows info
            if ((countOpenBrace == countClosedBrace) && (countClosedBrace != 0)) { end = true; break; }
            if ((openParen == closeParen) && (closeParen != 0)) { //process row 
                numRows++;
                vector<string> items;
                m->splitAtChar(nextRow, items, ','); //parse by comma, will return junk for metadata but we aren't using that anyway
                string part = items[0]; items.clear();
                m->splitAtChar(part, items, ':'); //split part we want containing the ids
                string name = items[1];
                
                //remove "" if needed
                int pos = name.find("\"");
                if (pos != string::npos) {
                    string newName = "";
                    for (int k = 0; k < name.length(); k++) {
                        if (name[k] != '\"') { newName += name[k]; }
                    }
                    name = newName;
                }
                names.push_back(name);
                nextRow = "";
                openParen = 0;
                closeParen = 0;
            }
        }
        
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "readRows");
		exit(1);
	}
}
//**********************************************************************************************************************
//designed for things like "type": "OTU table", returns type 
string SharedCommand::getTag(string& line) {
	try {
        bool inQuotes = false;
        string tag = "";
        char c = '\"';
        
        for (int i = 0; i < line.length(); i++) {
        
            //you want to ignore any ; until you reach the next '
			if ((line[i] == c) && (!inQuotes)) {  inQuotes = true;  } 
			else if ((line[i] == c) && (inQuotes)) {  
                inQuotes= false;  
                line = line.substr(i+1);
                return tag;
            } 
            
			if (inQuotes) {  if (line[i] != c) { tag += line[i]; }  }
        }
        
        return tag;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getInfo");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedCommand::createSharedFromListGroup(string filename) {
	try {
        ofstream out;
        m->openOutputFile(filename, out);
        
        GroupMap* groupMap = NULL;
        CountTable* countTable = NULL;
        if (groupfile != "") {
            groupMap = new GroupMap(groupfile);
        
            int groupError = groupMap->readMap();
            if (groupError == 1) { delete groupMap; return 0; }
            vector<string> allGroups = groupMap->getNamesOfGroups();
            m->setAllGroups(allGroups);
        }else{
            countTable = new CountTable();
            countTable->readTable(countfile, true, false);
        }
        
        if (m->control_pressed) { return 0; }
        
        pickedGroups = false;
        
        //if hte user has not specified any groups then use them all
        if (Groups.size() == 0) {
            if (groupfile != "") { Groups = groupMap->getNamesOfGroups();  }
            else {  Groups = countTable->getNamesOfGroups();  }
            m->setGroups(Groups);
        }else { pickedGroups = true; }
        
        //fill filehandles with neccessary ofstreams
        int i;
        ofstream* temp;
        for (i=0; i<Groups.size(); i++) {
            temp = new ofstream;
            filehandles[Groups[i]] = temp;
        }
        
        //set fileroot
        fileroot = outputDir + m->getRootName(m->getSimpleName(listfile));
        map<string, string> variables; 
		variables["[filename]"] = fileroot;
        //clears file before we start to write to it below
        for (int i=0; i<Groups.size(); i++) {
            variables["[group]"] = Groups[i];
            string rabundFIleName = getOutputFileName("rabund",variables);
            m->mothurRemove(rabundFIleName);
            outputNames.push_back(rabundFIleName);
            outputTypes["rabund"].push_back(rabundFIleName);
        }
        
        string errorOff = "no error";
        
        //if user provided an order file containing the order the shared file should be in read it
        //if (ordergroupfile != "") { readOrderFile(); }
        
        InputData input(listfile, "shared");
        SharedListVector* SharedList = input.getSharedListVector();
        string lastLabel = SharedList->getLabel();
        vector<SharedRAbundVector*> lookup; 
        
        if (m->control_pressed) { 
            delete SharedList; if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
            for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
            out.close(); m->mothurRemove(filename); 
            for (int i=0; i<Groups.size(); i++) {  variables["[group]"] = Groups[i];
                string rabundFIleName = getOutputFileName("rabund",variables);
                m->mothurRemove(rabundFIleName);		}
            return 0; 
        }
        
        //sanity check
        vector<string> namesSeqs;
        int numGroupNames = 0;
        if (m->groupMode == "group") { namesSeqs = groupMap->getNamesSeqs(); numGroupNames = groupMap->getNumSeqs(); }
        else { namesSeqs = countTable->getNamesOfSeqs(); numGroupNames = countTable->getNumUniqueSeqs(); }
        int error = ListGroupSameSeqs(namesSeqs, SharedList);
        
        if ((!pickedGroups) && (SharedList->getNumSeqs() != numGroupNames)) {  //if the user has not specified any groups and their files don't match exit with error
            m->mothurOut("Your group file contains " + toString(numGroupNames) + " sequences and list file contains " + toString(SharedList->getNumSeqs()) + " sequences. Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
            
            out.close(); m->mothurRemove(filename); //remove blank shared file you made
            
            //delete memory
            for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { delete it3->second; }
            delete SharedList; if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
            return 0; 
        }
        
        if (error == 1) { m->control_pressed = true; }
        
        //if user has specified groups make new groupfile for them
        if ((pickedGroups) && (m->groupMode == "group")) { //make new group file
            string groups = "";
            if (m->getNumGroups() < 4) {
                for (int i = 0; i < m->getNumGroups(); i++) {
                    groups += (m->getGroups())[i] + ".";
                }
            }else { groups = "merge"; }
            map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listfile));
            variables["[group]"] = groups;
            string newGroupFile = getOutputFileName("group",variables);
            outputTypes["group"].push_back(newGroupFile); 
            outputNames.push_back(newGroupFile);
            ofstream outGroups;
            m->openOutputFile(newGroupFile, outGroups);
            
            vector<string> names = groupMap->getNamesSeqs();
            string groupName;
            for (int i = 0; i < names.size(); i++) {
                groupName = groupMap->getGroup(names[i]);
                if (isValidGroup(groupName, m->getGroups())) {
                    outGroups << names[i] << '\t' << groupName << endl;
                }
            }
            outGroups.close();
        }
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;	
        
        while((SharedList != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            if (m->control_pressed) { 
                delete SharedList; if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
                for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
                out.close(); m->mothurRemove(filename); 
                for (int i=0; i<Groups.size(); i++) {  variables["[group]"] = Groups[i];
                    string rabundFIleName = getOutputFileName("rabund",variables);
                    m->mothurRemove(rabundFIleName);		}
                return 0; 
            }
            
            if(allLines == 1 || labels.count(SharedList->getLabel()) == 1){
                
                lookup = SharedList->getSharedRAbundVector();
                
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                if (pickedGroups) { //check for otus with no seqs in them
                    eliminateZeroOTUS(lookup);
                }
                
                if (m->control_pressed) { 
                    delete SharedList; if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
                    for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                    for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
                    out.close(); m->mothurRemove(filename); 
                    for (int i=0; i<Groups.size(); i++) {  variables["[group]"] = Groups[i];
                        string rabundFIleName = getOutputFileName("rabund",variables);
                        m->mothurRemove(rabundFIleName);		}
                    return 0; 
                }
                
                if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
                printSharedData(lookup, out); //prints info to the .shared file
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                
                processedLabels.insert(SharedList->getLabel());
                userLabels.erase(SharedList->getLabel());
            }
            
            if ((m->anyLabelsToProcess(SharedList->getLabel(), userLabels, errorOff) == true) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = SharedList->getLabel();
                
                delete SharedList;
                SharedList = input.getSharedListVector(lastLabel); //get new list vector to process
                
                lookup = SharedList->getSharedRAbundVector();
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                if (pickedGroups) { //check for otus with no seqs in them
                    eliminateZeroOTUS(lookup);
                }
                
                
                if (m->control_pressed) { 
                    delete SharedList; if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
                    for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                    for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
                    out.close(); m->mothurRemove(filename); 
                    for (int i=0; i<Groups.size(); i++) {  variables["[group]"] = Groups[i];
                        string rabundFIleName = getOutputFileName("rabund",variables);
                        m->mothurRemove(rabundFIleName);		}
                    return 0; 
                }
                
                if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
                printSharedData(lookup, out); //prints info to the .shared file
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                
                processedLabels.insert(SharedList->getLabel());
                userLabels.erase(SharedList->getLabel());
                
                //restore real lastlabel to save below
                SharedList->setLabel(saveLabel);
            }
            
            
            lastLabel = SharedList->getLabel();
            
            delete SharedList;
            SharedList = input.getSharedListVector(); //get new list vector to process
        }
        
        //output error messages about any remaining user labels
        set<string>::iterator it;
        bool needToRun = false;
        for (it = userLabels.begin(); it != userLabels.end(); it++) {  
            if (processedLabels.count(lastLabel) != 1) {
                needToRun = true;
            }
        }
        
        //run last label if you need to
        if (needToRun == true)  {
            if (SharedList != NULL) {	delete SharedList;	}
            SharedList = input.getSharedListVector(lastLabel); //get new list vector to process
            
            lookup = SharedList->getSharedRAbundVector();
            m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
            if (pickedGroups) { //check for otus with no seqs in them
                eliminateZeroOTUS(lookup);
            }
            
            if (m->control_pressed) { 
                if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
                for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;   }
                out.close(); m->mothurRemove(filename); 
                for (int i=0; i<Groups.size(); i++) {  variables["[group]"] = Groups[i];
                    string rabundFIleName = getOutputFileName("rabund",variables);
                    m->mothurRemove(rabundFIleName);		}
                return 0; 
            }
            
            if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
            printSharedData(lookup, out); //prints info to the .shared file
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
            delete SharedList;
        }
        
        out.close();
        
        for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
            delete it3->second;
        }
        
        if (groupMap != NULL) { delete groupMap; } if (countTable != NULL) { delete countTable; }
		
        if (m->control_pressed) { 
            m->mothurRemove(filename); 
            for (int i=0; i<Groups.size(); i++) {  variables["[group]"] = Groups[i];
                string rabundFIleName = getOutputFileName("rabund",variables);
                m->mothurRemove(rabundFIleName);		}
            return 0; 
        }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "createSharedFromListGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
void SharedCommand::printSharedData(vector<SharedRAbundVector*> thislookup, ofstream& out) {
	try {
		
		if (order.size() == 0) { //user has not specified an order so do aplabetically
			sort(thislookup.begin(), thislookup.end(), compareSharedRabunds);
			
			m->clearGroups();
			vector<string> Groups;
			
			//initialize bin values
			for (int i = 0; i < thislookup.size(); i++) {
				out << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() << '\t';
				thislookup[i]->print(out);
				
				Groups.push_back(thislookup[i]->getGroup());
				
				RAbundVector rav = thislookup[i]->getRAbundVector();
                map<string, string> variables; 
                variables["[filename]"] = fileroot;
                variables["[group]"] = thislookup[i]->getGroup();
				m->openOutputFileAppend(getOutputFileName("rabund",variables), *(filehandles[thislookup[i]->getGroup()]));
				rav.print(*(filehandles[thislookup[i]->getGroup()]));
				(*(filehandles[thislookup[i]->getGroup()])).close();
			}
			m->setGroups(Groups);
		}else{
			//create a map from groupName to each sharedrabund
			map<string, SharedRAbundVector*> myMap;
			map<string, SharedRAbundVector*>::iterator myIt;
			
			for (int i = 0; i < thislookup.size(); i++) {
				myMap[thislookup[i]->getGroup()] = thislookup[i];
			}
			
			m->clearGroups();
			vector<string> Groups;
			
			//loop through ordered list and print the rabund
			for (int i = 0; i < order.size(); i++) {
				myIt = myMap.find(order[i]);
				
				if(myIt != myMap.end()) { //we found it
					out << (myIt->second)->getLabel() << '\t' << (myIt->second)->getGroup() << '\t';
					(myIt->second)->print(out);
					
					Groups.push_back((myIt->second)->getGroup());
				
					RAbundVector rav = (myIt->second)->getRAbundVector();
                    map<string, string> variables; 
                    variables["[filename]"] = fileroot;
                    variables["[group]"] = (myIt->second)->getGroup();
					m->openOutputFileAppend(getOutputFileName("rabund",variables), *(filehandles[(myIt->second)->getGroup()]));
					rav.print(*(filehandles[(myIt->second)->getGroup()]));
					(*(filehandles[(myIt->second)->getGroup()])).close();
				}else{
					m->mothurOut("Can't find shared info for " + order[i] + ", skipping."); m->mothurOutEndLine();
				}
			}
			
			m->setGroups(Groups);
		
		}
 
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "printSharedData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedCommand::ListGroupSameSeqs(vector<string>& groupMapsSeqs, SharedListVector* SharedList) {
	try {
		int error = 0; 
       	
		set<string> groupNamesSeqs;
		for(int i = 0; i < groupMapsSeqs.size(); i++) {
			groupNamesSeqs.insert(groupMapsSeqs[i]);
		}
		
		//go through list and if group returns "not found" output it
		for (int i = 0; i < SharedList->getNumBins(); i++) {
			if (m->control_pressed) { return 0; } 
			
			string names = SharedList->get(i); 
			
			vector<string> listNames;
			m->splitAtComma(names, listNames);
			
			for (int j = 0; j < listNames.size(); j++) {
				int num = groupNamesSeqs.count(listNames[j]);
				
				if (num == 0) { 
                    error = 1; 
                    if (groupfile != "") { 
                        m->mothurOut("[ERROR]: " + listNames[j] + " is in your listfile and not in your groupfile. Please correct."); m->mothurOutEndLine();	} 
                    else{ m->mothurOut("[ERROR]: " + listNames[j] + " is in your listfile and not in your count file. Please correct."); m->mothurOutEndLine();	}
                }else { groupNamesSeqs.erase(listNames[j]); }
			}
		}
		
		for (set<string>::iterator itGroupSet = groupNamesSeqs.begin(); itGroupSet != groupNamesSeqs.end(); itGroupSet++) {
			error = 1; 
			m->mothurOut("[ERROR]: " + (*itGroupSet) + " is in your groupfile and not your listfile. Please correct."); m->mothurOutEndLine();
		}
				
		return error;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "ListGroupSameSeqs");
		exit(1);
	}
}
//**********************************************************************************************************************

SharedCommand::~SharedCommand(){
	//delete list;
	
	
}
//**********************************************************************************************************************
int SharedCommand::readOrderFile() {
	try {
		//remove old names
		order.clear();
		
		ifstream in;
		m->openInputFile(ordergroupfile, in);
		string thisGroup;
		
		while(!in.eof()){
			in >> thisGroup; m->gobble(in);
						
			order.push_back(thisGroup);
			
			if (m->control_pressed) { order.clear(); break; }
		}
		in.close();		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "readOrderFile");
		exit(1);
	}
}
//**********************************************************************************************************************

bool SharedCommand::isValidGroup(string groupname, vector<string> groups) {
	try {
		for (int i = 0; i < groups.size(); i++) {
			if (groupname == groups[i]) { return true; }
		}
		
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "isValidGroup");
		exit(1);
	}
}
/************************************************************/


