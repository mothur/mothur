/*
 *  sharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "makesharedcommand.h"
#include "counttable.h"

//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSharedRabunds(SharedRAbundVector* left, SharedRAbundVector* right){
    return (left->getGroup() < right->getGroup());
}
//**********************************************************************************************************************
vector<string> SharedCommand::setParameters(){
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "BiomListGroup", "BiomListGroup", "none","shared",false,false); parameters.push_back(pshared);
        CommandParameter pbiom("biom", "InputTypes", "", "", "BiomListGroup", "BiomListGroup", "none","shared",false,false); parameters.push_back(pbiom);
		CommandParameter plist("list", "InputTypes", "", "", "BiomListGroup", "BiomListGroup", "ListGroup","shared",false,false,true); parameters.push_back(plist);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "GroupCount", "none","",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "GroupCount", "ListGroup","",false,false,true); parameters.push_back(pgroup);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","group",false,false); parameters.push_back(pgroups);
        CommandParameter pzero("keepzeroes", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pzero);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["tshared"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["map"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
        
        abort = false; calledHelp = false; pickedGroups=false;
        allLines = true;

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
		helpString += "The make.shared command reads a list and group / count file or a biom file, or a shared file to convert or simply a count file and creates a shared file.\n";
		helpString += "The make.shared command parameters are list, group, biom, groups, count, shared and label. list and group or count are required unless a current file is available or you provide a biom file or you are converting a shared file.\n";
        helpString += "The count parameter allows you to provide a count file containing the group info for the list file. When the count file is provided without the list file, mothur will create a list and shared file for you.\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include, group names should be separated by dashes. ex. groups=A-B-C. Default is all groups in your groupfile.\n";
		helpString += "The label parameter is only valid with the list and group option and allows you to indicate which labels you want to include, label names should be separated by dashes. Default is all labels in your list file.\n";
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

        if (type == "shared") {  pattern = "[filename],shared-[filename],[distance],shared"; }
        else if (type == "tshared") {  pattern = "[filename],tshared-[filename],[distance],tshared"; }
        else if (type == "group") {  pattern = "[filename],[group],groups"; }
        else if (type == "list") {  pattern = "[filename],[distance],list"; }
        else if (type == "map") {  pattern = "[filename],map"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }

        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SharedCommand::SharedCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }

		else {
			 OptionParser parser(option, setParameters());
			 map<string, string> parameters = parser.getParameters();

			 ValidParameters validParameter;

			 //check for required parameters
			 listfile = validParameter.validFile(parameters, "list");
			 if (listfile == "not open") { listfile = ""; abort = true; }
			 else if (listfile == "not found") { listfile = "";  }
			 else { current->setListFile(listfile); }

             biomfile = validParameter.validFile(parameters, "biom");
             if (biomfile == "not open") { biomfile = ""; abort = true; }
             else if (biomfile == "not found") { biomfile = "";  }
             else { current->setBiomFile(biomfile); }
            
             sharedfile = validParameter.validFile(parameters, "shared");
             if (sharedfile == "not open") { sharedfile = ""; abort = true; }
             else if (sharedfile == "not found") { sharedfile = "";  }
             else { current->setSharedFile(sharedfile); }

			 ordergroupfile = validParameter.validFile(parameters, "ordergroup");
			 if (ordergroupfile == "not open") { abort = true; }
			 else if (ordergroupfile == "not found") { ordergroupfile = ""; }

			 groupfile = validParameter.validFile(parameters, "group");
			 if (groupfile == "not open") { groupfile = ""; abort = true; }
			 else if (groupfile == "not found") { groupfile = ""; }
			 else {  current->setGroupFile(groupfile); }

             countfile = validParameter.validFile(parameters, "count");
             if (countfile == "not open") { countfile = ""; abort = true; }
             else if (countfile == "not found") { countfile = ""; }
             else {
                 current->setCountFile(countfile);
                 CountTable temp;
                 if (!temp.testGroups(countfile)) {
                     m->mothurOut("\n[WARNING]: Your count file does not have group info, all reads will be assigned to mothurGroup.\n");
                     
                     temp.readTable(countfile, false, false); //dont read groups
                     map<string, int> seqs = temp.getNameMap();
                     
                     CountTable newCountTable;
                     newCountTable.addGroup("mothurGroup");
                     
                     for (map<string, int>::iterator it = seqs.begin(); it != seqs.end(); it++) {
                         vector<int> counts; counts.push_back(it->second);
                         newCountTable.push_back(it->first, counts);
                     }
                     
                     string newCountfileName = util.getRootName(countfile) + "mothurGroup" + util.getExtension(countfile);
                     newCountTable.printTable(newCountfileName);
                     
                     current->setCountFile(newCountfileName);
                     countfile = newCountfileName;
                     outputNames.push_back(newCountfileName);
                 }
             }

            if ((biomfile == "") && (listfile == "") && (countfile == "") && (sharedfile == "")) { //you must provide at least one of the following
				//is there are current file available for either of these?
				//give priority to list, then biom, then count
				listfile = current->getListFile();
				if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
				else {
					biomfile = current->getBiomFile();
                    if (biomfile != "") {  m->mothurOut("Using " + biomfile + " as input file for the biom parameter.\n"); }
					else {
                        countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                        else {
                            sharedfile = current->getSharedFile();
                            if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
                            else {
                                m->mothurOut("[ERROR]: No valid current files. You must provide a list, biom, shared or count file before you can use the make.shared command.\n");  abort = true;
                            }
                        }

					}
				}
			}
			else if ((biomfile != "") && (listfile != "")) { m->mothurOut("When executing a make.shared command you must enter ONLY ONE of the following: list or biom.\n"); abort = true; }

			if (listfile != "") {
				if ((groupfile == "") && (countfile == "")) {
					groupfile = current->getGroupFile();
					if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter.\n");  }
					else {
						countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter.\n"); }
                        else { m->mothurOut("[ERROR]: You need to provide a groupfile or countfile if you are going to use the list format.\n");  abort = true; }
					}
				}
			}


			 string groups = validParameter.valid(parameters, "groups");
			 if (groups == "not found") { groups = ""; }
			 else {
                 pickedGroups=true;
				 util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			 }

			 //check for optional parameter and set defaults
			 // ...at some point should added some additional type checking...
			 string label = validParameter.valid(parameters, "label");
			 if (label == "not found") { label = ""; }
			 else {
				 if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				 else { allLines = true;  }
			 }
            
            string temp = validParameter.valid(parameters, "keepzeroes");   if (temp == "not found"){    temp = "f";                }
            keepZeroes = util.isTrue(temp);
            
            if ((listfile == "") && (biomfile == "") && (countfile != "")) { //building a shared file from a count file, require label
                if (labels.size() == 0) { labels.insert("ASV"); }
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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

        if (listfile != "")         {  createSharedFromListGroup();     }
        else if (biomfile != "")    {  createSharedFromBiom();          }
        else if (sharedfile != "")  {  convertSharedFormat();           }
        else if ((listfile == "") && (countfile != "")) {  createSharedFromCount();      }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } }

		string currentName = "";
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
		}
        
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }

		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}

		m->mothurOut("\nOutput File Names:\n");
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
string SharedCommand::findFormat() {
    try {
        ifstream in; util.openInputFile(sharedfile, in);
        vector<string> headers; util.getline(in, headers);
        
        if (headers.size() > 4) { return "shared"; }
        else {
            if (headers.size() == 4) { //check to make sure this isn't a shared file with 1 OTU
                if (headers[3] == "abundance") { return "tshared"; }
            }else { m->mothurOut("[ERROR]: cannot determine format of shared file. Expected 4 or more columns, found " + toString(headers.size()) + "columns, please correct.\n"); m->setControl_pressed(true);  }
        }
        
        return "shared";
    }
    catch(exception& e) {
        m->errorOut(e, "SharedCommand", "findFormat");
        exit(1);
    }
}
//**********************************************************************************************************************
void SharedCommand::convertSharedFormat() {
    try {
        //getting output filename
        map<string, string> variables;
        if (outputdir == "") { outputdir += util.hasPath(sharedfile); }
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
        
        string tag = findFormat();
        
        if (m->getControl_pressed()) { return; }
        
        string sharedFilename = "";
        
        if (tag == "shared") { //converting shared to tshared
            tag = "tshared";
            sharedFilename = getOutputFileName(tag,variables);
            ofstream out; util.openOutputFile(sharedFilename, out);
            
            InputData input(sharedfile, "sharedfile", Groups);
            set<string> processedLabels;
            set<string> userLabels = labels;
            string lastLabel = "";
            
            SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            
            bool printHeaders = true;
            while (lookup != nullptr) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                lookup->printTidy(out, printHeaders, keepZeroes); delete lookup;
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
            out.close();
            
        }else { //tshared - converting tshared to shared
            tag = "shared";
            sharedFilename = getOutputFileName(tag,variables);
            ofstream out; util.openOutputFile(sharedFilename, out);
            
            ifstream inScan; util.openInputFile(sharedfile, inScan);
            
            util.getline(inScan); //read headers
            
            string label, group, otuName;
            int abundance;
            
            set<string> labels; set<string> groups; set<string> otuNames;
            
            while (!inScan.eof()) {
                
                if (m->getControl_pressed()) {  break; }
                
                inScan >> label >> group >> otuName >> abundance; util.gobble(inScan);
                
                labels.insert(label); groups.insert(group); otuNames.insert(otuName);
            }
            inScan.close();
            
            vector<string> oNames = util.mothurConvert(otuNames);
            vector<string> gNames = util.mothurConvert(groups); sort(gNames.begin(), gNames.end());
            int numGroups = gNames.size();
            int numOTUs = oNames.size();
        
            map<string, int> groupNameToIndex; //groupName -> index in otuAbunds
            for (int i = 0; i < numGroups; i++) { groupNameToIndex[gNames[i]] = i; }
            
            map<string, map<string, vector<int> > > sharedVectors;
            map<string, map<string, vector<int> > >::iterator itDistance;
            map<string, vector<int> >::iterator itOTU;
            map<string, int>::iterator itSample;
            for (set<string>::iterator it = labels.begin(); it != labels.end(); it++) { //for each distance
                
                map<string, vector<int> > otus; //otuName -> abunds
                
                for (int i = 0; i < numOTUs; i++) { //for each OTU, set all otus to 0
                    
                    vector<int> emptyOTU; emptyOTU.resize(numGroups, 0);
                    otus[oNames[i]] = emptyOTU; //add empty otu
                }
                sharedVectors[*it] = otus; //add empty vector
            }
            
            ifstream in; util.openInputFile(sharedfile, in);
            
            util.getline(in); //read headers
            
            while (!in.eof()) {
                
                if (m->getControl_pressed()) {  break; }
                
                in >> label >> group >> otuName >> abundance; util.gobble(in);
                
                itDistance = sharedVectors.find(label);
                
                if (itDistance != sharedVectors.end()) { //we have this label before - ie 0.03 or 0.05
                    
                    itOTU = (itDistance->second).find(otuName);
                    if (itOTU != (itDistance->second).end()) { //we have this otuName before - ie OTU0001 or OTU0234
                        
                        itSample = groupNameToIndex.find(group);
                        
                        if (itSample != groupNameToIndex.end()) { //we have this sample before - ie FD01 or FD03
                            
                            (itOTU->second)[itSample->second] = abundance;
                            
                        }else {
                            m->mothurOut("[ERROR]: Cannot find sample " + group + ", skipping.\n");
                        }
                    }else {
                        m->mothurOut("[ERROR]: Cannot find otu " + otuName + ", skipping.\n");
                    }
                }else {
                    m->mothurOut("[ERROR]: Cannot find label " + label + ", skipping.\n");
                }
            }
            in.close();
            
            bool printHeaders = true;
            //create sharedRabundVectors
            for (itDistance = sharedVectors.begin(); itDistance != sharedVectors.end(); itDistance++) { //for each distance
                
                //create empty shared vector with samples
                SharedRAbundVectors* shared = new SharedRAbundVectors();
                for (itSample = groupNameToIndex.begin(); itSample != groupNameToIndex.end(); itSample++) {
                    SharedRAbundVector* thisSample = new SharedRAbundVector();
                    thisSample->setGroup(itSample->first);
                    shared->push_back(thisSample);
                }
                
                shared->setLabels(itDistance->first); //set distance for shared vector
                m->mothurOut(itDistance->first+"\n");
                
                for (itOTU = (itDistance->second).begin(); itOTU != (itDistance->second).end(); itOTU++) { //for each OTU
                    shared->push_back(itOTU->second, itOTU->first); //add otus abundance
                }
                
                shared->eliminateZeroOTUS();
                shared->print(out, printHeaders);
                delete shared;
            }
            out.close();
        }
        
        outputNames.push_back(sharedFilename); outputTypes[tag].push_back(sharedFilename);
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedCommand", "convertSharedFormat");
        exit(1);
    }
}
//**********************************************************************************************************************
int SharedCommand::createSharedFromCount() {
    try {
        //getting output filename
        if (outputdir == "") { outputdir += util.hasPath(countfile); }
        string label = "ASV";
        if (labels.size() != 0) { label = *labels.begin(); }
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(countfile));
        variables["[distance]"] = "asv";
        string listFilename = getOutputFileName("list",variables);
        outputNames.push_back(listFilename); outputTypes["list"].push_back(listFilename);
        ofstream outlist; util.openOutputFile(listFilename, outlist);
        
        CountTable ct;  ct.readTable(countfile, true, false);
        map<string, int> counts = ct.getNameMap();
        
        ListVector list = ct.getListVector();
        list.setLabel(label);
        list.print(outlist, counts);
        outlist.close();
        
        listfile = listFilename;
        
        createSharedFromListGroup();
       
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedCommand", "createSharedFromCount");
        exit(1);
    }
}
//**********************************************************************************************************************
int SharedCommand::createSharedFromBiom() {
	try {
        //getting output filename
        string filename = biomfile;
		if (outputdir == "") { outputdir += util.hasPath(filename); }

        map<string, string> variables;
		variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(filename));
		filename = getOutputFileName("shared",variables);
		outputNames.push_back(filename); outputTypes["shared"].push_back(filename);

        ofstream out;
        util.openOutputFile(filename, out);

        /*{
            "id":"/Users/SarahsWork/Desktop/release/temp.job2.shared-unique",
            "format": "Biological Observation Matrix 0.9.1",
            "format_url": "http://biom-format.org",
            "type": "OTU table",
            "generated_by": "mothur1.24.0",
            "date": "Tue Apr 17 13:12:07 2012", */

        ifstream in; util.openInputFile(biomfile, in);

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
            if (m->getControl_pressed()) { break; }

            char c = in.get(); util.gobble(in);

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
        
        string biomType;
        map<string, string>::iterator it;
        it = fileLines.find("type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a type provided.\n"); }
        else {
            string thisLine = it->second;
            biomType = getTag(thisLine);
//            if ((biomType != "OTU table") && (biomType != "OTUtable") && (biomType != "Taxon table") && (biomType != "Taxontable")) { m->mothurOut("[ERROR]: " + biomType + " is not a valid biom type for mothur. Only types allowed are OTU table and Taxon table.\n"); m->setControl_pressed(true);  }
        }

        if (m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 0; }

        it = fileLines.find("matrix_type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a matrix_type provided.\n"); }
        else {
            string thisLine = it->second;
            matrixFormat = getTag(thisLine);
            if ((matrixFormat != "sparse") && (matrixFormat != "dense")) { m->mothurOut("[ERROR]: " + matrixFormat + " is not a valid biom matrix_type for mothur. Types allowed are sparse and dense.\n"); m->setControl_pressed(true); }
        }

        if (m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 0; }

        it = fileLines.find("matrix_element_type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a matrix_element_type provided.\n"); }
        else {
            string thisLine = it->second;
            matrixElementType = getTag(thisLine);
            if ((matrixElementType != "int") && (matrixElementType != "float")) { m->mothurOut("[ERROR]: " + matrixElementType + " is not a valid biom matrix_element_type for mothur. Types allowed are int and float.\n"); m->setControl_pressed(true); }
            if (matrixElementType == "float") { m->mothurOut("[WARNING]: the shared file only uses integers, any float values will be rounded down to the nearest integer.\n"); }
        }

        if (m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 0; }

        it = fileLines.find("rows");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a rows provided.\n"); }
        else {
            string thisLine = it->second;
            if ((biomType == "Taxon table") || (biomType == "Taxontable")) {
                string mapFilename = getOutputFileName("map",variables);
                outputNames.push_back(mapFilename); outputTypes["map"].push_back(mapFilename);
                ofstream outMap;
                util.openOutputFile(mapFilename, outMap);

                vector<string> taxonomies = readRows(thisLine, numRows);

                string snumBins = toString(numRows);
                for (int i = 0; i < numRows; i++) {

                    //if there is a bin label use it otherwise make one
                    string binLabel = "OTU";
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;

                    otuNames.push_back(binLabel);
                    outMap << otuNames[i] << '\t' << taxonomies[i] << endl;
                }
                outMap.close();
            }else{  otuNames = readRows(thisLine, numRows); }
        }

        if (m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 0; }

        it = fileLines.find("columns");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a columns provided.\n"); }
        else {
            string thisLine = it->second;

            //read sample names
            groupNames = readRows(thisLine, numCols);
            
            //if users selected groups, then remove the groups not wanted.
            if (Groups.size() == 0) { Groups = groupNames; }
            else { groupNames = Groups; }

            //set fileroot
            fileroot = outputdir + util.getRootName(util.getSimpleName(biomfile));
        }

        if (m->getControl_pressed()) {  out.close(); util.mothurRemove(filename); return 0; }

        it = fileLines.find("shape");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a shape provided.\n"); }
        else {
            string thisLine = it->second;
            getDims(thisLine, shapeNumRows, shapeNumCols);

            //check shape
            if (shapeNumCols != numCols) { m->mothurOut("[ERROR]: shape indicates " + toString(shapeNumCols) + " columns, but I only read " + toString(numCols) + " columns.\n"); m->setControl_pressed(true); }

            if (shapeNumRows != numRows) { m->mothurOut("[ERROR]: shape indicates " + toString(shapeNumRows) + " rows, but I only read " + toString(numRows) + " rows.\n"); m->setControl_pressed(true); }
        }

        if (m->getControl_pressed()) {  out.close(); util.mothurRemove(filename); return 0; }

        bool printHeaders = true;
        it = fileLines.find("data");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a data provided.\n"); }
        else {
            string thisLine = it->second;
            //read data
            SharedRAbundVectors* lookup = readData(matrixFormat, thisLine, matrixElementType, groupNames, otuNames.size());
            lookup->setOTUNames(otuNames);
            lookup->eliminateZeroOTUS();

            m->mothurOutEndLine(); m->mothurOut(lookup->getLabel()+"\n"); 
            printSharedData(lookup, out, printHeaders);
        }

        if (m->getControl_pressed()) {  util.mothurRemove(filename); return 0; }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "createSharedFromBiom");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundVectors* SharedCommand::readData(string matrixFormat, string line, string matrixElementType, vector<string>& groupNames, int numOTUs) {
	try {

        SharedRAbundVectors* lookup = new SharedRAbundVectors();
        
        //creates new sharedRAbunds
        for (int i = 0; i < groupNames.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(numOTUs); //sets all abunds to 0
            temp->setGroup(groupNames[i]);
            lookup->push_back(temp);
        }
        lookup->setLabels("userLabel");

        bool dataStart = false;
        bool inBrackets = false;
        string num = "";
        vector<int> nums;
        int otuCount = 0;
        for (int i = 0; i < line.length(); i++) {

            if (m->getControl_pressed()) { return lookup; }

            //look for opening [ to indicate data is starting
            if ((line[i] == '[') && (!dataStart)) { dataStart = true; i++;  if (!(i < line.length())) { break; } }
            else if ((line[i] == ']') && dataStart && (!inBrackets)) { break; } //we are done reading data

            if (dataStart) {
                if ((line[i] == '[') && (!inBrackets)) { inBrackets = true; i++;  if (!(i < line.length())) { break; } }
                else if ((line[i] == ']') && (inBrackets)) {
                    inBrackets = false;
                    int temp;
                    float temp2;
                    if (matrixElementType == "float") { util.mothurConvert(num, temp2); temp = floor(temp2); }
                    else { util.mothurConvert(num, temp); }
                    nums.push_back(temp);
                    num = "";

                    //save info to vectors
                    if (matrixFormat == "dense") {

                        //sanity check
                        if (nums.size() != lookup->getNumGroups()) { m->mothurOut("[ERROR]: trouble parsing OTU data.  OTU " + toString(otuCount) + " causing errors.\n"); m->setControl_pressed(true); }

                        //set abundances for this otu
                        //nums contains [abundSample0, abundSample1, abundSample2, ...] for current OTU
                        for (int j = 0; j < groupNames.size(); j++) { lookup->set(otuCount, nums[j], groupNames[j]); }

                        otuCount++;
                    }else {
                        //sanity check
                        if (nums.size() != 3) { m->mothurOut("[ERROR]: trouble parsing OTU data.\n"); m->setControl_pressed(true); }

                        //nums contains [otuNum, sampleNum, abundance]
                        lookup->set(nums[0], nums[2], groupNames[nums[1]]);
                    }
                    nums.clear();
                }

                if (inBrackets) {
                    if (line[i] == ',') {
                        int temp;
                        float temp2;
                        if (matrixElementType == "float") { util.mothurConvert(num, temp2); temp = floor(temp2); }
                        else { util.mothurConvert(num, temp); }
                        nums.push_back(temp);
                        num = "";
                    }else { if (!isspace(line[i])) { num += line[i]; }  }
                }
            }
        }

        if (pickedGroups) { lookup->eliminateZeroOTUS(); }

        return lookup;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "readData");
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
                util.mothurConvert(num, shapeNumCols);
                break;
            }

            if (inBar) {
                if (line[i] == ',') {
                    util.mothurConvert(num, shapeNumRows);
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

         "rows":[{"id": "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae", "metadata": null},
         {"id": "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae", "metadata": null}
         ....

         make look like above


         ],*/

        vector<string> names;
        int countOpenBrace = 0;
        int countClosedBrace = 0;
        int openParen = 0;
        int closeParen = 0;
        string nextRow = "";
        bool end = false;

        for (int i = 0; i < line.length(); i++) {

            if (m->getControl_pressed()) { return names; }

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
                util.splitAtChar(nextRow, items, ','); //parse by comma, will return junk for metadata but we aren't using that anyway
                string part = items[0]; items.clear();
                util.splitAtChar(part, items, ':'); //split part we want containing the ids
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
int SharedCommand::createSharedFromListGroup() {
	try {

        GroupMap* groupMap = nullptr;
        CountTable* countTable = nullptr;
        pickedGroups = false;
        if (groupfile != "") {
            groupMap = new GroupMap(groupfile);

            int groupError = groupMap->readMap();
            if (groupError == 1) { delete groupMap; return 0; }
            vector<string> allGroups = groupMap->getNamesOfGroups();
            if (Groups.size() == 0) { Groups = allGroups; }
            else { pickedGroups = true; }
        }else{
            countTable = new CountTable();
            countTable->readTable(countfile, true, false);
            vector<string> allGroups = countTable->getNamesOfGroups();
            if (Groups.size() == 0) { Groups = allGroups; }
            else { pickedGroups = true; }
        }
        int numGroups = Groups.size();
        if (m->getControl_pressed()) { return 0; }

        ofstream out;
        string filename = "";
        if (!pickedGroups) {
            string filename = listfile;
            if (outputdir == "") { outputdir += util.hasPath(filename); }

            map<string, string> variables;
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(filename));
            filename = getOutputFileName("shared",variables);
            outputNames.push_back(filename); outputTypes["shared"].push_back(filename);
            util.openOutputFile(filename, out);
        }

        //set fileroot
        fileroot = outputdir + util.getRootName(util.getSimpleName(listfile));
        map<string, string> variables;
		variables["[filename]"] = fileroot;
        string errorOff = "no error";

        InputData input(listfile, "shared", Groups);
        SharedListVector* SharedList = input.getSharedListVector();
        string lastLabel = SharedList->getLabel();
        SharedRAbundVectors* lookup;

        if (m->getControl_pressed()) {
            delete SharedList; if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }
            out.close(); if (!pickedGroups) { util.mothurRemove(filename); }
            return 0;
        }

        //sanity check
        vector<string> namesSeqs;
        int numGroupNames = 0;
        if (current->getGroupMode() == "group") { namesSeqs = groupMap->getNamesSeqs(); numGroupNames = groupMap->getNumSeqs(); }
        else { namesSeqs = countTable->getNamesOfSeqs(); numGroupNames = countTable->getNumUniqueSeqs(); }
        int error = ListGroupSameSeqs(namesSeqs, SharedList);

        if ((!pickedGroups) && (SharedList->getNumSeqs() != numGroupNames)) {  //if the user has not specified any groups and their files don't match exit with error
            m->mothurOut("Your group file contains " + toString(numGroupNames) + " sequences and list file contains " + toString(SharedList->getNumSeqs()) + " sequences. Please correct.\n");  m->setControl_pressed(true);

            out.close(); if (!pickedGroups) { util.mothurRemove(filename); } //remove blank shared file you made

            //delete memory
            delete SharedList; if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }
            return 0;
        }

        if (error == 1) { m->setControl_pressed(true); }

        //if user has specified groups make new groupfile for them
        if ((pickedGroups) && (current->getGroupMode() == "group")) { //make new group file
            string groups = "";
            if (numGroups < 4) {
                for (int i = 0; i < numGroups-1; i++) {
                    groups += Groups[i] + ".";
                }
                groups+=Groups[numGroups-1];
            }else { groups = "merge"; }
            map<string, string> variables;
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(listfile));
            variables["[group]"] = groups;
            string newGroupFile = getOutputFileName("group",variables);
            outputTypes["group"].push_back(newGroupFile);
            outputNames.push_back(newGroupFile);
            ofstream outGroups;
            util.openOutputFile(newGroupFile, outGroups);

            vector<string> names = groupMap->getNamesSeqs();
            string groupName;
            for (int i = 0; i < names.size(); i++) {
                groupName = groupMap->getGroup(names[i]);
                if (isValidGroup(groupName, Groups)) {
                    outGroups << names[i] << '\t' << groupName << endl;
                }
            }
            outGroups.close();
        }

        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        bool printHeaders = true;
    
        while((SharedList != nullptr) && ((allLines == 1) || (userLabels.size() != 0))) {
            if (m->getControl_pressed()) {
                delete SharedList; if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }
                if (!pickedGroups) { out.close(); util.mothurRemove(filename); }
                return 0;
            }

            if(allLines == 1 || labels.count(SharedList->getLabel()) == 1){

                lookup = SharedList->getSharedRAbundVector();

                m->mothurOut(lookup->getLabel()+"\n"); 

                if (m->getControl_pressed()) {
                    delete SharedList; if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }
                    delete lookup;
                    if (!pickedGroups) { out.close(); util.mothurRemove(filename); }
                    return 0;
                }

                //if picked groups must split the shared file by label
                if (pickedGroups) {
                    string filename = listfile;
                    if (outputdir == "") { outputdir += util.hasPath(filename); }

                    map<string, string> variables;
                    variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(filename));
                    variables["[distance]"] = lookup->getLabel();
                    filename = getOutputFileName("shared",variables);
                    outputNames.push_back(filename); outputTypes["shared"].push_back(filename);
                    ofstream out2;
                    util.openOutputFile(filename, out2);

                    lookup->eliminateZeroOTUS();
                    printSharedData(lookup, out2, printHeaders);
                    out2.close();
                }else {
                    printSharedData(lookup, out, printHeaders); //prints info to the .shared file
                }
                delete lookup;

                processedLabels.insert(SharedList->getLabel());
                userLabels.erase(SharedList->getLabel());
            }

            if ((util.anyLabelsToProcess(SharedList->getLabel(), userLabels, errorOff) ) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = SharedList->getLabel();

                delete SharedList;
                SharedList = input.getSharedListVector(lastLabel); //get new list vector to process

                lookup = SharedList->getSharedRAbundVector();
                m->mothurOut(lookup->getLabel()+"\n"); 

                if (m->getControl_pressed()) {
                    delete SharedList; if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }
                    delete lookup;
                    if (!pickedGroups) { out.close(); util.mothurRemove(filename); }
                    return 0;
                }

                //if picked groups must split the shared file by label
                if (pickedGroups) {
                    string filename = listfile;
                    if (outputdir == "") { outputdir += util.hasPath(filename); }

                    map<string, string> variables;
                    variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(filename));
                    variables["[distance]"] = lookup->getLabel();
                    filename = getOutputFileName("shared",variables);
                    outputNames.push_back(filename); outputTypes["shared"].push_back(filename);
                    ofstream out2;
                    util.openOutputFile(filename, out2);

                    lookup->eliminateZeroOTUS();
                    printSharedData(lookup, out2, printHeaders);
                    out2.close();
                }else {
                    printSharedData(lookup, out, printHeaders); //prints info to the .shared file
                }
                delete lookup;

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
        if (needToRun )  {
            if (SharedList != nullptr) {	delete SharedList;	}
            SharedList = input.getSharedListVector(lastLabel); //get new list vector to process

            lookup = SharedList->getSharedRAbundVector();
            m->mothurOut(lookup->getLabel()+"\n"); 

            if (m->getControl_pressed()) {
                if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }
                if (!pickedGroups) { out.close(); util.mothurRemove(filename); }
                return 0;
            }

            //if picked groups must split the shared file by label
            if (pickedGroups) {
                string filename = listfile;
                if (outputdir == "") { outputdir += util.hasPath(filename); }

                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(filename));
                variables["[distance]"] = lookup->getLabel();
                filename = getOutputFileName("shared",variables);
                outputNames.push_back(filename); outputTypes["shared"].push_back(filename);
                ofstream out2;
                util.openOutputFile(filename, out2);

                lookup->eliminateZeroOTUS();
                printSharedData(lookup, out2, printHeaders);
                out2.close();
            }else {
                printSharedData(lookup, out, printHeaders); //prints info to the .shared file
            }
            delete lookup;
            delete SharedList;
        }
        
        if (!pickedGroups) { out.close(); }

        if (groupMap != nullptr) { delete groupMap; } if (countTable != nullptr) { delete countTable; }

        if (m->getControl_pressed()) {
            if (!pickedGroups) { util.mothurRemove(filename); }
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
void SharedCommand::printSharedData(SharedRAbundVectors*& thislookup, ofstream& out, bool& printHeaders) {
	try {

		if (order.size() == 0) { //user has not specified an order so do aplabetically
            thislookup->print(out, printHeaders);
		}else{
			//create a map from groupName to each sharedrabund
			map<string, SharedRAbundVector*> myMap;
			map<string, SharedRAbundVector*>::iterator myIt;
            vector<SharedRAbundVector*> data = thislookup->getSharedRAbundVectors();

			for (int i = 0; i < data.size(); i++) { myMap[data[i]->getGroup()] = data[i]; }

			
			vector<string> Groups;

			//loop through ordered list and print the rabund
			for (int i = 0; i < order.size(); i++) {
				myIt = myMap.find(order[i]);

				if(myIt != myMap.end()) { //we found it
					out << (myIt->second)->getLabel() << '\t' << (myIt->second)->getGroup() << '\t';
					(myIt->second)->print(out);

					Groups.push_back((myIt->second)->getGroup());
				}else{
					m->mothurOut("Can't find shared info for " + order[i] + ", skipping.\n"); 
				}
			}

            for (int i = 0; i < data.size(); i++) { delete data[i]; } data.clear();
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
			if (m->getControl_pressed()) { return 0; }

			string names = SharedList->get(i);

			vector<string> listNames;
			util.splitAtComma(names, listNames);

			for (int j = 0; j < listNames.size(); j++) {
				int num = groupNamesSeqs.count(listNames[j]);

				if (num == 0) {
                    error = 1;
                    if (groupfile != "") {
                        m->mothurOut("[ERROR]: " + listNames[j] + " is in your listfile and not in your groupfile. Please correct.\n"); 	}
                    else{ m->mothurOut("[ERROR]: " + listNames[j] + " is in your listfile and not in your count file. Please correct.\n"); 	}
                }else {
                    groupNamesSeqs.erase(listNames[j]);
                    
                }
			}
		}

		for (set<string>::iterator itGroupSet = groupNamesSeqs.begin(); itGroupSet != groupNamesSeqs.end(); itGroupSet++) {
			error = 1;
			m->mothurOut("[ERROR]: " + (*itGroupSet) + " is in your groupfile and not your listfile. Please correct.\n"); 
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

		ifstream in; util.openInputFile(ordergroupfile, in);
		string thisGroup;

		while(!in.eof()){
			in >> thisGroup; util.gobble(in);

			order.push_back(thisGroup);

			if (m->getControl_pressed()) { order.clear(); break; }
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
