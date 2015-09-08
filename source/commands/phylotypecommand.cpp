/*
 *  phylotypecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/20/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "phylotypecommand.h"
#include "phylotree.h"
#include "listvector.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> PhylotypeCommand::setParameters(){	
	try {
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","list-rabund-sabund",false,true,true); parameters.push_back(ptaxonomy);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName","rabund-sabund",false,false,true); parameters.push_back(pname);
		CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pcutoff("cutoff", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PhylotypeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The phylotype command reads a taxonomy file and outputs a .list, .rabund and .sabund file. \n";
		helpString += "The phylotype command parameter options are taxonomy, name, count, cutoff and label. The taxonomy parameter is required.\n";
		helpString += "The cutoff parameter allows you to specify the level you want to stop at.  The default is the highest level in your taxonomy file. \n";
		helpString += "For example: taxonomy = Bacteria;Bacteroidetes-Chlorobi;Bacteroidetes; - cutoff=2, would truncate the taxonomy to Bacteria;Bacteroidetes-Chlorobi; \n";
		helpString += "For the cutoff parameter levels count up from the root of the phylotree. This enables you to look at the grouping down to a specific resolution, say the genus level.\n";
		helpString += "The label parameter allows you to specify which level you would like, and are separated by dashes.  The default all levels in your taxonomy file. \n";
		helpString += "For the label parameter, levels count down from the root to keep the output similar to mothur's other commands which report information from finer resolution to coarser resolutions.\n";
		helpString += "The phylotype command should be in the following format: \n";
		helpString += "phylotype(taxonomy=yourTaxonomyFile, cutoff=yourCutoff, label=yourLabels) \n";
		helpString += "Eaxample: phylotype(taxonomy=amazon.taxonomy, cutoff=5, label=1-3-5).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string PhylotypeCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; }
        else if (type == "rabund") {  pattern = "[filename],[clustertag],rabund"; }
        else if (type == "sabund") {  pattern = "[filename],[clustertag],sabund"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PhylotypeCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
PhylotypeCommand::PhylotypeCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "PhylotypeCommand");
		exit(1);
	}
}
/**********************************************************************************************************************/
PhylotypeCommand::PhylotypeCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			taxonomyFileName = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyFileName == "not found") { 
				taxonomyFileName = m->getTaxonomyFile(); 
				if (taxonomyFileName != "") {  m->mothurOut("Using " + taxonomyFileName + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("No valid current files. taxonomy is a required parameter."); m->mothurOutEndLine(); 
					abort = true; 
				}
			}else if (taxonomyFileName == "not open") { taxonomyFileName = ""; abort = true; }	
			else { m->setTaxonomyFile(taxonomyFileName); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { readNamesFile(); m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(taxonomyFileName); //if user entered a file with a path then preserve it	
			}
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			string temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, cutoff); 
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; allLines = 1; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(taxonomyFileName);
                    parser.getNameFile(files);
                }
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "PhylotypeCommand");
		exit(1);
	}
}
/**********************************************************************************************************************/

int PhylotypeCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//reads in taxonomy file and makes all the taxonomies the same length 
		//by appending the last taxon to a given taxonomy as many times as needed to 
		//make it as long as the longest taxonomy in the file 
		TaxEqualizer* taxEqual = new TaxEqualizer(taxonomyFileName, cutoff, outputDir);
		
		if (m->control_pressed) { delete taxEqual; return 0; }
		
		string equalizedTaxFile = taxEqual->getEqualizedTaxFile();
		
		delete taxEqual;
		
		//build taxonomy tree from equalized file
		PhyloTree* tree = new PhyloTree(equalizedTaxFile);
		vector<int> leaves = tree->getGenusNodes();
		
		//store leaf nodes in current map
		for (int i = 0; i < leaves.size(); i++)		{	currentNodes[leaves[i]] = leaves[i];	}
		
		bool done = false;
		if (tree->get(leaves[0]).parent == -1) {  m->mothurOut("Empty Tree"); m->mothurOutEndLine();	done = true;	}
		
		if (m->control_pressed) { delete tree; return 0; }
		
		ofstream outList, outRabund, outSabund;
        map<string, string> variables;
        string fileroot = outputDir + m->getRootName(m->getSimpleName(taxonomyFileName));
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = "tx";
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        
        map<string, int> counts;
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	outSabund);
            m->openOutputFile(rabundFileName,	outRabund);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
        }else {
            CountTable ct;
            ct.readTable(countfile, false, false);
            counts = ct.getNameMap();
        }
		m->openOutputFile(listFileName,	outList);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);

        
		int count = 1;		
		//start at leaves of tree and work towards root, processing the labels the user wants
		while((!done) && ((allLines == 1) || (labels.size() != 0))) {
		
			string level = toString(count); 
			count++;
			
			if (m->control_pressed) { 
				if (countfile == "") { outRabund.close(); outSabund.close(); } outList.close();
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
				delete tree; return 0; 
			}
			
			//is this a level the user want output for
			if(allLines == 1 || labels.count(level) == 1){	
				
				//output level
				m->mothurOut(level); m->mothurOutEndLine();
				
				ListVector list;
				list.setLabel(level);
                
				//go through nodes and build listvector 
				for (itCurrent = currentNodes.begin(); itCurrent != currentNodes.end(); itCurrent++) {
			
					//get parents
					TaxNode node = tree->get(itCurrent->first);
					parentNodes[node.parent] = node.parent;
					
					vector<string> names = node.accessions;
					
					//make the names compatable with listvector
					string name = "";
					for (int i = 0; i < names.size(); i++) {  
                        
                        if (names[i] != "unknown") {
                            if (namefile != "") {	
                                map<string, string>::iterator itNames = namemap.find(names[i]);  //make sure this name is in namefile
                                
                                if (itNames != namemap.end()) {  name += namemap[names[i]] + ",";   } //you found it in namefile
                                else { m->mothurOut("[ERROR]: " + names[i] + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true;  }
                                
                            }else{   name += names[i] + ",";	}
                        }
					}
                    
                    if (m->control_pressed) { break; }
                    
					name = name.substr(0, name.length()-1);  //rip off extra ','
					//add bin to list vector
					if (name != "") { list.push_back(name); } //caused by unknown
				}	
				
				//print listvector
                if (!m->printedListHeaders) { list.printHeaders(outList); }
                if (countfile == "") { list.print(outList);  }
                else { list.print(outList, counts);  }
                
                if (countfile == "") {
                    //print rabund
                    list.getRAbundVector().print(outRabund);
                    //print sabund
                    list.getSAbundVector().print(outSabund);
                }
				labels.erase(level);
				
			}else {
				
				//just get parents
				for (itCurrent = currentNodes.begin(); itCurrent != currentNodes.end(); itCurrent++) {
					int parent = tree->get(itCurrent->first).parent;
					parentNodes[parent] = parent;
				}
			}
			
			//move up a level
			currentNodes = parentNodes;
			parentNodes.clear();
			
			//have we reached the rootnode
			if (tree->get(currentNodes.begin()->first).parent == -1)  {  done = true;  }
		}
			
		outList.close();
        if (countfile == "") {
            outSabund.close();
            outRabund.close();
        }
		
		delete tree;
		
		if (m->control_pressed) { 
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }
			return 0; 
		}
		
		//set list file as new current listfile
		string current = "";
		itTypes = outputTypes.find("list");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
		}
		
		//set rabund file as new current rabundfile
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		//set sabund file as new current sabundfile
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSabundFile(current); }
		}		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "execute");
		exit(1);
	}
}
/*****************************************************************/
int PhylotypeCommand::readNamesFile() {
	try {
				
		ifstream in;
		m->openInputFile(namefile, in);
		
		string first, second;
		map<string, string>::iterator itNames;
		
		while(!in.eof()) {
			in >> first >> second; m->gobble(in);
			
			itNames = namemap.find(first);
			if (itNames == namemap.end()) {  
				namemap[first] = second; 
			}else {  m->mothurOut(first + " has already been seen in namefile, disregarding names file."); m->mothurOutEndLine(); in.close(); namemap.clear(); namefile = ""; return 1; }			
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "readNamesFile");
		exit(1);
	}
}

/**********************************************************************************************************************/
