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

//**********************************************************************************************************************
vector<string> PhylotypeCommand::getValidParameters(){	
	try {
		string Array[] =  {"taxonomy","cutoff","label","name","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
PhylotypeCommand::PhylotypeCommand(){	
	try {
		abort = true; calledHelp = true; 
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
//**********************************************************************************************************************
vector<string> PhylotypeCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"taxonomy"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PhylotypeCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "getRequiredFiles");
		exit(1);
	}
}
/**********************************************************************************************************************/
PhylotypeCommand::PhylotypeCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"taxonomy","cutoff","label","name","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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
			}

			taxonomyFileName = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyFileName == "not found") { 
				m->mothurOut("taxonomy is a required parameter for the phylotype command."); 
				m->mothurOutEndLine();
				abort = true; 
			}else if (taxonomyFileName == "not open") { abort = true; }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { readNamesFile(); }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(taxonomyFileName); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "-1"; }
			convert(temp, cutoff); 
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; allLines = 1; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "PhylotypeCommand");
		exit(1);
	}
}
/**********************************************************************************************************************/

void PhylotypeCommand::help(){
	try {
		m->mothurOut("The phylotype command reads a taxonomy file and outputs a .list, .rabund and .sabund file. \n");
		m->mothurOut("The phylotype command parameter options are taxonomy, cutoff and label. The taxonomy parameter is required.\n");
		m->mothurOut("The cutoff parameter allows you to specify the level you want to stop at.  The default is the highest level in your taxonomy file. \n");
		m->mothurOut("For example: taxonomy = Bacteria;Bacteroidetes-Chlorobi;Bacteroidetes; - cutoff=2, would truncate the taxonomy to Bacteria;Bacteroidetes-Chlorobi; \n");
		m->mothurOut("For the cutoff parameter levels count up from the root of the phylotree. This enables you to look at the grouping down to a specific resolution, say the genus level.\n");
		m->mothurOut("The label parameter allows you to specify which level you would like, and are separated by dashes.  The default all levels in your taxonomy file. \n");
		m->mothurOut("For the label parameter, levels count down from the root to keep the output similiar to mothur's other commands which report information from finer resolution to coarser resolutions.\n");
		m->mothurOut("The phylotype command should be in the following format: \n");
		m->mothurOut("phylotype(taxonomy=yourTaxonomyFile, cutoff=yourCutoff, label=yourLabels) \n");
		m->mothurOut("Eaxample: phylotype(taxonomy=amazon.taxonomy, cutoff=5, label=1-3-5).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "PhylotypeCommand", "help");
		exit(1);
	}
}
/**********************************************************************************************************************/

PhylotypeCommand::~PhylotypeCommand(){}

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
		
		string fileroot = outputDir + m->getRootName(m->getSimpleName(taxonomyFileName));
		
		ofstream outList;
		string outputListFile = fileroot + "tx.list";
		m->openOutputFile(outputListFile, outList);
		ofstream outSabund;
		string outputSabundFile = fileroot + "tx.sabund";
		m->openOutputFile(outputSabundFile, outSabund);
		ofstream outRabund;
		string outputRabundFile = fileroot + "tx.rabund";
		m->openOutputFile(outputRabundFile, outRabund);
		
		outputNames.push_back(outputListFile); outputTypes["list"].push_back(outputListFile);
		outputNames.push_back(outputSabundFile); outputTypes["sabund"].push_back(outputSabundFile);
		outputNames.push_back(outputRabundFile); outputTypes["rabund"].push_back(outputRabundFile);
		
		int count = 1;		
		//start at leaves of tree and work towards root, processing the labels the user wants
		while((!done) && ((allLines == 1) || (labels.size() != 0))) {
		
			string level = toString(count); 
			count++;
			
			if (m->control_pressed) { 
				outRabund.close(); outSabund.close(); outList.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
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
						if (namefile != "") {	
							map<string, string>::iterator itNames = namemap.find(names[i]);  //make sure this name is in namefile
		
							if (itNames != namemap.end()) {  name += namemap[names[i]] + ",";   } //you found it in namefile
							else { m->mothurOut(names[i] + " is not in your namefile, please correct."); m->mothurOutEndLine(); exit(1);  }
							
						}else{   name += names[i] + ",";	}
					}
					name = name.substr(0, name.length()-1);  //rip off extra ','
					
					//add bin to list vector
					list.push_back(name);
				}	
				
				//print listvector
				list.print(outList);
				//print rabund
				list.getRAbundVector().print(outRabund);
				//print sabund
				list.getSAbundVector().print(outSabund);
			
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
		outSabund.close();
		outRabund.close();	
		
		delete tree;
		
		if (m->control_pressed) { 
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
			return 0; 
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
