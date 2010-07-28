/*
 *  classifyotucommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "classifyotucommand.h"
#include "phylotree.h"


//**********************************************************************************************************************
ClassifyOtuCommand::ClassifyOtuCommand(string option)  {
	try{
		abort = false;
		allLines = 1;
		labels.clear();
				
		//allow user to run help
		if (option == "help") { 
			help(); abort = true;
		} else {
			//valid paramters for this command
			string Array[] =  {"list","label","name","taxonomy","cutoff","probs","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
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
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") { m->mothurOut("list is a required parameter for the classify.otu command."); m->mothurOutEndLine(); abort = true; }
			else if (listfile == "not open") { abort = true; }	
			
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not found") {  m->mothurOut("taxonomy is a required parameter for the classify.otu command."); m->mothurOutEndLine(); abort = true; }
			else if (taxfile == "not open") { abort = true; }	
	
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; allLines = 1;  }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			string temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "51"; }
			convert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "probs", false);					if (temp == "not found"){	temp = "true";			}
			probs = isTrue(temp);
			
			
			if ((cutoff < 51) || (cutoff > 100)) { m->mothurOut("cutoff must be above 50, and no greater than 100."); m->mothurOutEndLine(); abort = true;  }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "ClassifyOtuCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void ClassifyOtuCommand::help(){
	try {
		m->mothurOut("The classify.otu command parameters are list, taxonomy, name, cutoff, label and probs.  The taxonomy and list parameters are required.\n");
		m->mothurOut("The name parameter allows you add a names file with your taxonomy file.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like a output files created for, and is separated by dashes.\n");
		m->mothurOut("The default value for label is all labels in your inputfile.\n");
		m->mothurOut("The cutoff parameter allows you to specify a consensus confidence threshold for your taxonomy.  The default is 51, meaning 51%. Cutoff cannot be below 51.\n");
		m->mothurOut("The probs parameter shuts off the outputting of the consensus confidence results. The default is true, meaning you want the confidence to be shown.\n");
		m->mothurOut("The classify.otu command should be in the following format: classify.otu(taxonomy=yourTaxonomyFile, list=yourListFile, name=yourNamesFile, label=yourLabels).\n");
		m->mothurOut("Example classify.otu(taxonomy=abrecovery.silva.full.taxonomy, list=abrecovery.fn.list, label=0.10).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

ClassifyOtuCommand::~ClassifyOtuCommand(){}

//**********************************************************************************************************************

int ClassifyOtuCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		//if user gave a namesfile then use it
		if (namefile != "") {	readNamesFile();	}
		
		//read taxonomy file and save in map for easy access in building bin trees
		readTaxonomyFile();
		
		if (m->control_pressed) { return 0; }
		
		input = new InputData(listfile, "list");
		list = input->getListVector();
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }  return 0; }
	
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (allLines == 1 || labels.count(list->getLabel()) == 1){
			
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					process(list);
					if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } delete input; delete list; return 0; }
										
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list;
					list = input->getListVector(lastLabel);
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					process(list);
				
					
					if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } delete input; delete list; return 0; }
										
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();
	
			delete list;
			list = input->getListVector();
		}
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + (*it)); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
			m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
			
			process(list);
			delete list;
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } delete input; delete list; return 0; }
		}
		
		delete input;  
				
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int ClassifyOtuCommand::readNamesFile() {
	try {
		
		ifstream inNames;
		openInputFile(namefile, inNames);
		
		string name, names;
	
		while(inNames){
			inNames >> name;			//read from first column  A
			inNames >> names;		//read from second column  A,B,C,D
			gobble(inNames);
			
			//parse names into vector
			vector<string> theseNames;
			splitAtComma(names, theseNames);

			for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = name;  }
			
			if (m->control_pressed) { inNames.close(); nameMap.clear(); return 0; }
		}
		inNames.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "readNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClassifyOtuCommand::readTaxonomyFile() {
	try {
		
		ifstream in;
		openInputFile(taxfile, in);
		
		string name, tax;
	
		while(!in.eof()){
			in >> name >> tax;		
			gobble(in);
			
			//are there confidence scores, if so remove them
			if (tax.find_first_of('(') != -1) {  removeConfidences(tax);	}
			
			taxMap[name] = tax;
			
			if (m->control_pressed) { in.close(); taxMap.clear(); return 0; }
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "readTaxonomyFile");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifyOtuCommand::findConsensusTaxonomy(int bin, ListVector* thisList, int& size) {
	try{
		string conTax = "";
		vector<string> names;
		map<string, string>::iterator it;
		map<string, string>::iterator it2;

		//parse names into vector
		string binnames = thisList->get(bin);
		splitAtComma(binnames, names);

		//create a tree containing sequences from this bin
		PhyloTree* phylo = new PhyloTree();
		
		size = 0;
		for (int i = 0; i < names.size(); i++) {
		
			//if namesfile include the names
			if (namefile != "") {
				//is this sequence in the name file - namemap maps seqName -> repSeqName
				it2 = nameMap.find(names[i]);
				
				if (it2 == nameMap.end()) { //this name is not in name file, skip it
					m->mothurOut(names[i] + " is not in your name file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
					
					//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
					it = taxMap.find(it2->second);
			
					if (it == taxMap.end()) { //this name is not in taxonomy file, skip it
						if (names[i] != it->second) { m->mothurOut(names[i] + " is represented by " +  it2->second + " and is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine(); }
						else {  m->mothurOut(names[i] + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine(); }
					}else{
				
						//add seq to tree
						phylo->addSeqToTree(names[i], it->second);
						size++;
					}
				}
				
			}else{
				//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
				it = taxMap.find(names[i]);
		
				if (it == taxMap.end()) { //this name is not in taxonomy file, skip it
					m->mothurOut(names[i] + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
					//add seq to tree
					phylo->addSeqToTree(names[i], it->second);
					size++;
				}
			}

			
			if (m->control_pressed) { delete phylo; return conTax; }
			
		}
		
		//build tree
		phylo->assignHeirarchyIDs(0);
		
		TaxNode currentNode = phylo->get(0);
		
		//at each level
		while (currentNode.children.size() != 0) { //you still have more to explore
		
			TaxNode bestChild;
			int bestChildSize = 0;
			
			//go through children
			for (map<string, int>::iterator itChild = currentNode.children.begin(); itChild != currentNode.children.end(); itChild++) {
				
				TaxNode temp = phylo->get(itChild->second);
				
				//select child with largest accesions - most seqs assigned to it
				if (temp.accessions.size() > bestChildSize) {
					bestChild = phylo->get(itChild->second);
					bestChildSize = temp.accessions.size();
				}
				
			}
				
			//is this taxonomy above cutoff
			int consensusConfidence = ceil((bestChildSize / (float) size) * 100);
			
			if (consensusConfidence >= cutoff) { //if yes, add it
				if (probs) {
					conTax += bestChild.name + "(" + toString(consensusConfidence) + ");";
				}else{
					conTax += bestChild.name + ";";
				}
			}else{ //if no, quit
				break;
			}
			
			//move down a level
			currentNode = bestChild;
		}
		
				
		if (conTax == "") {  conTax = "unclassified;";  }
		
		delete phylo;	
		
		return conTax;
			
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "findConsensusTaxonomy");
		exit(1);
	}
}

//**********************************************************************************************************************
int ClassifyOtuCommand::process(ListVector* processList) {
	try{
		string conTax;
		int size;
		
		//create output file
		if (outputDir == "") { outputDir += hasPath(listfile); }
				
		ofstream out;
		string outputFile = outputDir + getRootName(getSimpleName(listfile)) + processList->getLabel() + ".cons.taxonomy";
		openOutputFile(outputFile, out);
		outputNames.push_back(outputFile);
		
		//for each bin in the list vector
		for (int i = 0; i < processList->getNumBins(); i++) {
			conTax  = findConsensusTaxonomy(i, processList, size);
			
			if (m->control_pressed) { out.close();  return 0; }
			
			//output to new names file
			out << (i+1) << '\t' << size << '\t' << conTax << endl;
		}

		out.close();
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
void ClassifyOtuCommand::removeConfidences(string& tax) {
	try {
		
		string taxon;
		string newTax = "";
		
		while (tax.find_first_of(';') != -1) {
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'));
			
			int pos = taxon.find_first_of('(');
			if (pos != -1) {
				taxon = taxon.substr(0, pos); //rip off confidence 
			}
			
			taxon += ";";
			
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			newTax += taxon;
		}
		
		tax = newTax;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "removeConfidences");
		exit(1);
	}
}
//**********************************************************************************************************************


