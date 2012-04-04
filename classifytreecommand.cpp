//
//  classifytreecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 2/20/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "classifytreecommand.h"
#include "phylotree.h"

//**********************************************************************************************************************
vector<string> ClassifyTreeCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "", "", "none",false,true); parameters.push_back(ptree);
        CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "", "", "none",false,true); parameters.push_back(ptaxonomy);
        CommandParameter pname("name", "InputTypes", "", "", "", "", "none",false,false); parameters.push_back(pname);
        CommandParameter pgroup("group", "InputTypes", "", "", "", "", "none",false,false); parameters.push_back(pgroup);
        CommandParameter pcutoff("cutoff", "Number", "", "51", "", "", "",false,true); parameters.push_back(pcutoff);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifyTreeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The classify.tree command reads a tree and taxonomy file and output the consensus taxonomy for each node on the tree. \n";
		helpString += "If you provide a group file, the concensus for each group will also be provided. \n";
		helpString += "The new tree contains labels at each internal node.  The label is the node number so you can relate the tree to the summary file.\n";
		helpString += "The summary file lists the concensus taxonomy for the descendants of each node.\n";
		helpString += "The classify.tree command parameters are tree, group, name and taxonomy. The tree and taxonomy files are required.\n";
        helpString += "The cutoff parameter allows you to specify a consensus confidence threshold for your taxonomy.  The default is 51, meaning 51%. Cutoff cannot be below 51.\n";
        helpString += "The classify.tree command should be used in the following format: classify.tree(tree=test.tre, group=test.group, taxonomy=test.taxonomy)\n";
		helpString += "Note: No spaces between parameter labels (i.e. tree), '=' and parameters (i.e.yourTreefile).\n"; 
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
ClassifyTreeCommand::ClassifyTreeCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "ClassifyTreeCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ClassifyTreeCommand::ClassifyTreeCommand(string option)  {
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
			
			m->runParse = true;
			m->clearGroups();
			m->clearAllGroups();
			m->Treenames.clear();
			m->names.clear();
			
			vector<string> tempOutNames;
			outputTypes["tree"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
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
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
			}
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
            
			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { treefile = ""; abort = true; }
			else if (treefile == "not found") { treefile = ""; 
                treefile = m->getTreeFile(); 
                if (treefile != "") {  m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
                else { m->mothurOut("No valid current files. You must provide a tree file."); m->mothurOutEndLine(); abort = true; }
            }else { m->setTreeFile(treefile); }	
            
            taxonomyfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyfile == "not open") { taxonomyfile = ""; abort = true; }
			else if (taxonomyfile == "not found") { taxonomyfile = ""; 
                taxonomyfile = m->getTaxonomyFile(); 
                if (taxonomyfile != "") {  m->mothurOut("Using " + taxonomyfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
                else { m->mothurOut("No valid current files. You must provide a taxonomy file."); m->mothurOutEndLine(); abort = true; }
            }else { m->setTaxonomyFile(taxonomyfile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            
            string temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "51"; }
			m->mothurConvert(temp, cutoff); 
			
			if ((cutoff < 51) || (cutoff > 100)) { m->mothurOut("cutoff must be above 50, and no greater than 100."); m->mothurOutEndLine(); abort = true;  }
            
            if (namefile == "") {
				vector<string> files; files.push_back(treefile);
				parser.getNameFile(files);
			}
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "ClassifyTreeCommand");		
		exit(1);
	}
}
//**********************************************************************************************************************

int ClassifyTreeCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		int start = time(NULL);
        
		/***************************************************/
		//    reading tree info							   //
		/***************************************************/
        m->setTreeFile(treefile);
        if (groupfile != "") {
			//read in group map info.
			tmap = new TreeMap(groupfile);
			tmap->readMap();
		}else{ //fake out by putting everyone in one group
			Tree* tree = new Tree(treefile); delete tree;  //extracts names from tree to make faked out groupmap
			tmap = new TreeMap();
			
			for (int i = 0; i < m->Treenames.size(); i++) { tmap->addSeq(m->Treenames[i], "Group1"); }
		}
		
		if (namefile != "") { readNamesFile(); }
		
		read = new ReadNewickTree(treefile);
		int readOk = read->read(tmap); 
		
		if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine(); delete tmap; delete read; return 0; }
		
		read->AssembleTrees();
		vector<Tree*> T = read->getTrees();
        Tree* outputTree = T[0]; 
		delete read;
		
		//make sure all files match
		//if you provide a namefile we will use the numNames in the namefile as long as the number of unique match the tree names size.
		int numNamesInTree;
		if (namefile != "")  {  
			if (numUniquesInName == m->Treenames.size()) {  numNamesInTree = nameMap.size();  }
			else {   numNamesInTree = m->Treenames.size();  }
		}else {  numNamesInTree = m->Treenames.size();  }
		
		
		//output any names that are in group file but not in tree
		if (numNamesInTree < tmap->getNumSeqs()) {
			for (int i = 0; i < tmap->namesOfSeqs.size(); i++) {
				//is that name in the tree?
				int count = 0;
				for (int j = 0; j < m->Treenames.size(); j++) {
					if (tmap->namesOfSeqs[i] == m->Treenames[j]) { break; } //found it
					count++;
				}
				
				if (m->control_pressed) { 
					delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
					m->clearGroups();
					return 0;
				}
				
				//then you did not find it so report it 
				if (count == m->Treenames.size()) { 
					//if it is in your namefile then don't remove
					map<string, string>::iterator it = nameMap.find(tmap->namesOfSeqs[i]);
					
					if (it == nameMap.end()) {
						m->mothurOut(tmap->namesOfSeqs[i] + " is in your groupfile and not in your tree. It will be disregarded."); m->mothurOutEndLine();
						tmap->removeSeq(tmap->namesOfSeqs[i]);
						i--; //need this because removeSeq removes name from namesOfSeqs
					}
				}
			}
		}
                        
        if (m->control_pressed) { delete outputTree; delete tmap;  return 0; }
		
        readTaxonomyFile();
        
        
        /***************************************************/
        //		get concensus taxonomies                    //
        /***************************************************/
        getClassifications(outputTree);
        delete outputTree; delete tmap;
			
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} return 0; }
		
		//set tree file as new current treefile
		if (treefile != "") {
			string current = "";
			itTypes = outputTypes.find("tree");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTreeFile(current); }
			}
		}
		
		m->mothurOutEndLine(); m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to find the concensus taxonomies."); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "execute");	
		exit(1);
	}
}
//**********************************************************************************************************************
//traverse tree finding concensus taxonomy at each node
//label node with a number to relate to output summary file
//report all concensus taxonomies to file 
int ClassifyTreeCommand::getClassifications(Tree*& T){
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(treefile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(treefile)) + "taxonomy.summary";
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//print headings
		out << "TreeNode\t";
		if (groupfile != "") { out << "Group\t"; } 
        out << "NumRep\tTaxonomy" << endl; 
		
		string treeOutputDir = outputDir;
		if (outputDir == "") {  treeOutputDir += m->hasPath(treefile);  }
		string outputTreeFileName = treeOutputDir + m->getRootName(m->getSimpleName(treefile)) + "taxonomy.tre";
		
		//create a map from tree node index to names of descendants, save time later
		map<int, map<string, set<string> > > nodeToDescendants; //node# -> (groupName -> groupMembers)
		for (int i = 0; i < T->getNumNodes(); i++) {
			if (m->control_pressed) { return 0; }
			
			nodeToDescendants[i] = getDescendantList(T, i, nodeToDescendants);
		}
		
		//for each node
		for (int i = T->getNumLeaves(); i < T->getNumNodes(); i++) {
			
			if (m->control_pressed) { out.close(); return 0; }
            
			string tax = "not classifed";
            int size;
            if (groupfile != "") {
                for (map<string, set<string> >::iterator itGroups = nodeToDescendants[i].begin(); itGroups != nodeToDescendants[i].end(); itGroups++) {
                    if (itGroups->first != "AllGroups") {
                        tax = getTaxonomy(itGroups->second, size);
                        out << (i+1) << '\t' << itGroups->first << '\t' << size << '\t' << tax << endl;
                    }
                }
            }else {
                string group = "AllGroups";
                tax = getTaxonomy(nodeToDescendants[i][group], size);
                out << (i+1) << '\t' << size << '\t' << tax << endl;
            }
           			
			T->tree[i].setLabel((i+1));
		}
		out.close();
        
		ofstream outTree;
		m->openOutputFile(outputTreeFileName, outTree);
		outputNames.push_back(outputTreeFileName); outputTypes["tree"].push_back(outputTreeFileName);
		T->print(outTree, "both");
		outTree.close();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "GetConcensusTaxonomies");	
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifyTreeCommand::getTaxonomy(set<string> names, int& size) {
	try{
		string conTax = "";
        size = 0;
		        
		//create a tree containing sequences from this bin
		PhyloTree* phylo = new PhyloTree();
		
		for (set<string>::iterator it = names.begin(); it != names.end(); it++) {
            
            
			//if namesfile include the names
			if (namefile != "") {
                
				//is this sequence in the name file - namemap maps seqName -> repSeqName
				map<string, string>::iterator it2 = nameMap.find(*it);
				
				if (it2 == nameMap.end()) { //this name is not in name file, skip it
					m->mothurOut((*it) + " is not in your name file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
					
					//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
					map<string, string>::iterator itTax = taxMap.find((it2->second));
                    
					if (itTax == taxMap.end()) { //this name is not in taxonomy file, skip it
                        
						if ((*it) != (it2->second)) { m->mothurOut((*it) + " is represented by " +  it2->second + " and is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine(); }
						else {  m->mothurOut((*it) + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine(); }
					}else{
						//add seq to tree
                        int num = nameCount[(*it)]; // we know its there since we found it in nameMap
						for (int i = 0; i < num; i++) {  phylo->addSeqToTree((*it)+toString(i), it2->second);  }
                        size += num;
					}
				}
				
			}else{
				//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
				map<string, string>::iterator itTax = taxMap.find((*it));
                
				if (itTax == taxMap.end()) { //this name is not in taxonomy file, skip it
					m->mothurOut((*it) + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
					//add seq to tree
					phylo->addSeqToTree((*it), itTax->second);
                    size++;
				}
			}
            
			if (m->control_pressed) { delete phylo; return conTax; }
			
		}
		
		//build tree
		phylo->assignHeirarchyIDs(0);
		
		TaxNode currentNode = phylo->get(0);
		int myLevel = 0; 	
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
                conTax += bestChild.name + "(" + toString(consensusConfidence) + ");";
				myLevel++;
			}else{ //if no, quit
				break;
			}
			
			//move down a level
			currentNode = bestChild;
		}
		
		if (myLevel != phylo->getMaxLevel()) {
			while (myLevel != phylo->getMaxLevel()) {
				conTax += "unclassified;";
				myLevel++;
			}
		}		
		if (conTax == "") {  conTax = "no_consensus;";  }
		
		delete phylo;	
        
        return conTax;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "getTaxonomy");
		exit(1);
	}
}

//**********************************************************************************************************************
map<string, set<string> > ClassifyTreeCommand::getDescendantList(Tree*& T, int i, map<int, map<string, set<string> > > descendants){
	try {
		map<string ,set<string> > names;
		
		map<string ,set<string> >::iterator it;
        map<string ,set<string> >::iterator it2;
		
		int lc = T->tree[i].getLChild();
		int rc = T->tree[i].getRChild();
		
		if (lc == -1) { //you are a leaf your only descendant is yourself
            string group = tmap->getGroup(T->tree[i].getName());
            set<string> mynames; mynames.insert(T->tree[i].getName());
            names[group] = mynames; //mygroup -> me
            names["AllGroups"] = mynames;
		}else{ //your descedants are the combination of your childrens descendants
			names = descendants[lc];
			for (it = descendants[rc].begin(); it != descendants[rc].end(); it++) {
                it2 = names.find(it->first); //do we already have this group
                if (it2 == names.end()) { //nope, so add it
                    names[it->first] = it->second;
                }else {
                    for (set<string>::iterator it3 = (it->second).begin(); it3 != (it->second).end(); it3++) {
                        names[it->first].insert(*it3);
                    }
                }
				
			}
		}
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "getDescendantList");	
		exit(1);
	}
}
//**********************************************************************************************************************
int ClassifyTreeCommand::readTaxonomyFile() {
	try {
		
		ifstream in;
		m->openInputFile(taxonomyfile, in);
		
		string name, tax;
        
		while(!in.eof()){
			in >> name >> tax;		
			m->gobble(in);
			
			//are there confidence scores, if so remove them
			if (tax.find_first_of('(') != -1) {  m->removeConfidences(tax);	}
			
			taxMap[name] = tax;
			
			if (m->control_pressed) { in.close(); taxMap.clear(); return 0; }
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "readTaxonomyFile");
		exit(1);
	}
}

/*****************************************************************/
int ClassifyTreeCommand::readNamesFile() {
	try {
		ifstream inNames;
		m->openInputFile(namefile, inNames);
		
		string name, names;
        
		while(!inNames.eof()){
			inNames >> name;			//read from first column  A
			inNames >> names;		//read from second column  A,B,C,D
			m->gobble(inNames);
			
			//parse names into vector
			vector<string> theseNames;
			m->splitAtComma(names, theseNames);
            
			for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = name;  }
            nameCount[name] = theseNames.size();
			
			if (m->control_pressed) { inNames.close(); nameMap.clear(); return 0; }
		}
		inNames.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "readNamesFile");
		exit(1);
	}
}

/*****************************************************************/


