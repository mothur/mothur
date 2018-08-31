//
//  classifytreecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 2/20/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "classifytreecommand.h"
#include "phylotree.h"
#include "treereader.h"

//**********************************************************************************************************************
vector<string> ClassifyTreeCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "", "", "none","tree-summary",false,true,true); parameters.push_back(ptree);
        CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "", "", "none","",false,true,true); parameters.push_back(ptaxonomy);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter pmethod("output", "Multiple", "node-taxon", "node", "", "", "","",false,false); parameters.push_back(pmethod);
        CommandParameter pcutoff("cutoff", "Number", "", "51", "", "", "","",false,true); parameters.push_back(pcutoff);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
        helpString += "The count parameter allows you add a count file so you can have the summary totals broken up by group.\n";
		helpString += "The summary file lists the concensus taxonomy for the descendants of each node.\n";
		helpString += "The classify.tree command parameters are tree, group, name, count and taxonomy. The tree and taxonomy files are required.\n";
        helpString += "The cutoff parameter allows you to specify a consensus confidence threshold for your taxonomy.  The default is 51, meaning 51%. Cutoff cannot be below 51.\n";
        helpString += "The output parameter allows you to specify whether you want the tree node number displayed on the tree, or the taxonomy displayed. Default=node. Options are node or taxon.\n";
        helpString += "The classify.tree command should be used in the following format: classify.tree(tree=test.tre, group=test.group, taxonomy=test.taxonomy)\n";
		 
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyTreeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifyTreeCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],taxonomy.summary"; } //makes file like: amazon.0.03.fasta
        else if (type == "tree") {  pattern = "[filename],taxonomy.tre"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClassifyTreeCommand", "getOutputPattern");
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
			
			vector<string> tempOutNames;
			outputTypes["tree"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
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
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
            
			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { treefile = ""; abort = true; }
			else if (treefile == "not found") { treefile = ""; 
                treefile = current->getTreeFile(); 
                if (treefile != "") {  m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
                else { m->mothurOut("No valid current files. You must provide a tree file."); m->mothurOutEndLine(); abort = true; }
            }else { current->setTreeFile(treefile); }	
            
            taxonomyfile = validParameter.validFile(parameters, "taxonomy");
			if (taxonomyfile == "not open") { taxonomyfile = ""; abort = true; }
			else if (taxonomyfile == "not found") { taxonomyfile = ""; 
                taxonomyfile = current->getTaxonomyFile(); 
                if (taxonomyfile != "") {  m->mothurOut("Using " + taxonomyfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
                else { m->mothurOut("No valid current files. You must provide a taxonomy file."); m->mothurOutEndLine(); abort = true; }
            }else { current->setTaxonomyFile(taxonomyfile); }	
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
            
            string temp = validParameter.valid(parameters, "cutoff");			if (temp == "not found") { temp = "51"; }
			util.mothurConvert(temp, cutoff); 
			
			if ((cutoff < 51) || (cutoff > 100)) { m->mothurOut("cutoff must be above 50, and no greater than 100."); m->mothurOutEndLine(); abort = true;  }
            
            output = validParameter.valid(parameters, "output");
            if (output == "not found") { output = "node"; }
            
            if ((output == "node") || (output == "taxon")) {
            }else { m->mothurOut("[ERROR]: " + output + "is not a valid output option.  Valid output options are node or taxon.\n");  abort = true; }
            
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(treefile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		long start = time(NULL);
        
		/***************************************************/
		//    reading tree info							   //
		/***************************************************/
        current->setTreeFile(treefile);
        
        TreeReader* reader = new TreeReader(treefile, groupfile, namefile);
        vector<Tree*> T = reader->getTrees();
        CountTable* tmap = T[0]->getCountTable();
        Tree* outputTree = T[0];
        delete reader;

        if (namefile != "") { util.readNames(namefile, nameMap, nameCount); }
                        
        if (m->getControl_pressed()) { delete tmap;  delete outputTree;  return 0; }
		
        util.readTax(taxonomyfile, taxMap, true);
        
        /***************************************************/
        //		get concensus taxonomies                    //
        /***************************************************/
        getClassifications(outputTree);
        delete outputTree; delete tmap;
			
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0; }
		
		//set tree file as new current treefile
		if (treefile != "") {
			string currentName = "";
			itTypes = outputTypes.find("tree");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTreeFile(currentName); }
			}
		}
		
		m->mothurOutEndLine(); m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to find the concensus taxonomies."); m->mothurOutEndLine();
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(treefile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(treefile));
		string outputFileName = getOutputFileName("summary", variables);
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//print headings
		out << "TreeNode\t";
		if (groupfile != "") { out << "Group\t"; } 
        out << "NumRep\tTaxonomy" << endl; 
		
		string treeOutputDir = outputDir;
		if (outputDir == "") {  treeOutputDir += util.hasPath(treefile);  }
        variables["[filename]"] = treeOutputDir + util.getRootName(util.getSimpleName(treefile));
		string outputTreeFileName = getOutputFileName("tree", variables);
		
		//create a map from tree node index to names of descendants, save time later
		map<int, map<string, set<string> > > nodeToDescendants; //node# -> (groupName -> groupMembers)
		for (int i = 0; i < T->getNumNodes(); i++) {
            if (m->getControl_pressed()) { return 0; }
			
			nodeToDescendants[i] = getDescendantList(T, i, nodeToDescendants);
		}
		
		//for each node
		for (int i = T->getNumLeaves(); i < T->getNumNodes(); i++) {
			
			if (m->getControl_pressed()) { out.close(); return 0; }
            
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
           	
            if (output == "node") {  T->tree[i].setLabel(toString(i+1));  }
            else {
                string cleanedTax = tax;
                util.removeConfidences(cleanedTax);
                for (int j = 0; j < cleanedTax.length(); j++) {
                    //special chars to trees -  , ) ( ; [ ] :
                    if ((cleanedTax[j] == ',') || (cleanedTax[j] == '(') || (cleanedTax[j] == ')') || (cleanedTax[j] == ';') || (cleanedTax[j] == ':') || (cleanedTax[j] == ']') || (cleanedTax[j] == '[')) {
                        cleanedTax[j] = '_'; //change any special chars to _ so the tree can be read by tree readers
                    }
                }
                cout << tax << '\t' << cleanedTax << endl;
                T->tree[i].setLabel(cleanedTax);
            }
            
		}
		out.close();
        
		ofstream outTree;
		util.openOutputFile(outputTreeFileName, outTree);
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
						for (int i = 0; i < num; i++) {  phylo->addSeqToTree((*it)+toString(i), itTax->second);  }
                        size += num;
					}
				}
				
			}else{
				//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
				map<string, string>::iterator itTax = taxMap.find((*it));
                
				if (itTax == taxMap.end()) { //this name is not in taxonomy file, skip it
					m->mothurOut((*it) + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
					if (countfile != "") {
                        int numDups = ct->getNumSeqs((*it)); 
                        for (int j = 0; j < numDups; j++) {  phylo->addSeqToTree((*it), itTax->second);  }
                        size += numDups;
                    }else{
                        //add seq to tree
                        phylo->addSeqToTree((*it), itTax->second);
                        size++;  
                    }				}
			}
            
			if (m->getControl_pressed()) { delete phylo; return conTax; }
			
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
       // TreeMap* tmap = T->getTreeMap();
		
		if (lc == -1) { //you are a leaf your only descendant is yourself
            vector<string> groups = T->tree[i].getGroup();
            set<string> mynames; mynames.insert(T->tree[i].getName());
            for (int j = 0; j < groups.size(); j++) { names[groups[j]] = mynames;   } //mygroup -> me
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
/*****************************************************************/


