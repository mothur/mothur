/*
 *  treegroupscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treegroupscommand.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"


//**********************************************************************************************************************

TreeGroupCommand::TreeGroupCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		lines.clear();
		labels.clear();
		Groups.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { validCalculator = new ValidCalculators(); help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"line","label","calc","groups", "phylip", "column", "name", "precision","cutoff"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  format = "phylip"; 	}
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  format = "column";	}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else {  globaldata->setNameFile(namefile);	}
			
			format = globaldata->getFormat();
			
			//error checking on files			
			if ((globaldata->getSharedFile() == "") && ((phylipfile == "") && (columnfile == "")))	{ mothurOut("You must run the read.otu command or provide a distance file before running the tree.shared command."); mothurOutEndLine(); abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { mothurOut("When running the tree.shared command with a distance file you may not use both the column and the phylip parameters."); mothurOutEndLine(); abort = true; }
			
			if (columnfile != "") {
				if (namefile == "") {  mothurOut("You need to provide a namefile if you are going to use the column format."); mothurOutEndLine(); abort = true; }
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter.validFile(parameters, "line", false);				
			if (line == "not found") { line = "";  }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//make sure user did not use both the line and label parameters
			if ((line != "") && (label != "")) { mothurOut("You cannot use both the line and label parameters at the same time. "); mothurOutEndLine(); abort = true; }
			//if the user has not specified any line or labels use the ones from read.otu
			else if((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
				
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "precision", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));

				
			if (abort == false) {
			
				validCalculator = new ValidCalculators();
				
				if (format == "sharedfile") {
					int i;
					for (i=0; i<Estimators.size(); i++) {
						if (validCalculator->isValidCalculator("treegroup", Estimators[i]) == true) { 
							if (Estimators[i] == "jabund") { 	
								treeCalculators.push_back(new JAbund());
							}else if (Estimators[i] == "sorabund") { 
								treeCalculators.push_back(new SorAbund());
							}else if (Estimators[i] == "jclass") { 
								treeCalculators.push_back(new Jclass());
							}else if (Estimators[i] == "sorclass") { 
								treeCalculators.push_back(new SorClass());
							}else if (Estimators[i] == "jest") { 
								treeCalculators.push_back(new Jest());
							}else if (Estimators[i] == "sorest") { 
								treeCalculators.push_back(new SorEst());
							}else if (Estimators[i] == "thetayc") { 
								treeCalculators.push_back(new ThetaYC());
							}else if (Estimators[i] == "thetan") { 
								treeCalculators.push_back(new ThetaN());
							}else if (Estimators[i] == "morisitahorn") { 
								treeCalculators.push_back(new MorHorn());
							}else if (Estimators[i] == "braycurtis") { 
								treeCalculators.push_back(new BrayCurtis());
							}
						}
					}
				}
			}	
		}

	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "TreeGroupCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void TreeGroupCommand::help(){
	try {
		mothurOut("The tree.shared command creates a .tre to represent the similiarity between groups or sequences.\n");
		mothurOut("The tree.shared command can only be executed after a successful read.otu command or by providing a distance file.\n");
		mothurOut("The tree.shared command parameters are groups, calc, phylip, column, name, cutoff, precision, line and label.  You may not use line and label at the same time.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n");
		mothurOut("The group names are separated by dashes. The line and label allow you to select what distance levels you would like trees created for, and are also separated by dashes.\n");
		mothurOut("The phylip or column parameter are required if you do not run the read.otu command first, and only one may be used.  If you use a column file the name filename is required. \n");
		mothurOut("If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.\n");
		mothurOut("The tree.shared command should be in the following format: tree.shared(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels).\n");
		mothurOut("Example tree.shared(groups=A-B-C, line=1-3-5, calc=jabund-sorabund).\n");
		mothurOut("The default value for groups is all the groups in your groupfile.\n");
		mothurOut("The default value for calc is jclass-thetayc.\n");
		mothurOut("The tree.shared command outputs a .tre file for each calculator you specify at each distance you choose.\n");
		validCalculator->printCalc("treegroup", cout);
		mothurOut("Or the tree.shared command can be in the following format: tree.shared(phylip=yourPhylipFile).\n");
		mothurOut("Example tree.shared(phylip=abrecovery.dist).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

TreeGroupCommand::~TreeGroupCommand(){
	if (abort == false) {
		
		if (format == "sharedfile") { delete read;  delete input; globaldata->ginput = NULL;}
		else { delete readMatrix;  delete matrix; delete list; }
		delete tmap;
		delete validCalculator;
	}
	
}

//**********************************************************************************************************************

int TreeGroupCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		if (format == "sharedfile") {
			//if the users entered no valid calculators don't execute command
			if (treeCalculators.size() == 0) { mothurOut("You have given no valid calculators."); mothurOutEndLine(); return 0; }

			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup[0]->getLabel();
			
			if (lookup.size() < 2) { mothurOut("You have not provided enough valid groups.  I cannot run the command."); mothurOutEndLine(); return 0; }
			
			//used in tree constructor 
			globaldata->runParse = false;
			
			//create tree file
			makeSimsShared();
		}else{
			//read in dist file
			filename = globaldata->inputFileName;
		
			if (format == "column") { readMatrix = new ReadColumnMatrix(filename); }	
			else if (format == "phylip") { readMatrix = new ReadPhylipMatrix(filename); }
				
			readMatrix->setCutoff(cutoff);
	
			if(namefile != ""){	
				nameMap = new NameAssignment(namefile);
				nameMap->readMap();
			}
			else{
				nameMap = NULL;
			}
	
			readMatrix->read(nameMap);
			list = readMatrix->getListVector();
			matrix = readMatrix->getMatrix();

			//make treemap
			tmap = new TreeMap();
			tmap->makeSim(list);
			globaldata->gTreemap = tmap;
			
			globaldata->Groups = tmap->namesOfGroups;
		
			//clear globaldatas old tree names if any
			globaldata->Treenames.clear();
		
			//fills globaldatas tree names
			globaldata->Treenames = globaldata->Groups;
			
			//used in tree constructor 
			globaldata->runParse = false;
			
			makeSimsDist();

			//create a new filename
			outputFile = getRootName(globaldata->inputFileName) + "tre";	
				
			createTree();
			mothurOut("Tree complete. "); mothurOutEndLine();
		}
				
		//reset groups parameter
		globaldata->Groups.clear();  

		return 0;
	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

void TreeGroupCommand::createTree(){
	try {
		//create tree
		t = new Tree();
		
		//do merges and create tree structure by setting parents and children
		//there are numGroups - 1 merges to do
		for (int i = 0; i < (numGroups - 1); i++) {
			float largest = -1000.0;

			int row, column;
			//find largest value in sims matrix by searching lower triangle
			for (int j = 1; j < simMatrix.size(); j++) {
				for (int k = 0; k < j; k++) {
					if (simMatrix[j][k] > largest) {  largest = simMatrix[j][k]; row = j; column = k;  }
				}
			}

			//set non-leaf node info and update leaves to know their parents
			//non-leaf
			t->tree[numGroups + i].setChildren(index[row], index[column]);
			
			//parents
			t->tree[index[row]].setParent(numGroups + i);
			t->tree[index[column]].setParent(numGroups + i);
			
			//blength = distance / 2;
			float blength = ((1.0 - largest) / 2);
			
			//branchlengths
			t->tree[index[row]].setBranchLength(blength - t->tree[index[row]].getLengthToLeaves());
			t->tree[index[column]].setBranchLength(blength - t->tree[index[column]].getLengthToLeaves());
			
			//set your length to leaves to your childs length plus branchlength
			t->tree[numGroups + i].setLengthToLeaves(t->tree[index[row]].getLengthToLeaves() + t->tree[index[row]].getBranchLength());
			
			
			//update index 
			index[row] = numGroups+i;
			index[column] = numGroups+i;
			
			//remove highest value that caused the merge.
			simMatrix[row][column] = -1000.0;
			simMatrix[column][row] = -1000.0;
			
			//merge values in simsMatrix
			for (int n = 0; n < simMatrix.size(); n++)	{
				//row becomes merge of 2 groups
				simMatrix[row][n] = (simMatrix[row][n] + simMatrix[column][n]) / 2;
				simMatrix[n][row] = simMatrix[row][n];
				//delete column
				simMatrix[column][n] = -1000.0;
				simMatrix[n][column] = -1000.0;
			}
		}
		
		//adjust tree to make sure root to tip length is .5
		int root = t->findRoot();
		t->tree[root].setBranchLength((0.5 - t->tree[root].getLengthToLeaves()));
		
		//assemble tree
		t->assembleTree();
		
		//print newick file
		t->createNewickFile(outputFile);
		
		//delete tree
		delete t;
	
	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "createTree");
		exit(1);
	}
}
/***********************************************************/
void TreeGroupCommand::printSims(ostream& out) {
	try {
		
		//output column headers
		//out << '\t';
		//for (int i = 0; i < lookup.size(); i++) {	out << lookup[i]->getGroup() << '\t';		}
		//out << endl;
		
		
		for (int m = 0; m < simMatrix.size(); m++)	{
			//out << lookup[m]->getGroup() << '\t';
			for (int n = 0; n < simMatrix.size(); n++)	{
				out << simMatrix[m][n] << '\t'; 
			}
			out << endl;
		}

	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "printSims");
		exit(1);
	}
}
/***********************************************************/
void TreeGroupCommand::makeSimsDist() {
	try {
		numGroups = list->size();
		
		//initialize index
		index.clear();
		for (int g = 0; g < numGroups; g++) {	index[g] = g;	}
		
		//initialize simMatrix
		simMatrix.clear();
		simMatrix.resize(numGroups);
		for (int m = 0; m < simMatrix.size(); m++)	{
			for (int j = 0; j < simMatrix.size(); j++)	{
				simMatrix[m].push_back(0.0);
			}
		}
		
		//go through sparse matrix and fill sims
		//go through each cell in the sparsematrix
		for(MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++){
			//similairity = -(distance-1)
			simMatrix[currentCell->row][currentCell->column] = -(currentCell->dist -1.0);	
			simMatrix[currentCell->column][currentCell->row] = -(currentCell->dist -1.0);				
		}


	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "makeSimsDist");
		exit(1);
	}
}

/***********************************************************/
void TreeGroupCommand::makeSimsShared() {
	try {
		int count = 1;	
	
		//clear globaldatas old tree names if any
		globaldata->Treenames.clear();
		
		//fills globaldatas tree names
		globaldata->Treenames = globaldata->Groups;
		
		//create treemap class from groupmap for tree class to use
		tmap = new TreeMap();
		tmap->makeSim(globaldata->gGroupmap);
		globaldata->gTreemap = tmap;
		
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;

		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
			if(allLines == 1 || lines.count(count) == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);

				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				process(lookup);
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}

			lastLabel = lookup[0]->getLabel();			
			
			//get next line to process
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
				needToRun = true;
			}else {
				mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {		delete lookup[i]; }		} 
			lookup = input->getSharedRAbundVectors(lastLabel);

			mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
			process(lookup);
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 	
		}
		
		for(int i = 0 ; i < treeCalculators.size(); i++) {  delete treeCalculators[i]; }
	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "makeSimsShared");
		exit(1);
	}
}

/***********************************************************/
void TreeGroupCommand::process(vector<SharedRAbundVector*> thisLookup) {
	try{
				EstOutput data;
				vector<SharedRAbundVector*> subset;
				numGroups = thisLookup.size();
				
				//for each calculator												
				for(int i = 0 ; i < treeCalculators.size(); i++) {
					//initialize simMatrix
					simMatrix.clear();
					simMatrix.resize(numGroups);
					for (int m = 0; m < simMatrix.size(); m++)	{
						for (int j = 0; j < simMatrix.size(); j++)	{
							simMatrix[m].push_back(0.0);
						}
					}
		
					//initialize index
					index.clear();
					for (int g = 0; g < numGroups; g++) {	index[g] = g;	}
		
					//create a new filename
					outputFile = getRootName(globaldata->inputFileName) + treeCalculators[i]->getName() + "." + thisLookup[0]->getLabel() + ".tre";				
												
					for (int k = 0; k < thisLookup.size(); k++) { 
						for (int l = k; l < thisLookup.size(); l++) {
							if (k != l) { //we dont need to similiarity of a groups to itself
								//get estimated similarity between 2 groups
								
								subset.clear(); //clear out old pair of sharedrabunds
								//add new pair of sharedrabunds
								subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
								
								data = treeCalculators[i]->getValues(subset); //saves the calculator outputs
								//save values in similarity matrix
								simMatrix[k][l] = data[0];
								simMatrix[l][k] = data[0];
							}
						}
					}
				
					//creates tree from similarity matrix and write out file
					createTree();
				}

	}
	catch(exception& e) {
		errorOut(e, "TreeGroupCommand", "process");
		exit(1);
	}
}
/***********************************************************/

	

