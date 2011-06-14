/*
 *  treegroupscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treegroupscommand.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "sharedjackknife.h"
#include "whittaker.h"
#include "odum.h"
#include "canberra.h"
#include "structeuclidean.h"
#include "structchord.h"
#include "hellinger.h"
#include "manhattan.h"
#include "structpearson.h"
#include "soergel.h"
#include "spearman.h"
#include "structkulczynski.h"
#include "structchi2.h"
#include "speciesprofile.h"
#include "hamming.h"
#include "gower.h"
#include "memchi2.h"
#include "memchord.h"
#include "memeuclidean.h"
#include "mempearson.h"

//**********************************************************************************************************************
vector<string> TreeGroupCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none",false,false); parameters.push_back(pshared);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "ColumnName",false,false); parameters.push_back(pcolumn);		
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);		
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pcalc("calc", "Multiple", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-sharednseqs-ochiai-anderberg-kulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-whittaker-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-hamming-structchi2-gower-memchi2-memchord-memeuclidean-mempearson", "jclass-thetayc", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter poutput("output", "Multiple", "lt-square", "lt", "", "", "",false,false); parameters.push_back(poutput);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string TreeGroupCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The tree.shared command creates a .tre to represent the similiarity between groups or sequences.\n";
		helpString += "The tree.shared command parameters are shared, groups, calc, phylip, column, name, cutoff, precision and label.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n";
		helpString += "The group names are separated by dashes. The label allow you to select what distance levels you would like trees created for, and are also separated by dashes.\n";
		helpString += "The phylip or column parameter are required if you do not provide a sharedfile, and only one may be used.  If you use a column file the name filename is required. \n";
		helpString += "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.\n";
		helpString += "The tree.shared command should be in the following format: tree.shared(groups=yourGroups, calc=yourCalcs, label=yourLabels).\n";
		helpString += "Example tree.shared(groups=A-B-C, calc=jabund-sorabund).\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The default value for calc is jclass-thetayc.\n";
		helpString += "The tree.shared command outputs a .tre file for each calculator you specify at each distance you choose.\n";
		helpString += validCalculator.printCalc("treegroup");
		helpString += "Or the tree.shared command can be in the following format: tree.shared(phylip=yourPhylipFile).\n";
		helpString += "Example tree.shared(phylip=abrecovery.dist).\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
TreeGroupCommand::TreeGroupCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "TreeGroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

TreeGroupCommand::TreeGroupCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["tree"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  inputfile = phylipfile;  format = "phylip"; m->setPhylipFile(phylipfile);	}
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  inputfile = columnfile; format = "column";	m->setColumnFile(columnfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  inputfile = sharedfile; format = "sharedfile";	m->setSharedFile(sharedfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			if ((phylipfile == "") && (columnfile == "") && (sharedfile == "")) { 
				//is there are current file available for either of these?
				//give priority to shared, then column, then phylip
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") {  inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					columnfile = m->getColumnFile(); 
					if (columnfile != "") { inputfile = columnfile; format = "column";  m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
					else { 
						phylipfile = m->getPhylipFile(); 
						if (phylipfile != "") { inputfile = phylipfile;  format = "phylip";  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("No valid current files. You must provide a shared, phylip or column file."); m->mothurOutEndLine(); 
							abort = true;
						}
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When running the tree.shared command with a distance file you may not use both the column and the phylip parameters."); m->mothurOutEndLine(); abort = true; }
			
			if (columnfile != "") {
				if (namefile == "") { 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a namefile if you are going to use the column format."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				}
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->Groups = Groups;
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}

			string temp;
			temp = validParameter.validFile(parameters, "precision", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(inputfile); //if user entered a file with a path then preserve it	
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "TreeGroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

TreeGroupCommand::~TreeGroupCommand(){
	if (abort == false) {
		if (format == "sharedfile") {  delete input; }
		else { delete readMatrix;  delete matrix; delete list; }
		delete tmap;  
	}
	
}

//**********************************************************************************************************************

int TreeGroupCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if (format == "sharedfile") {
			
			ValidCalculators validCalculator;
		
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("treegroup", Estimators[i]) == true) { 
					if (Estimators[i] == "sharedsobs") { 
						treeCalculators.push_back(new SharedSobsCS());
					}else if (Estimators[i] == "sharedchao") { 
						treeCalculators.push_back(new SharedChao1());
					}else if (Estimators[i] == "sharedace") { 
						treeCalculators.push_back(new SharedAce());
					}else if (Estimators[i] == "jabund") { 	
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
					}else if (Estimators[i] == "kstest") { 
						treeCalculators.push_back(new KSTest());
					}else if (Estimators[i] == "sharednseqs") { 
						treeCalculators.push_back(new SharedNSeqs());
					}else if (Estimators[i] == "ochiai") { 
						treeCalculators.push_back(new Ochiai());
					}else if (Estimators[i] == "anderberg") { 
						treeCalculators.push_back(new Anderberg());
					}else if (Estimators[i] == "kulczynski") { 
						treeCalculators.push_back(new Kulczynski());
					}else if (Estimators[i] == "kulczynskicody") { 
						treeCalculators.push_back(new KulczynskiCody());
					}else if (Estimators[i] == "lennon") { 
						treeCalculators.push_back(new Lennon());
					}else if (Estimators[i] == "morisitahorn") { 
						treeCalculators.push_back(new MorHorn());
					}else if (Estimators[i] == "braycurtis") { 
						treeCalculators.push_back(new BrayCurtis());
					}else if (Estimators[i] == "whittaker") { 
						treeCalculators.push_back(new Whittaker());
					}else if (Estimators[i] == "odum") { 
						treeCalculators.push_back(new Odum());
					}else if (Estimators[i] == "canberra") { 
						treeCalculators.push_back(new Canberra());
					}else if (Estimators[i] == "structeuclidean") { 
						treeCalculators.push_back(new StructEuclidean());
					}else if (Estimators[i] == "structchord") { 
						treeCalculators.push_back(new StructChord());
					}else if (Estimators[i] == "hellinger") { 
						treeCalculators.push_back(new Hellinger());
					}else if (Estimators[i] == "manhattan") { 
						treeCalculators.push_back(new Manhattan());
					}else if (Estimators[i] == "structpearson") { 
						treeCalculators.push_back(new StructPearson());
					}else if (Estimators[i] == "soergel") { 
						treeCalculators.push_back(new Soergel());
					}else if (Estimators[i] == "spearman") { 
						treeCalculators.push_back(new Spearman());
					}else if (Estimators[i] == "structkulczynski") { 
						treeCalculators.push_back(new StructKulczynski());
					}else if (Estimators[i] == "speciesprofile") { 
						treeCalculators.push_back(new SpeciesProfile());
					}else if (Estimators[i] == "hamming") { 
						treeCalculators.push_back(new Hamming());
					}else if (Estimators[i] == "structchi2") { 
						treeCalculators.push_back(new StructChi2());
					}else if (Estimators[i] == "gower") { 
						treeCalculators.push_back(new Gower());
					}else if (Estimators[i] == "memchi2") { 
						treeCalculators.push_back(new MemChi2());
					}else if (Estimators[i] == "memchord") { 
						treeCalculators.push_back(new MemChord());
					}else if (Estimators[i] == "memeuclidean") { 
						treeCalculators.push_back(new MemEuclidean());
					}else if (Estimators[i] == "mempearson") { 
						treeCalculators.push_back(new MemPearson());
					}
				}
			}
			
			//if the users entered no valid calculators don't execute command
			if (treeCalculators.size() == 0) { m->mothurOut("You have given no valid calculators."); m->mothurOutEndLine(); return 0; }
			
			input = new InputData(sharedfile, "sharedfile");
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup[0]->getLabel();
			
			if (lookup.size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command."); m->mothurOutEndLine(); return 0; }
			
			//used in tree constructor 
			m->runParse = false;
			
			//create treemap class from groupmap for tree class to use
			tmap = new TreeMap();
			tmap->makeSim(m->namesOfGroups);
			
			//clear globaldatas old tree names if any
			m->Treenames.clear();
			
			//fills globaldatas tree names
			m->Treenames = m->Groups;
		
			if (m->control_pressed) { return 0; }
			
			//create tree file
			makeSimsShared();
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
		}else{
			//read in dist file
			filename = inputfile;
		
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
			
			if (m->control_pressed) { return 0; }
			
			tmap->makeSim(list);
			
			m->Groups = tmap->namesOfGroups;
		
			//clear globaldatas old tree names if any
			m->Treenames.clear();
		
			//fills globaldatas tree names
			m->Treenames = m->Groups;
			
			//used in tree constructor 
			m->runParse = false;
			
			if (m->control_pressed) { return 0; }
			
			makeSimsDist();
			
			if (m->control_pressed) { return 0; }

			//create a new filename
			outputFile = outputDir + m->getRootName(m->getSimpleName(inputfile)) + "tre";	
			outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
				
			createTree();
			
			if (m->control_pressed) { return 0; }

			m->mothurOut("Tree complete. "); m->mothurOutEndLine();
			
		}
				
		//reset groups parameter
		m->Groups.clear(); 
		
		//set tree file as new current treefile
		string current = "";
		itTypes = outputTypes.find("tree");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTreeFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int TreeGroupCommand::createTree(){
	try {
		//create tree
		t = new Tree(tmap);
		
		//do merges and create tree structure by setting parents and children
		//there are numGroups - 1 merges to do
		for (int i = 0; i < (numGroups - 1); i++) {
			float largest = -1000.0;
			
			if (m->control_pressed) { delete t; return 1; }
			
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
		
		if (m->control_pressed) { delete t; return 1; }
		
		//print newick file
		t->createNewickFile(outputFile);
		
		//delete tree
		delete t;
		
		if (m->control_pressed) { remove(outputFile.c_str()); outputNames.pop_back(); return 1; }
		
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "createTree");
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
		m->errorOut(e, "TreeGroupCommand", "printSims");
		exit(1);
	}
}
/***********************************************************/
int TreeGroupCommand::makeSimsDist() {
	try {
		numGroups = list->size();
		
		//initialize index
		index.clear();
		for (int g = 0; g < numGroups; g++) {	index[g] = g;	}
		
		//initialize simMatrix
		simMatrix.clear();
		simMatrix.resize(numGroups);
		for (int k = 0; k < simMatrix.size(); k++)	{
			for (int j = 0; j < simMatrix.size(); j++)	{
				simMatrix[k].push_back(0.0);
			}
		}
		
		//go through sparse matrix and fill sims
		//go through each cell in the sparsematrix
		for(MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++){
			//similairity = -(distance-1)
			simMatrix[currentCell->row][currentCell->column] = -(currentCell->dist -1.0);	
			simMatrix[currentCell->column][currentCell->row] = -(currentCell->dist -1.0);	
			
			if (m->control_pressed) { return 1; }
			
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "makeSimsDist");
		exit(1);
	}
}

/***********************************************************/
int TreeGroupCommand::makeSimsShared() {
	try {
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } for(int i = 0 ; i < treeCalculators.size(); i++) {  delete treeCalculators[i]; } return 1; }
		
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
			
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}

			lastLabel = lookup[0]->getLabel();			
			
			//get next line to process
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } for(int i = 0 ; i < treeCalculators.size(); i++) {  delete treeCalculators[i]; } return 1; }

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
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {		delete lookup[i]; }		} 
			lookup = input->getSharedRAbundVectors(lastLabel);

			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			process(lookup);
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 	
		}
		
		for(int i = 0 ; i < treeCalculators.size(); i++) {  delete treeCalculators[i]; }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "makeSimsShared");
		exit(1);
	}
}

/***********************************************************/
int TreeGroupCommand::process(vector<SharedRAbundVector*> thisLookup) {
	try{
				EstOutput data;
				vector<SharedRAbundVector*> subset;
				numGroups = thisLookup.size();
				
				//for each calculator												
				for(int i = 0 ; i < treeCalculators.size(); i++) {
					//initialize simMatrix
					simMatrix.clear();
					simMatrix.resize(numGroups);
					for (int k = 0; k < simMatrix.size(); k++)	{
						for (int j = 0; j < simMatrix.size(); j++)	{
							simMatrix[k].push_back(0.0);
						}
					}
		
					//initialize index
					index.clear();
					for (int g = 0; g < numGroups; g++) {	index[g] = g;	}
		
					//create a new filename
					outputFile = outputDir + m->getRootName(m->getSimpleName(inputfile)) + treeCalculators[i]->getName() + "." + thisLookup[0]->getLabel() + ".tre";				
					outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
												
					for (int k = 0; k < thisLookup.size(); k++) { 
						for (int l = k; l < thisLookup.size(); l++) {
							if (k != l) { //we dont need to similiarity of a groups to itself
								//get estimated similarity between 2 groups
								
								subset.clear(); //clear out old pair of sharedrabunds
								//add new pair of sharedrabunds
								subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
								
								//if this calc needs all groups to calculate the pair load all groups
								if (treeCalculators[i]->getNeedsAll()) { 
									//load subset with rest of lookup for those calcs that need everyone to calc for a pair
									for (int w = 0; w < thisLookup.size(); w++) {
										if ((w != k) && (w != l)) { subset.push_back(thisLookup[w]); }
									}
								}
								
								data = treeCalculators[i]->getValues(subset); //saves the calculator outputs
						//cout << thisLookup[k]->getGroup() << '\t' << thisLookup[l]->getGroup() << '\t' << (1.0 - data[0]) << endl;
								if (m->control_pressed) { return 1; }
								
								//save values in similarity matrix
								simMatrix[k][l] = -(data[0]-1.0);
								simMatrix[l][k] = -(data[0]-1.0);
							}
						}
					}
					
					//createdistance file from simMatrix
					/*string o = outputDir + m->getRootName(m->getSimpleName(globaldata->inputFileName)) + treeCalculators[i]->getName() + "." + thisLookup[0]->getLabel() + ".dist";
					ofstream outDist;
					m->openOutputFile(o, outDist);
					outDist << simMatrix.size() << endl;
					for (int k = 0; k < simMatrix.size(); k++) {
						outDist << thisLookup[k]->getGroup() << '\t';
						for (int l = 0; l < k; l++) {
							outDist << (1.0-simMatrix[k][l]) << '\t';
						}
						outDist << endl;
					}
					outDist.close();*/

					
					if (m->control_pressed) { return 1; }
					//creates tree from similarity matrix and write out file
					createTree();
					
					if (m->control_pressed) { return 1; }
				}
				
				return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "process");
		exit(1);
	}
}
/***********************************************************/

	

