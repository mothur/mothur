/*
 *  treegroupscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treegroupscommand.h"
#include "subsample.h"
#include "consensus.h"

//**********************************************************************************************************************
vector<string> TreeGroupCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none","tree",false,false,true); parameters.push_back(pshared);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none","tree",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName","",false,false); parameters.push_back(pname);
		CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "countcolumn","",false,false); parameters.push_back(pcount);
        CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "ColumnName-countcolumn","tree",false,false); parameters.push_back(pcolumn);		
        CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","",false,false); parameters.push_back(psubsample);
        CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);		
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pcalc("calc", "Multiple", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-sharednseqs-ochiai-anderberg-kulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-whittaker-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-hamming-structchi2-gower-memchi2-memchord-memeuclidean-mempearson-jsd-rjsd", "jclass-thetayc", "", "", "","",true,false,true); parameters.push_back(pcalc);
		
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
//CommandParameter poutput("output", "Multiple", "lt-square", "lt", "", "", "",false,false); parameters.push_back(poutput);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The tree.shared command creates a .tre to represent the similarity between groups or sequences.\n";
		helpString += "The tree.shared command parameters are shared, groups, calc, phylip, column, name, cutoff, precision, processors, subsample, iters and label.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n";
		helpString += "The group names are separated by dashes. The label allow you to select what distance levels you would like trees created for, and are also separated by dashes.\n";
		helpString += "The phylip or column parameter are required if you do not provide a sharedfile, and only one may be used.  If you use a column file the name filename is required. \n";
		helpString += "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.\n";
		helpString += "The tree.shared command should be in the following format: tree.shared(groups=yourGroups, calc=yourCalcs, label=yourLabels).\n";
        helpString += "The iters parameter allows you to choose the number of times you would like to run the subsample.\n";
        helpString += "The subsample parameter allows you to enter the size pergroup of the sample or you can set subsample=T and mothur will use the size of your smallest group. The subsample parameter may only be used with a shared file.\n";
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
string TreeGroupCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "tree") {  pattern = "[filename],[calc],[distance],[tag],tre-[filename],tre"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "TreeGroupCommand", "getOutputPattern");
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
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
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
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
				if ((namefile == "") && (countfile == "")){ 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = m->getCountTableFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the column format."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
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
				m->setGroups(Groups);
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
			m->mothurConvert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "10"; }
			m->mothurConvert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
            temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
            
            temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, iters); 
            
            temp = validParameter.validFile(parameters, "subsample", false);		if (temp == "not found") { temp = "F"; }
			if (m->isNumeric1(temp)) { m->mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (m->isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
            
            if (subsample == false) { iters = 1; }
            
            if (subsample && (format != "sharedfile")) { m->mothurOut("[ERROR]: the subsample parameter can only be used with a shared file.\n"); abort=true; }
            
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
		else { delete list; }
		delete ct;  
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
                    }else if (Estimators[i] == "jsd") {
                        treeCalculators.push_back(new JSD());
                    }else if (Estimators[i] == "rjsd") {
                        treeCalculators.push_back(new RJSD());
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
			ct = new CountTable();
            set<string> nameMap;
            map<string, string> groupMap;
            set<string> gps;
            for (int i = 0; i < m->getAllGroups().size(); i++) { 
                nameMap.insert(m->getAllGroups()[i]); 
                gps.insert(m->getAllGroups()[i]); 
                groupMap[m->getAllGroups()[i]] = m->getAllGroups()[i];
            }
            ct->createTable(nameMap, groupMap, gps);
			
			//clear globaldatas old tree names if any
			m->Treenames.clear();
			
			//fills globaldatas tree names
			//m->Treenames = m->getGroups();
            for (int k = 0; k < lookup.size(); k++) {
                m->Treenames.push_back(lookup[k]->getGroup());
            }
		
			if (m->control_pressed) { return 0; }
			
			//create tree file
			makeSimsShared();
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
		}else{
			//read in dist file
			filename = inputfile;
            
            ReadMatrix* readMatrix;
			if (format == "column") { readMatrix = new ReadColumnMatrix(filename); }	
			else if (format == "phylip") { readMatrix = new ReadPhylipMatrix(filename); }
				
			readMatrix->setCutoff(cutoff);
	
            ct = NULL;
            nameMap = NULL;
            if(namefile != ""){	
                nameMap = new NameAssignment(namefile);
                nameMap->readMap();
                readMatrix->read(nameMap);
            }else if (countfile != "") {
                ct = new CountTable();
                ct->readTable(countfile, true, false);
                readMatrix->read(ct);
            }else {
                readMatrix->read(nameMap);
            }

			list = readMatrix->getListVector();
			SparseDistanceMatrix* dMatrix = readMatrix->getDMatrix();
            
            //clear globaldatas old tree names if any
			m->Treenames.clear();
            
			//make treemap
            if (ct != NULL) { delete ct; }
			ct = new CountTable();
            set<string> nameMap;
            map<string, string> groupMap;
            set<string> gps;
            for (int i = 0; i < list->getNumBins(); i++) {
                string bin = list->get(i);
                nameMap.insert(bin); 
                gps.insert(bin); 
                groupMap[bin] = bin;
                m->Treenames.push_back(bin);
            }
            ct->createTable(nameMap, groupMap, gps);
			
			vector<string> namesGroups = ct->getNamesOfGroups();
			m->setGroups(namesGroups);
			
			//used in tree constructor 
			m->runParse = false;
			
			if (m->control_pressed) { return 0; }
			
			vector< vector<double> > matrix = makeSimsDist(dMatrix);
            delete readMatrix;
            delete dMatrix;
			
			if (m->control_pressed) { return 0; }

			//create a new filename
            map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
			string outputFile = getOutputFileName("tree",variables);	
			outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
				
			Tree* newTree = createTree(matrix);
            
            if (newTree != NULL) {  writeTree(outputFile, newTree); delete newTree; }
			
			if (m->control_pressed) { return 0; }

			m->mothurOut("Tree complete. "); m->mothurOutEndLine();
			
		}
				
		//reset groups parameter
		m->clearGroups(); 
		
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

Tree* TreeGroupCommand::createTree(vector< vector<double> >& simMatrix){
	try {
		//create tree
		t = new Tree(ct, simMatrix);
        
        if (m->control_pressed) { delete t; t = NULL; return t; }
		
        //assemble tree
		t->assembleTree();

		return t;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "createTree");
		exit(1);
	}
}
/***********************************************************/
int TreeGroupCommand::writeTree(string out, Tree* T) {
	try {
		
        //print newick file
		t->createNewickFile(out);
		
        if (m->control_pressed) { m->mothurRemove(out); outputNames.pop_back(); return 1; }
        
        return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "printSims");
		exit(1);
	}
}

/***********************************************************/
void TreeGroupCommand::printSims(ostream& out, vector< vector<double> >& simMatrix) {
	try {
		
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
vector< vector<double> > TreeGroupCommand::makeSimsDist(SparseDistanceMatrix* matrix) {
	try {
		numGroups = list->size();
		
		//initialize simMatrix
		vector< vector<double> > simMatrix;
		simMatrix.resize(numGroups);
		for (int k = 0; k < simMatrix.size(); k++)	{
			for (int j = 0; j < simMatrix.size(); j++)	{
				simMatrix[k].push_back(0.0);
			}
		}
		
		//go through sparse matrix and fill sims
		//go through each cell in the sparsematrix
        for (int i = 0; i < matrix->seqVec.size(); i++) {
            for (int j = 0; j < matrix->seqVec[i].size(); j++) {
                
                //already checked everyone else in row
                if (i < matrix->seqVec[i][j].index) {   
                    simMatrix[i][matrix->seqVec[i][j].index] = -(matrix->seqVec[i][j].dist -1.0);	
                    simMatrix[matrix->seqVec[i][j].index][i] = -(matrix->seqVec[i][j].dist -1.0);	
			
                    if (m->control_pressed) { return simMatrix; }
                }
            }
		}

		return simMatrix;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "makeSimsDist");
		exit(1);
	}
}

/***********************************************************/
int TreeGroupCommand::makeSimsShared() {
	try {
        
        if (subsample) { 
            if (subsampleSize == -1) { //user has not set size, set size = smallest samples size
                subsampleSize = lookup[0]->getNumSeqs();
                for (int i = 1; i < lookup.size(); i++) {
                    int thisSize = lookup[i]->getNumSeqs();
                    
                    if (thisSize < subsampleSize) {	subsampleSize = thisSize;	}
                }
            }else {
                m->clearGroups();
                Groups.clear();
                m->Treenames.clear();
                vector<SharedRAbundVector*> temp;
                for (int i = 0; i < lookup.size(); i++) {
                    if (lookup[i]->getNumSeqs() < subsampleSize) { 
                        m->mothurOut(lookup[i]->getGroup() + " contains " + toString(lookup[i]->getNumSeqs()) + ". Eliminating."); m->mothurOutEndLine();
                        delete lookup[i];
                    }else { 
                        Groups.push_back(lookup[i]->getGroup()); 
                        temp.push_back(lookup[i]);
                        m->Treenames.push_back(lookup[i]->getGroup());
                    }
                } 
                lookup = temp;
                m->setGroups(Groups);
            }
            
            if (lookup.size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command."); m->mothurOutEndLine(); m->control_pressed = true; return 0; }
        }
        numGroups = lookup.size();
        
        //sanity check to make sure processors < numComparisions
        int numDists = 0;
        for(int i=0;i<numGroups;i++){
            for(int j=0;j<i;j++){
                numDists++;
                if (numDists > processors) { break; }
            }
        }
        if (numDists < processors) { processors = numDists; }
        
		lines.resize(processors);
		for (int i = 0; i < processors; i++) {
			lines[i].start = int (sqrt(float(i)/float(processors)) * numGroups);
			lines[i].end = int (sqrt(float(i+1)/float(processors)) * numGroups);
		}	
        
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
		vector< vector< vector<seqDist> > > calcDistsTotals;  //each iter, one for each calc, then each groupCombos dists. this will be used to make .dist files
        vector< vector<seqDist>  > calcDists; calcDists.resize(treeCalculators.size()); 		
        
        for (int thisIter = 0; thisIter < iters; thisIter++) {
            
            vector<SharedRAbundVector*> thisItersLookup = thisLookup;
            
            if (subsample) {
                SubSample sample;
                vector<string> tempLabels; //dont need since we arent printing the sampled sharedRabunds
                
                //make copy of lookup so we don't get access violations
                vector<SharedRAbundVector*> newLookup;
                for (int k = 0; k < thisItersLookup.size(); k++) {
                    SharedRAbundVector* temp = new SharedRAbundVector();
                    temp->setLabel(thisItersLookup[k]->getLabel());
                    temp->setGroup(thisItersLookup[k]->getGroup());
                    newLookup.push_back(temp);
                }
                
                //for each bin
                for (int k = 0; k < thisItersLookup[0]->getNumBins(); k++) {
                    if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
                    for (int j = 0; j < thisItersLookup.size(); j++) { newLookup[j]->push_back(thisItersLookup[j]->getAbundance(k), thisItersLookup[j]->getGroup()); }
                }
                
                tempLabels = sample.getSample(newLookup, subsampleSize);
                thisItersLookup = newLookup;
            }
            
            if(processors == 1){
                driver(thisItersLookup, 0, numGroups, calcDists);
            }else{
                int process = 1;
                vector<int> processIDS;
                bool recalc = false;
                
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                //loop through and create all the processes you want
                while (process != processors) {
                    pid_t pid = fork();
                    
                    if (pid > 0) {
                        processIDS.push_back(pid); 
                        process++;
                    }else if (pid == 0){
                        
                        driver(thisItersLookup, lines[process].start, lines[process].end, calcDists);   
                        
                        string tempdistFileName = m->getRootName(m->getSimpleName(sharedfile)) + m->mothurGetpid(process) + ".dist";
                        ofstream outtemp;
                        m->openOutputFile(tempdistFileName, outtemp);
                        
                        for (int i = 0; i < calcDists.size(); i++) {
                            outtemp << calcDists[i].size() << endl;
                            
                            for (int j = 0; j < calcDists[i].size(); j++) {
                                outtemp << calcDists[i][j].seq1 << '\t' << calcDists[i][j].seq2 << '\t' << calcDists[i][j].dist << endl;
                            }
                        }
                        outtemp.close();
                        
                        exit(0);
                    }else { 
                        m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                        for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                        //wait to die
                        for (int i=0;i<processIDS.size();i++) {
                            int temp = processIDS[i];
                            wait(&temp);
                        }
                        m->control_pressed = false;
                        for (int i=0;i<processIDS.size();i++) {
                            m->mothurRemove(m->getRootName(m->getSimpleName(sharedfile)) + (toString(processIDS[i]) + ".dist"));
                        }
                        recalc = true;
                        break;
                    }
                }
                
                if (recalc) {
                    //test line, also set recalc to true.
                    //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(m->getRootName(m->getSimpleName(sharedfile)) + (toString(processIDS[i]) + ".dist"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
                    
                    /******************************************************/
                    //comparison breakup to be used by different processes later
                    lines.clear();
                    numGroups = thisLookup.size();
                    lines.resize(processors);
                    for (int i = 0; i < processors; i++) {
                        lines[i].start = int (sqrt(float(i)/float(processors)) * numGroups);
                        lines[i].end = int (sqrt(float(i+1)/float(processors)) * numGroups);
                    }
                    /******************************************************/
                    
                    calcDists.clear(); calcDists.resize(treeCalculators.size());
                    processIDS.resize(0);
                    process = 1;
                    
                    //loop through and create all the processes you want
                    while (process != processors) {
                        pid_t pid = fork();
                        
                        if (pid > 0) {
                            processIDS.push_back(pid);
                            process++;
                        }else if (pid == 0){
                            
                            driver(thisItersLookup, lines[process].start, lines[process].end, calcDists);
                            
                            string tempdistFileName = m->getRootName(m->getSimpleName(sharedfile)) + m->mothurGetpid(process) + ".dist";
                            ofstream outtemp;
                            m->openOutputFile(tempdistFileName, outtemp);
                            
                            for (int i = 0; i < calcDists.size(); i++) {
                                outtemp << calcDists[i].size() << endl;
                                
                                for (int j = 0; j < calcDists[i].size(); j++) {
                                    outtemp << calcDists[i][j].seq1 << '\t' << calcDists[i][j].seq2 << '\t' << calcDists[i][j].dist << endl;
                                }
                            }
                            outtemp.close();
                            
                            exit(0);
                        }else {
                            m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                            for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                            exit(0);
                        }
                    }                    
                }

                //parent do your part
                driver(thisItersLookup, lines[0].start, lines[0].end, calcDists);   
                
                //force parent to wait until all the processes are done
                for (int i = 0; i < processIDS.size(); i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                
                for (int i = 0; i < processIDS.size(); i++) {
                    string tempdistFileName = m->getRootName(m->getSimpleName(sharedfile)) + toString(processIDS[i]) +  ".dist";
                    ifstream intemp;
                    m->openInputFile(tempdistFileName, intemp);
                    
                    for (int k = 0; k < calcDists.size(); k++) {
                        int size = 0;
                        intemp >> size; m->gobble(intemp);
                        
                        for (int j = 0; j < size; j++) {
                            int seq1 = 0;
                            int seq2 = 0;
                            float dist = 1.0;
                            
                            intemp >> seq1 >> seq2 >> dist;   m->gobble(intemp);
                            
                            seqDist tempDist(seq1, seq2, dist);
                            calcDists[k].push_back(tempDist);
                        }
                    }
                    intemp.close();
                    m->mothurRemove(tempdistFileName);
                }
#else
                //////////////////////////////////////////////////////////////////////////////////////////////////////
                //Windows version shared memory, so be careful when passing variables through the treeSharedData struct. 
                //Above fork() will clone, so memory is separate, but that's not the case with windows, 
                //Taking advantage of shared memory to pass results vectors.
                //////////////////////////////////////////////////////////////////////////////////////////////////////
                
                vector<treeSharedData*> pDataArray; 
                DWORD   dwThreadIdArray[processors-1];
                HANDLE  hThreadArray[processors-1]; 
                
                //Create processor worker threads.
                for( int i=1; i<processors; i++ ){
                    
                    //make copy of lookup so we don't get access violations
                    vector<SharedRAbundVector*> newLookup;
                    for (int k = 0; k < thisItersLookup.size(); k++) {
                        SharedRAbundVector* temp = new SharedRAbundVector();
                        temp->setLabel(thisItersLookup[k]->getLabel());
                        temp->setGroup(thisItersLookup[k]->getGroup());
                        newLookup.push_back(temp);
                    }
                    
                    //for each bin
                    for (int k = 0; k < thisItersLookup[0]->getNumBins(); k++) {
                        if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
                        for (int j = 0; j < thisItersLookup.size(); j++) { newLookup[j]->push_back(thisItersLookup[j]->getAbundance(k), thisItersLookup[j]->getGroup()); }
                    }
                    
                    // Allocate memory for thread data.
                    treeSharedData* tempSum = new treeSharedData(m, lines[i].start, lines[i].end, Estimators, newLookup);
                    pDataArray.push_back(tempSum);
                    processIDS.push_back(i);
                    
                    hThreadArray[i-1] = CreateThread(NULL, 0, MyTreeSharedThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
                }
                
                //parent do your part
                driver(thisItersLookup, lines[0].start, lines[0].end, calcDists);   
                
                //Wait until all threads have terminated.
                WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
                
                //Close all thread handles and free memory allocations.
                for(int i=0; i < pDataArray.size(); i++){
                    if (pDataArray[i]->count != (pDataArray[i]->end-pDataArray[i]->start)) {
                        m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end-pDataArray[i]->start) + " groups assigned to it, quitting. \n"); m->control_pressed = true; 
                    }
                    for (int j = 0; j < pDataArray[i]->thisLookup.size(); j++) {  delete pDataArray[i]->thisLookup[j];  } 
                    
                    for (int k = 0; k < calcDists.size(); k++) {
                        int size = pDataArray[i]->calcDists[k].size();
                        for (int j = 0; j < size; j++) {    calcDists[k].push_back(pDataArray[i]->calcDists[k][j]);    }
                    }
                    
                    CloseHandle(hThreadArray[i]);
                    delete pDataArray[i];
                }
                
#endif
            }
            
            calcDistsTotals.push_back(calcDists);
            
            if (subsample) {  
                
                //clean up memory
                for (int i = 0; i < thisItersLookup.size(); i++) { delete thisItersLookup[i]; }
                thisItersLookup.clear();
                for (int i = 0; i < calcDists.size(); i++) {  calcDists[i].clear(); }
            }
            
            if (m->debug) {  m->mothurOut("[DEBUG]: iter = " + toString(thisIter) + ".\n"); }
		}
        
		if (m->debug) {  m->mothurOut("[DEBUG]: done with iters.\n"); }
            
        if (iters != 1) {
            //we need to find the average distance and standard deviation for each groups distance
            vector< vector<seqDist>  > calcAverages = m->getAverages(calcDistsTotals);  
            
            if (m->debug) {  m->mothurOut("[DEBUG]: found averages.\n"); }
            
            //create average tree for each calc
            for (int i = 0; i < calcDists.size(); i++) {
                vector< vector<double> > matrix; //square matrix to represent the distance
                matrix.resize(thisLookup.size());
                for (int k = 0; k < thisLookup.size(); k++) {  matrix[k].resize(thisLookup.size(), 0.0); }
                
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    int row = calcAverages[i][j].seq1;
                    int column = calcAverages[i][j].seq2;
                    float dist = calcAverages[i][j].dist;
                    
                    matrix[row][column] = dist;
                    matrix[column][row] = dist;
                }
                
                //create a new filename
                map<string, string> variables; 
                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
                variables["[calc]"] = treeCalculators[i]->getName();
                variables["[distance]"] = thisLookup[0]->getLabel();
                variables["[tag]"] = "ave";
                string outputFile = getOutputFileName("tree",variables);				
                outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
                
                //creates tree from similarity matrix and write out file
                Tree* newTree = createTree(matrix);
                if (newTree != NULL) { writeTree(outputFile, newTree); }                
            }
            
            if (m->debug) {  m->mothurOut("[DEBUG]: done averages trees.\n"); }
            
            //create all trees for each calc and find their consensus tree
            for (int i = 0; i < calcDists.size(); i++) {
                if (m->control_pressed) { break; }
                
                //create a new filename
                //create a new filename
                map<string, string> variables; 
                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
                variables["[calc]"] = treeCalculators[i]->getName();
                variables["[distance]"] = thisLookup[0]->getLabel();
                variables["[tag]"] = "all";
                string outputFile = getOutputFileName("tree",variables);				
                outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
                
                ofstream outAll;
                m->openOutputFile(outputFile, outAll);
                
                vector<Tree*> trees; 
                for (int myIter = 0; myIter < iters; myIter++) {
                    
                    if(m->control_pressed) { break; }
                    
                    //initialize matrix
                    vector< vector<double> > matrix; //square matrix to represent the distance
                    matrix.resize(thisLookup.size());
                    for (int k = 0; k < thisLookup.size(); k++) {  matrix[k].resize(thisLookup.size(), 0.0); }
                    
                    for (int j = 0; j < calcDistsTotals[myIter][i].size(); j++) {
                        int row = calcDistsTotals[myIter][i][j].seq1;
                        int column = calcDistsTotals[myIter][i][j].seq2;
                        double dist = calcDistsTotals[myIter][i][j].dist;
                       
                        matrix[row][column] = dist;
                        matrix[column][row] = dist;
                    }
                    
                    //creates tree from similarity matrix and write out file
                    Tree* newTree = createTree(matrix);
                    if (newTree != NULL) { 
                        newTree->print(outAll);
                        trees.push_back(newTree);
                    }
                }
                outAll.close();
                if (m->control_pressed) { for (int k = 0; k < trees.size(); k++) { delete trees[k]; } }
                
                if (m->debug) {  m->mothurOut("[DEBUG]: done all trees.\n"); }
                
                Consensus consensus;
                //clear old tree names if any
                m->Treenames.clear(); m->Treenames = m->getGroups(); //may have changed if subsample eliminated groups
                Tree* conTree = consensus.getTree(trees);
                
                if (m->debug) {  m->mothurOut("[DEBUG]: done cons tree.\n"); }
                
                //create a new filename
                variables["[tag]"] = "cons";
                string conFile = getOutputFileName("tree",variables);				
                outputNames.push_back(conFile); outputTypes["tree"].push_back(conFile); 
                ofstream outTree;
                m->openOutputFile(conFile, outTree);
                
                if (conTree != NULL) { conTree->print(outTree, "boot"); delete conTree; }
            }

        }else {
            
            for (int i = 0; i < calcDists.size(); i++) {
                if (m->control_pressed) { break; }
                
                //initialize matrix
                vector< vector<double> > matrix; //square matrix to represent the distance
                matrix.resize(thisLookup.size());
                for (int k = 0; k < thisLookup.size(); k++) {  matrix[k].resize(thisLookup.size(), 0.0); }
                
                for (int j = 0; j < calcDists[i].size(); j++) {
                    int row = calcDists[i][j].seq1;
                    int column = calcDists[i][j].seq2;
                    double dist = calcDists[i][j].dist;
                    
                    matrix[row][column] = dist;
                    matrix[column][row] = dist;
                }
                
                //create a new filename
                map<string, string> variables; 
                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
                variables["[calc]"] = treeCalculators[i]->getName();
                variables["[distance]"] = thisLookup[0]->getLabel();
                variables["[tag]"] = "";
                string outputFile = getOutputFileName("tree",variables);					
                outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
                
                //creates tree from similarity matrix and write out file
                Tree* newTree = createTree(matrix);
                if (newTree != NULL) { writeTree(outputFile, newTree); delete newTree; }
            }
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
int TreeGroupCommand::driver(vector<SharedRAbundVector*> thisLookup, int start, int end, vector< vector<seqDist> >& calcDists) { 
	try {
		vector<SharedRAbundVector*> subset;
		for (int k = start; k < end; k++) { // pass cdd each set of groups to compare
			
			for (int l = 0; l < k; l++) {
				
				if (k != l) { //we dont need to similiarity of a groups to itself
					subset.clear(); //clear out old pair of sharedrabunds
					//add new pair of sharedrabunds
					subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
					
					for(int i=0;i<treeCalculators.size();i++) {
						
						//if this calc needs all groups to calculate the pair load all groups
						if (treeCalculators[i]->getNeedsAll()) { 
							//load subset with rest of lookup for those calcs that need everyone to calc for a pair
							for (int w = 0; w < thisLookup.size(); w++) {
								if ((w != k) && (w != l)) { subset.push_back(thisLookup[w]); }
							}
						}
						
						vector<double> tempdata = treeCalculators[i]->getValues(subset); //saves the calculator outputs
						
						if (m->control_pressed) { return 1; }
						
						seqDist temp(l, k, -(tempdata[0]-1.0));
						calcDists[i].push_back(temp);
					}
				}
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeGroupCommand", "driver");
		exit(1);
	}
}
/***********************************************************/

	

