/*
 *  treegroupscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treesharedcommand.h"
#include "subsample.h"
#include "consensus.h"

//**********************************************************************************************************************
vector<string> TreeSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none","tree",false,false,true); parameters.push_back(pshared);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none","tree",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName","",false,false); parameters.push_back(pname);
		CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "countcolumn","",false,false); parameters.push_back(pcount);
        CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "ColumnName-countcolumn","tree",false,false); parameters.push_back(pcolumn);		
        CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","",false,false); parameters.push_back(psubsample);
        CommandParameter pwithreplacement("withreplacement", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pwithreplacement);
        CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);		
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pcalc("calc", "Multiple", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-sharednseqs-ochiai-anderberg-kulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-whittaker-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-hamming-structchi2-gower-memchi2-memchord-memeuclidean-mempearson-jsd-rjsd", "jclass-thetayc", "", "", "","",true,false,true); parameters.push_back(pcalc);
		
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["tree"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string TreeSharedCommand::getHelpString(){	
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
        helpString += "The withreplacement parameter allows you to indicate you want to subsample your data allowing for the same read to be included multiple times. Default=f. \n";
		helpString += "Example tree.shared(groups=A-B-C, calc=jabund-sorabund).\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The default value for calc is jclass-thetayc.\n";
		helpString += "The tree.shared command outputs a .tre file for each calculator you specify at each distance you choose.\n";
		helpString += validCalculator.printCalc("treegroup");
		helpString += "Or the tree.shared command can be in the following format: tree.shared(phylip=yourPhylipFile).\n";
		helpString += "Example tree.shared(phylip=abrecovery.dist).\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string TreeSharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "tree") {  pattern = "[filename],[calc],[distance],[tag],tre-[filename],tre"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "TreeSharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
TreeSharedCommand::TreeSharedCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
			phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  inputfile = phylipfile;  format = "phylip"; current->setPhylipFile(phylipfile);	}
			
			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  inputfile = columnfile; format = "column";	current->setColumnFile(columnfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  inputfile = sharedfile; format = "sharedfile";	current->setSharedFile(sharedfile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
			if ((phylipfile == "") && (columnfile == "") && (sharedfile == "")) { 
				//is there are current file available for either of these?
				//give priority to shared, then column, then phylip
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") {  inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 
					columnfile = current->getColumnFile(); 
					if (columnfile != "") { inputfile = columnfile; format = "column";  m->mothurOut("Using " + columnfile + " as input file for the column parameter.\n");  }
					else { 
						phylipfile = current->getPhylipFile(); 
						if (phylipfile != "") { inputfile = phylipfile;  format = "phylip";  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter.\n");  }
						else { 
							m->mothurOut("No valid current files. You must provide a shared, phylip or column file.\n");
							abort = true;
						}
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When running the tree.shared command with a distance file you may not use both the column and the phylip parameters.\n");  abort = true; }
			
			if (columnfile != "") {
				if ((namefile == "") && (countfile == "")){ 
					namefile = current->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter.\n");  }
					else { 
						countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                        else { 
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the column format.\n");
                            abort = true; 
                        }	
					}	
				}
			}

			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
				
			calc = validParameter.valid(parameters, "calc");			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			util.splitAtDash(calc, Estimators);
			if (util.inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}

			string temp;
			temp = validParameter.valid(parameters, "precision");			if (temp == "not found") { temp = "100"; }
			util.mothurConvert(temp, precision); 
			
			temp = validParameter.valid(parameters, "cutoff");			if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
            temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
            
            temp = validParameter.valid(parameters, "iters");			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
            
            temp = validParameter.valid(parameters, "subsample");		if (temp == "not found") { temp = "F"; }
			if (util.isNumeric1(temp)) { util.mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (util.isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
            
            if (!subsample) { iters = 1; }
            
            temp = validParameter.valid(parameters, "withreplacement");		if (temp == "not found"){	temp = "f";		}
            withReplacement = util.isTrue(temp);
            
            if (subsample && (format != "sharedfile")) { m->mothurOut("[ERROR]: the subsample parameter can only be used with a shared file.\n"); abort=true; }
            
			if (outputdir == ""){ outputdir += util.hasPath(inputfile);  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "TreeSharedCommand", "TreeSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

TreeSharedCommand::~TreeSharedCommand(){}

//**********************************************************************************************************************

int TreeSharedCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (format == "sharedfile") {
			InputData input(sharedfile, "sharedfile", Groups);
            set<string> processedLabels;
            set<string> userLabels = labels;
            string lastLabel = "";
            
			SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
			
            if (subsample) {
                if (subsampleSize == -1) { //user has not set size, set size = smallest samples size
                    subsampleSize = lookup->getNumSeqsSmallestGroup();
                }else {
                    lookup->removeGroups(subsampleSize);
                    Groups = lookup->getNamesGroups();
                    Treenames = Groups;
                }
                
                if (lookup->size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command.\n");  m->setControl_pressed(true); return 0; }
            }
            numGroups = lookup->size();
            
            if (numGroups < 2) { m->mothurOut("[ERROR]: You have not provided enough valid groups.  I cannot run the command.\n");   return 0; }
			
			//create treemap class from groupmap for tree class to use
			CountTable ct;
            set<string> nameMap; map<string, string> groupMap; set<string> gps;
            for (int i = 0; i < Groups.size(); i++) {
                nameMap.insert(Groups[i]); gps.insert(Groups[i]); groupMap[Groups[i]] = Groups[i];
            }
            ct.createTable(nameMap, groupMap, gps);
			
			//fills tree names with shared files groups
			Treenames = lookup->getNamesGroups();
            
			if (m->getControl_pressed()) { return 0; }
            
            while (lookup != nullptr) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                createProcesses(lookup, ct); delete lookup;
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
			
			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  }  return 0; }
		}else{
			//read in dist file
			filename = inputfile;
            
            ReadMatrix* readMatrix;
			if (format == "column") { readMatrix = new ReadColumnMatrix(filename); }	
			else if (format == "phylip") { readMatrix = new ReadPhylipMatrix(filename); }
				
			readMatrix->setCutoff(cutoff);
	
            ListVector* list;
            if(namefile != ""){	
                NameAssignment* nameMap = new NameAssignment(namefile);
                nameMap->readMap();
                readMatrix->read(nameMap);
                list = readMatrix->getListVector();
                delete nameMap;
            }else if (countfile != "") {
                CountTable* ct = new CountTable();
                ct->readTable(countfile, true, false);
                readMatrix->read(ct);
                list = readMatrix->getListVector();
                delete ct;
            }else { NameAssignment* nameMap = nullptr; readMatrix->read(nameMap); list = readMatrix->getListVector(); }

			SparseDistanceMatrix* dMatrix = readMatrix->getDMatrix();
			Treenames.clear();
            
			//make treemap
			CountTable ct;
            set<string> nameMap;
            map<string, string> groupMap;
            set<string> gps;
            for (int i = 0; i < list->getNumBins(); i++) {
                string bin = list->get(i);
                nameMap.insert(bin); 
                gps.insert(bin); 
                groupMap[bin] = bin;
                Treenames.push_back(bin);
            }
            ct.createTable(nameMap, groupMap, gps);
			vector<string> namesGroups = ct.getNamesOfGroups();
			
			if (m->getControl_pressed()) { return 0; }
			
			vector< vector<double> > matrix = makeSimsDist(dMatrix, list->getNumBins());
            delete readMatrix; delete dMatrix;
			
			if (m->getControl_pressed()) { return 0; }

			//create a new filename
            map<string, string> variables; 
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
			string outputFile = getOutputFileName("tree",variables);	
			outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
				//printSims(cout, matrix, Treenames);
            Tree* newTree = new Tree(&ct, matrix, Treenames);
            if (m->getControl_pressed()) { delete newTree; newTree = nullptr; }
            else { newTree->assembleTree(); }
 
            if (newTree != nullptr) {  newTree->createNewickFile(outputFile);  delete newTree; }
			
			if (m->getControl_pressed()) { return 0; } m->mothurOut("Tree complete.\n");
		}
				
		//set tree file as new current treefile
		string currentName = "";
		itTypes = outputTypes.find("tree");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTreeFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeSharedCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
void TreeSharedCommand::printSims(ostream& out, vector< vector<double> >& simMatrix, vector<string> groupNames) {
    try {
        
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        out << simMatrix.size() << endl;
        for (int b = 0; b < simMatrix.size(); b++)	{
            out << groupNames[b];
            for (int n = 0; n < b; n++)	{
                out  << '\t' << simMatrix[b][n];
            }
            out << endl;
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "TreeSharedCommand", "printSims");
        exit(1);
    }
}
/***********************************************************/
vector< vector<double> > TreeSharedCommand::makeSimsDist(SparseDistanceMatrix* matrix, int numGroups) {
	try {
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
			
                    if (m->getControl_pressed()) { return simMatrix; }
                }
            }
		}

		return simMatrix;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeSharedCommand", "makeSimsDist");
		exit(1);
	}
}
/**************************************************************************************************/
int driverTreeShared(vector<SharedRAbundVector*>& thisLookup, vector< vector<seqDist> >& calcDists, vector<Calculator*> treeCalculators, MothurOut* m) {
    try {
        vector<SharedRAbundVector*> subset;
        
        for (int k = 0; k < thisLookup.size(); k++) { // pass cdd each set of groups to compare
            
            for (int l = 0; l < k; l++) {
                
                if (k != l) { //we dont need to similarity of a groups to itself
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
                        
                        if (m->getControl_pressed()) { return 1; }
                        
                        seqDist temp(l, k, tempdata[0]);
                        calcDists[i].push_back(temp);
                    }
                }
            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TreeSharedCommand", "driverTreeShared");
        exit(1);
    }
}
/**************************************************************************************************/
struct treeSharedData {
    SharedRAbundVectors* thisLookup;
    vector< vector< vector<seqDist> > > calcDistsTotals;  //each iter, one for each calc, then each groupCombos dists. this will be used to make .dist files
    vector< vector< vector<double> > > matrices; //for each calculator a square matrix to represent the distances, only filled by main thread
    vector<string>  Estimators;
    long long numIters;
    MothurOut* m;
    int count, subsampleSize;
    bool subsample, withReplacement;
    
    treeSharedData(){}
    treeSharedData(long long st, bool su, bool wr, int subsize, vector<string> est, SharedRAbundVectors* lu) {
        m = MothurOut::getInstance();
        numIters = st;
        Estimators = est;
        thisLookup = lu;
        count = 0;
        subsample = su;
        withReplacement = wr;
        subsampleSize = subsize;
    }
};
/***********************************************************/
int process(treeSharedData* params) {
    try{
        
        ValidCalculators validCalculator;
        vector<Calculator*> treeCalculators;
        for (int i=0; i<params->Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("treegroup", params->Estimators[i]) ) {
                if (params->Estimators[i] == "sharedsobs") {
                    treeCalculators.push_back(new SharedSobsCS());
                }else if (params->Estimators[i] == "sharedchao") {
                    treeCalculators.push_back(new SharedChao1());
                }else if (params->Estimators[i] == "sharedace") {
                    treeCalculators.push_back(new SharedAce());
                }else if (params->Estimators[i] == "jabund") {
                    treeCalculators.push_back(new JAbund());
                }else if (params->Estimators[i] == "sorabund") {
                    treeCalculators.push_back(new SorAbund());
                }else if (params->Estimators[i] == "jclass") {
                    treeCalculators.push_back(new Jclass());
                }else if (params->Estimators[i] == "sorclass") {
                    treeCalculators.push_back(new SorClass());
                }else if (params->Estimators[i] == "jest") {
                    treeCalculators.push_back(new Jest());
                }else if (params->Estimators[i] == "sorest") {
                    treeCalculators.push_back(new SorEst());
                }else if (params->Estimators[i] == "thetayc") {
                    treeCalculators.push_back(new ThetaYC());
                }else if (params->Estimators[i] == "thetan") {
                    treeCalculators.push_back(new ThetaN());
                }else if (params->Estimators[i] == "kstest") {
                    treeCalculators.push_back(new KSTest());
                }else if (params->Estimators[i] == "sharednseqs") {
                    treeCalculators.push_back(new SharedNSeqs());
                }else if (params->Estimators[i] == "ochiai") {
                    treeCalculators.push_back(new Ochiai());
                }else if (params->Estimators[i] == "anderberg") {
                    treeCalculators.push_back(new Anderberg());
                }else if (params->Estimators[i] == "kulczynski") {
                    treeCalculators.push_back(new Kulczynski());
                }else if (params->Estimators[i] == "kulczynskicody") {
                    treeCalculators.push_back(new KulczynskiCody());
                }else if (params->Estimators[i] == "lennon") {
                    treeCalculators.push_back(new Lennon());
                }else if (params->Estimators[i] == "morisitahorn") {
                    treeCalculators.push_back(new MorHorn());
                }else if (params->Estimators[i] == "braycurtis") {
                    treeCalculators.push_back(new BrayCurtis());
                }else if (params->Estimators[i] == "whittaker") {
                    treeCalculators.push_back(new Whittaker());
                }else if (params->Estimators[i] == "odum") {
                    treeCalculators.push_back(new Odum());
                }else if (params->Estimators[i] == "canberra") {
                    treeCalculators.push_back(new Canberra());
                }else if (params->Estimators[i] == "structeuclidean") {
                    treeCalculators.push_back(new StructEuclidean());
                }else if (params->Estimators[i] == "structchord") {
                    treeCalculators.push_back(new StructChord());
                }else if (params->Estimators[i] == "hellinger") {
                    treeCalculators.push_back(new Hellinger());
                }else if (params->Estimators[i] == "manhattan") {
                    treeCalculators.push_back(new Manhattan());
                }else if (params->Estimators[i] == "structpearson") {
                    treeCalculators.push_back(new StructPearson());
                }else if (params->Estimators[i] == "soergel") {
                    treeCalculators.push_back(new Soergel());
                }else if (params->Estimators[i] == "spearman") {
                    treeCalculators.push_back(new Spearman());
                }else if (params->Estimators[i] == "structkulczynski") {
                    treeCalculators.push_back(new StructKulczynski());
                }else if (params->Estimators[i] == "speciesprofile") {
                    treeCalculators.push_back(new SpeciesProfile());
                }else if (params->Estimators[i] == "hamming") {
                    treeCalculators.push_back(new Hamming());
                }else if (params->Estimators[i] == "structchi2") {
                    treeCalculators.push_back(new StructChi2());
                }else if (params->Estimators[i] == "gower") {
                    treeCalculators.push_back(new Gower());
                }else if (params->Estimators[i] == "memchi2") {
                    treeCalculators.push_back(new MemChi2());
                }else if (params->Estimators[i] == "memchord") {
                    treeCalculators.push_back(new MemChord());
                }else if (params->Estimators[i] == "memeuclidean") {
                    treeCalculators.push_back(new MemEuclidean());
                }else if (params->Estimators[i] == "mempearson") {
                    treeCalculators.push_back(new MemPearson());
                }else if (params->Estimators[i] == "jsd") {
                    treeCalculators.push_back(new JSD());
                }else if (params->Estimators[i] == "rjsd") {
                    treeCalculators.push_back(new RJSD());
                }
            }
        }
        
        //if the users entered no valid calculators don't execute command
        if (treeCalculators.size() == 0) { params->m->mothurOut("You have given no valid calculators.\n");  return 0; }
        
        params->Estimators.clear();
        for (int i=0; i<treeCalculators.size(); i++) { params->Estimators.push_back(treeCalculators[i]->getName()); }
        
        vector< vector<seqDist>  > calcDists; calcDists.resize(treeCalculators.size());
        SubSample sample; 
        for (int thisIter = 0; thisIter < params->numIters; thisIter++) {
            
            SharedRAbundVectors* thisItersLookup = new SharedRAbundVectors(*params->thisLookup);
            vector<string> namesOfGroups = thisItersLookup->getNamesGroups();
            
            if (params->subsample) {
                if (params->withReplacement)    { sample.getSampleWithReplacement(thisItersLookup, params->subsampleSize);  }
                else                            { sample.getSample(thisItersLookup, params->subsampleSize);                 }
            }
            
            vector<SharedRAbundVector*> thisItersRabunds = thisItersLookup->getSharedRAbundVectors();
            vector<string> thisItersGroupNames = params->thisLookup->getNamesGroups();
            
            driverTreeShared(thisItersRabunds, calcDists, treeCalculators, params->m);
            
            for (int i = 0; i < thisItersRabunds.size(); i++) { delete thisItersRabunds[i]; }
            
            if (params->subsample){
                if((thisIter+1) % 100 == 0){	params->m->mothurOutJustToScreen(toString(thisIter+1)+"\n"); 		}
                params->calcDistsTotals.push_back(calcDists);
                for (int i = 0; i < calcDists.size(); i++) {
                    for (int j = 0; j < calcDists[i].size(); j++) {
                        if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: Results: iter = " + toString(thisIter) + ", " + thisItersGroupNames[calcDists[i][j].seq1] + " - " + thisItersGroupNames[calcDists[i][j].seq2] + " distance = " + toString(calcDists[i][j].dist) + ".\n");  }
                    }
                }
            }else { //print results for whole dataset
                for (int i = 0; i < calcDists.size(); i++) {
                    if (params->m->getControl_pressed()) { break; }
                    
                    //initialize matrix
                    vector< vector<double> > matrix; //square matrix to represent the distance
                    matrix.resize(thisItersLookup->size());
                    for (int k = 0; k < thisItersLookup->size(); k++) {  matrix[k].resize(thisItersLookup->size(), 0.0); }
                    
                    for (int j = 0; j < calcDists[i].size(); j++) {
                        int row = calcDists[i][j].seq1;
                        int column = calcDists[i][j].seq2;
                        double dist = calcDists[i][j].dist;
                        
                        matrix[row][column] = -(dist-1.0);
                        matrix[column][row] = -(dist-1.0);
                    }
                    params->matrices.push_back(matrix);
                }
            }
            for (int i = 0; i < calcDists.size(); i++) {  calcDists[i].clear(); }
            delete thisItersLookup;
        }
        if((params->numIters) % 100 != 0){	params->m->mothurOutJustToScreen(toString(params->numIters)+"\n"); 		}
        for (int i=0; i<treeCalculators.size(); i++) { delete treeCalculators[i]; }
        
        return 0;
    }
    catch(exception& e) {
        params->m->errorOut(e, "TreeSharedCommand", "process");
        exit(1);
    }
}
/***********************************************************/
int TreeSharedCommand::createProcesses(SharedRAbundVectors*& thisLookup, CountTable& ct){
    try {
        
        vector<string> groupNames = thisLookup->getNamesGroups();
        Treenames = groupNames; //may have changed if subsample eliminated groups
        
        vector<int> lines;
        if (processors > (iters+1)) { processors = iters+1; }
        
        //figure out how many sequences you have to process
        int numItersPerProcessor = (iters+1) / processors;
        for (int i = 0; i < processors; i++) {
            if(i == (processors - 1)){	numItersPerProcessor = (iters+1) - i * numItersPerProcessor; 	}
            lines.push_back(numItersPerProcessor);
        }
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<treeSharedData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            
            //make copy of lookup so we don't get access violations
            SharedRAbundVectors* newLookup = new SharedRAbundVectors(*thisLookup);
            treeSharedData* dataBundle = new treeSharedData(lines[i+1], subsample, withReplacement, subsampleSize, Estimators, newLookup);
            
            data.push_back(dataBundle);
            workerThreads.push_back(new std::thread(process, dataBundle));
        }
        
        //make copy of lookup so we don't get access violations
        SharedRAbundVectors* newLookup = new SharedRAbundVectors(*thisLookup);
        treeSharedData* dataBundle = new treeSharedData(lines[0], subsample, withReplacement, subsampleSize, Estimators, newLookup);
        process(dataBundle);
        delete newLookup;
        
        Estimators.clear(); Estimators = dataBundle->Estimators;
        vector< vector< vector<seqDist> > > calcDistsTotals = dataBundle->calcDistsTotals;
        vector< vector< vector<double> > > matrices = dataBundle->matrices;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            //get calcDistsTotal info - one entry per iter
            for (int j = 0; j < data[i]->calcDistsTotals.size(); j++) { calcDistsTotals.push_back(data[i]->calcDistsTotals[j]); }
            
            delete data[i]->thisLookup;
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        
        if (subsample) {
            //we need to find the average distance and standard deviation for each groups distance
            vector< vector<seqDist>  > calcAverages = util.getAverages(calcDistsTotals);
            
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: found averages.\n"); }
            
            //create average tree for each calc
            for (int i = 0; i < Estimators.size(); i++) {
                vector< vector<double> > matrix; //square matrix to represent the distance
                matrix.resize(thisLookup->size());
                for (int k = 0; k < thisLookup->size(); k++) {  matrix[k].resize(thisLookup->size(), 0.0); }
                
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    int row = calcAverages[i][j].seq1;
                    int column = calcAverages[i][j].seq2;
                    float dist = calcAverages[i][j].dist;
                    
                    matrix[row][column] = -(dist-1.0); //-(matrix->seqVec[i][j].dist -1.0)
                    matrix[column][row] = -(dist-1.0);
                }
                //printSims(cout, matrix, Treenames);
                //create a new filename
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
                variables["[calc]"] = Estimators[i];
                variables["[distance]"] = thisLookup->getLabel();
                variables["[tag]"] = "ave";
                string outputFile = getOutputFileName("tree",variables);
                outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
                
                //creates tree from similarity matrix and write out file
                Tree* newTree = new Tree(&ct, matrix, Treenames);
                if (m->getControl_pressed()) { delete newTree; newTree = nullptr; }
                else { newTree->assembleTree(); }
                if (newTree != nullptr) { newTree->createNewickFile(outputFile);  delete newTree; }
            }
            
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: done averages trees.\n"); }
            
            //create all trees for each calc and find their consensus tree
            for (int i = 0; i < Estimators.size(); i++) {
                if (m->getControl_pressed()) { break; }
                
                //create a new filename
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
                variables["[calc]"] = Estimators[i];
                variables["[distance]"] = thisLookup->getLabel();
                variables["[tag]"] = "all";
                string outputFile = getOutputFileName("tree",variables);
                outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
                
                ofstream outAll;
                util.openOutputFile(outputFile, outAll);
                
                vector<Tree*> trees;
                for (int myIter = 0; myIter < iters; myIter++) {
                    
                    if(m->getControl_pressed()) { break; }
                    
                    //initialize matrix
                    vector< vector<double> > matrix; //square matrix to represent the distance
                    matrix.resize(thisLookup->size());
                    for (int k = 0; k < thisLookup->size(); k++) {  matrix[k].resize(thisLookup->size(), 0.0); }
                    
                    for (int j = 0; j < calcDistsTotals[myIter][i].size(); j++) {
                        int row = calcDistsTotals[myIter][i][j].seq1;
                        int column = calcDistsTotals[myIter][i][j].seq2;
                        double dist = calcDistsTotals[myIter][i][j].dist;
                        
                        matrix[row][column] = -(dist-1.0);
                        matrix[column][row] = -(dist-1.0);
                    }
                    
                    //creates tree from similarity matrix and write out file
                    Tree* newTree = new Tree(&ct, matrix, Treenames);
                    if (m->getControl_pressed()) { delete newTree; newTree = nullptr; }
                    else { newTree->assembleTree(); }
                    if (newTree != nullptr) {
                        newTree->print(outAll);
                        trees.push_back(newTree);
                    }
                }
                outAll.close();
                if (m->getControl_pressed()) { for (int k = 0; k < trees.size(); k++) { delete trees[k]; } }
                
                if (m->getDebug()) {  m->mothurOut("[DEBUG]: done all trees.\n"); }
                
                Consensus consensus;
                Tree* conTree = consensus.getTree(trees);
                
                if (m->getDebug()) {  m->mothurOut("[DEBUG]: done cons tree.\n"); }
                
                //create a new filename
                variables["[tag]"] = "cons";
                string conFile = getOutputFileName("tree",variables);
            
                outputNames.push_back(conFile); outputTypes["tree"].push_back(conFile);
                ofstream outTree;
                util.openOutputFile(conFile, outTree);
                
                if (conTree != nullptr) { conTree->print(outTree, "boot"); delete conTree; }
            }
        }else {
            for (int i = 0; i < matrices.size(); i++) {
                if (m->getControl_pressed()) { break; }
                
                //initialize matrix
                vector< vector<double> > matrix = matrices[i]; //square matrix to represent the distance
                
                //create a new filename
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
                variables["[calc]"] = Estimators[i];
                variables["[distance]"] = thisLookup->getLabel();
                variables["[tag]"] = "";
                string outputFile = getOutputFileName("tree",variables);
                outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
                
                //creates tree from similarity matrix and write out file
                Tree* newTree = new Tree(&ct, matrix, Treenames);
                if (m->getControl_pressed()) { delete newTree; newTree = nullptr; }
                else { newTree->assembleTree(); }
                if (newTree != nullptr) { newTree->createNewickFile(outputFile);  delete newTree; }
            }
        }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "TreeSharedCommand", "createProcesses");
        exit(1);
    }
}
/***********************************************************/

	

