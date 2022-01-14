/*
 *  unifracweightedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unifracweightedcommand.h"
#include "consensus.h"
#include "subsample.h"
#include "treereader.h"

//**********************************************************************************************************************
vector<string> UnifracWeightedCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none","weighted-wsummary",false,true,true); parameters.push_back(ptree);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","",false,false); parameters.push_back(psubsample);
        CommandParameter pwithreplacement("withreplacement", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pwithreplacement);
        CommandParameter pconsensus("consensus", "Boolean", "", "F", "", "", "","tree",false,false); parameters.push_back(pconsensus);
        CommandParameter prandom("random", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prandom);
		CommandParameter pdistance("distance", "Multiple", "column-lt-square-phylip", "column", "", "", "","phylip-column",false,false); parameters.push_back(pdistance);
		CommandParameter proot("root", "Boolean", "F", "", "", "", "","",false,false); parameters.push_back(proot);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> tempOutNames;
        outputTypes["weighted"] = tempOutNames;
        outputTypes["wsummary"] = tempOutNames;
        outputTypes["phylip"] = tempOutNames;
        outputTypes["column"] = tempOutNames;
        outputTypes["tree"] = tempOutNames;
        
        abort = false; calledHelp = false;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string UnifracWeightedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The unifrac.weighted command parameters are tree, group, name, count, groups, iters, distance, processors, root, subsample, consensus and random.  tree parameter is required unless you have valid current tree file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		helpString += "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n";
		helpString += "The distance parameter allows you to create a distance file from the results. The default is false.\n";
		helpString += "The random parameter allows you to shut off the comparison to random trees. The default is false, meaning don't compare your trees with randomly generated trees.\n";
		helpString += "The root parameter allows you to include the entire root in your calculations. The default is false, meaning stop at the root for this comparision instead of the root of the entire tree.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
        helpString += "The subsample parameter allows you to enter the size pergroup of the sample or you can set subsample=T and mothur will use the size of your smallest group. The subsample parameter may only be used with a group file.\n";
        helpString += "The withreplacement parameter allows you to indicate you want to subsample your data allowing for the same read to be included multiple times. Default=f. \n";
        helpString += "The consensus parameter allows you to indicate you would like trees built from distance matrices created with the results, as well as a consensus tree built from these trees. Default=F.\n";
        helpString += "The unifrac.weighted command should be in the following format: unifrac.weighted(groups=yourGroups, iters=yourIters).\n";
		helpString += "Example unifrac.weighted(groups=A-B-C, iters=500).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and iters is 1000.\n";
		helpString += "The unifrac.weighted command output two files: .weighted and .wsummary their descriptions are in the manual.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string UnifracWeightedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        if (type == "weighted")            {  pattern = "[filename],weighted-[filename],[tag],weighted";   }
        else if (type == "wsummary")        {  pattern = "[filename],[tag],wsummary";   }
        else if (type == "phylip")           {  pattern = "[filename],[tag],[tag2],dist";   }
        else if (type == "column")           {  pattern = "[filename],[tag],[tag2],dist";   }
        else if (type == "tree")             {  pattern = "[filename],[tag],[tag2],tre";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "UnifracWeightedCommand", "getOutputPattern");
        exit(1);
    }
}
/***********************************************************/
UnifracWeightedCommand::UnifracWeightedCommand(string option) : Command() {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters=parser.getParameters();
			
			ValidParameters validParameter;
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { treefile = ""; abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = current->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter.\n");  }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required.\n");  abort = true; }								
			}else { current->setTreeFile(treefile); }	
			
			//check for required parameters
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
			
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n");  abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count.\n");  abort=true;
            }

					if (outputdir == ""){    outputdir = util.hasPath(treefile);	}
			
																	
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
				
			itersString = validParameter.valid(parameters, "iters");			if (itersString == "not found") { itersString = "1000"; }
			util.mothurConvert(itersString, iters); 
			
			string temp = validParameter.valid(parameters, "distance");
			if (temp == "not found") { phylip = false; outputForm = ""; }
			else{
                if (temp=="phylip") { temp = "lt"; }
				if ((temp == "lt") || (temp == "column") || (temp == "square")) {  phylip = true;  outputForm = temp; }
				else { m->mothurOut("Options for distance are: lt, square, or column. Using lt.\n");  phylip = true; outputForm = "lt"; }
			}
			
			temp = validParameter.valid(parameters, "random");				if (temp == "not found") { temp = "F"; }
			random = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "root");					if (temp == "not found") { temp = "F"; }
			includeRoot = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
            processors = current->setProcessors(temp);
            
            temp = validParameter.valid(parameters, "subsample");		if (temp == "not found") { temp = "F"; }
			if (util.isNumeric1(temp)) { util.mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (util.isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
			
            if (!subsample) { subsampleIters = 0;   }
            else { subsampleIters = iters;          }
            
            temp = validParameter.valid(parameters, "withreplacement");		if (temp == "not found"){	temp = "f";		}
            withReplacement = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "consensus");					if (temp == "not found") { temp = "F"; }
			consensus = util.isTrue(temp);
            
			if (subsample && random) {  m->mothurOut("[ERROR]: random must be false, if subsample=t.\n"); abort=true;  } 
			if (countfile == "") { if (subsample && (groupfile == "")) {  m->mothurOut("[ERROR]: if subsample=t, a group file must be provided.\n"); abort=true;  } }
            else {  
                CountTable testCt; 
                if ((!testCt.testGroups(countfile)) && (subsample)) {
                    m->mothurOut("[ERROR]: if subsample=t, a count file with group info must be provided.\n"); abort=true;  
                }
            }
            if (subsample && (!phylip)) { phylip=true; outputForm = "lt"; }
            if (consensus && (!subsample)) { m->mothurOut("[ERROR]: you cannot use consensus without subsample.\n"); abort=true; }
            
		}
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "UnifracWeightedCommand");
		exit(1);
	}
}
/***********************************************************/
int UnifracWeightedCommand::execute() {
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        long start = time(NULL);
        
        TreeReader* reader;
        if (countfile == "") { reader = new TreeReader(treefile, groupfile, namefile); }
        else { reader = new TreeReader(treefile, countfile); }
        vector<Tree*> T = reader->getTrees();
        CountTable* ct; ct = T[0]->getCountTable();
        if ((Groups.size() == 0) || (Groups.size() < 2)) {  Groups = ct->getNamesOfGroups();  } //must have at least 2 groups to compare
        delete reader;
       
        if (m->getControl_pressed()) {  delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } return 0; }
		
        map<string, string> variables;
		vector<string> nameGroups = ct->getNamesOfGroups();
        if (Groups.size() < 2) { m->mothurOut("[ERROR]: You cannot run unifrac.weighted with less than 2 groups, aborting.\n"); delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } return 0; }
        if (m->getControl_pressed()) {  delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } return 0; }
        
        if ((Groups.size() == 0) || (Groups.size() < 2)) {  Groups = ct->getNamesOfGroups();  } //must have at least 2 groups to compare
   
        //set or check size
        if (subsample) {
            //user has not set size, set size = smallest samples size
            if (subsampleSize == -1) {
                subsampleSize = ct->getNumSeqsSmallestGroup();
                m->mothurOut("\nSetting subsample size to " + toString(subsampleSize) + ".\n\n");
            }else { //eliminate any too small groups
                vector<string> newGroups = Groups;
                Groups.clear();
                for (int i = 0; i < newGroups.size(); i++) {
                    int thisSize = ct->getGroupCount(newGroups[i]);
                    
                    if (thisSize >= subsampleSize) {    Groups.push_back(newGroups[i]);	}
                    else {   m->mothurOut("You have selected a size that is larger than "+newGroups[i]+" number of sequences, removing "+newGroups[i]+".\n"); }
                }
            }
        }
        
        vector<string> groupComb; util.getCombos(groupComb, Groups, numComp); //here in case some groups are removed by subsample
        
        if (numComp < processors) { processors = numComp; m->mothurOut("Reducing processors to " + toString(numComp) + ".\n"); }
        if (consensus && (numComp < 2)) { m->mothurOut("consensus can only be used with numComparisions greater than 1, setting consensus=f.\n"); consensus=false; }
        
        Weighted weighted(includeRoot, Groups);
        
        for (int i = 0; i < T.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            vector<double> WScoreSig;  //tree weighted score signifigance when compared to random trees - percentage of random trees with that score or lower.
            vector< vector<double> > uScores;  uScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
            vector<double> userData; userData.resize(numComp,0);  //weighted score info for user tree. data[0] = weightedscore AB, data[1] = weightedscore AC...
            vector<double> randomData; randomData.resize(numComp,0); //weighted score info for random trees. data[0] = weightedscore AB, data[1] = weightedscore AC...
            
            userData = weighted.getValues(T[i], processors, outputdir); //userData[0] = weightedscore
            if (m->getControl_pressed()) { break; }
            
            if (phylip) {	createPhylipFile((i+1), userData);          }
            if (random) { runRandomCalcs(T[i], ct, userData, (i+1), WScoreSig, groupComb);    }
            printWSummaryFile((i+1), userData, WScoreSig, groupComb);
            
            if (m->getControl_pressed()) { break; }
            
            //subsample loop
            vector< vector<double> > calcDistsTotals;  //each iter, each groupCombos dists. this will be used to make .dist files
            SubSample sample;
            for (int thisIter = 0; thisIter < subsampleIters; thisIter++) { //subsampleIters=0, if subsample=f.
                if (m->getControl_pressed()) { break; }

                //uses method of setting groups to doNotIncludeMe
                //copy to preserve old one - would do this in subsample but memory cleanup becomes messy.
                CountTable* newCt = new CountTable();
                int sampleTime = time(NULL);
                
                Tree* subSampleTree;
                if (withReplacement)    { subSampleTree = sample.getSampleWithReplacement(T[i], ct, newCt, subsampleSize, Groups);  }
                else                    { subSampleTree = sample.getSample(T[i], ct, newCt, subsampleSize, Groups);                 }
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: iter " + toString(thisIter) + " took " + toString(time(NULL) - sampleTime) + " seconds to sample tree.\n"); }
                
                //call new weighted function
                vector<double> iterData; iterData.resize(numComp,0);
                Weighted thisWeighted(includeRoot, Groups);
                iterData = thisWeighted.getValues(subSampleTree, processors, outputdir); //userData[0] = weightedscore
                
                //save data to make ave dist, std dist
                calcDistsTotals.push_back(iterData);
                
                delete newCt; delete subSampleTree;
                
                if((thisIter+1) % 100 == 0){	m->mothurOutJustToScreen(toString(thisIter+1)+"\n"); 	}
            }
            
            if (m->getControl_pressed()) { break; }
            
            if (subsample) {  getAverageSTDMatrices(calcDistsTotals, i);    }
            if (consensus) {  getConsensusTrees(calcDistsTotals, i);        }
        }
		delete ct;  for (int i = 0; i < T.size(); i++) { delete T[i]; }
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } return 0; }
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run unifrac.weighted.\n"); 
		
		//set phylip file as new current phylipfile
		string currentName = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setPhylipFile(currentName); }
		}
		
		//set column file as new current columnfile
		itTypes = outputTypes.find("column");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setColumnFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n");
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int UnifracWeightedCommand::getAverageSTDMatrices(vector< vector<double> >& dists, int treeNum) {
	try {
        vector<double> averages = util.getAverages(dists);
        vector<double> stdDev = util.getStandardDeviation(dists, averages);
        
        int numGroups = Groups.size();
        vector< vector<double> > avedists;
        for (int i = 0; i < numGroups; i++) {
            vector<double> temp; temp.resize(numGroups, 0.0);
            avedists.push_back(temp);
        }
        
        //make matrix with scores in it
        vector< vector<double> > stddists;	//stddists.resize(m->getNumGroups());
        for (int i = 0; i < numGroups; i++) {
            vector<double> temp; for (int j = 0; j < numGroups; j++) { temp.push_back(0.0); }
            stddists.push_back(temp);
        }
        
        //flip it so you can print it
        int count = 0;
        for (int r=0; r< numGroups; r++) {
            for (int l = 0; l < r; l++) {
                avedists[r][l] = averages[count];
                avedists[l][r] = averages[count];
                stddists[r][l] = stdDev[count];
                stddists[l][r] = stdDev[count];
                count++;
            }
        }
        
        map<string, string> variables; 
		variables["[filename]"] = outputdir + util.getSimpleName(treefile);
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "weighted.ave";
        string aveFileName = getOutputFileName("phylip",variables);
        if (outputForm != "column") { outputNames.push_back(aveFileName); outputTypes["phylip"].push_back(aveFileName);  }
        else { outputNames.push_back(aveFileName); outputTypes["column"].push_back(aveFileName);  }
        ofstream out; util.openOutputFile(aveFileName, out);
        
        variables["[tag2]"] = "weighted.std";
        string stdFileName = getOutputFileName("phylip",variables);
        if (outputForm != "column") { outputNames.push_back(stdFileName); outputTypes["phylip"].push_back(stdFileName); }
        else { outputNames.push_back(stdFileName); outputTypes["column"].push_back(stdFileName); }        
        ofstream outStd; util.openOutputFile(stdFileName, outStd);
        
        if ((outputForm == "lt") || (outputForm == "square")) {
            //output numSeqs
            out << numGroups << endl;
            outStd << numGroups << endl;
        }
        
        //output to file
        for (int r=0; r< numGroups; r++) {
            string name = Groups[r];
            if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } } //pad with spaces to make compatible
            
            if (outputForm == "lt") {
                out << name; outStd << name;
                for (int l = 0; l < r; l++) {	out  << '\t' << avedists[r][l];  outStd   << '\t' << stddists[r][l];}  //output distances
                out << endl;  outStd << endl;
            }else if (outputForm == "square") {
                out << name; outStd << name;
                for (int l = 0; l < numGroups; l++) {	out  << '\t' << avedists[r][l]; outStd  << '\t' << stddists[r][l]; }  //output distances
                out << endl; outStd << endl;
            }else{
                for (int l = 0; l < r; l++) {	
                    string otherName = Groups[l];
                    if (otherName.length() < 10) { while (otherName.length() < 10) {  otherName += " ";  } }  //pad with spaces to make compatible
                    
                    out  << name << '\t' << otherName  << '\t' << avedists[r][l] << endl;  //output distances
                    outStd  << name << '\t' << otherName  << '\t' << stddists[r][l] << endl;
                }
            }
        }
        out.close(); outStd.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "getAverageSTDMatrices");
		exit(1);
	}
}
/**************************************************************************************************/
int UnifracWeightedCommand::getConsensusTrees(vector< vector<double> >& dists, int treeNum) {
	try {
        ///create treemap class from groupmap for tree class to use
        CountTable newCt;
        set<string> nameMap;
        map<string, string> groupMap;
        set<string> gps;
        int numGroups = Groups.size();
        for (int i = 0; i < numGroups; i++) {
            nameMap.insert(Groups[i]);
            gps.insert(Groups[i]);
            groupMap[Groups[i]] = Groups[i];
        }
        newCt.createTable(nameMap, groupMap, gps);
        
        vector<Tree*> newTrees = buildTrees(dists, treeNum, newCt); //also creates .all.tre file containing the trees created

        if (m->getControl_pressed()) { return 0; }
        
        Consensus con;
        Tree* conTree = con.getTree(newTrees);
        
        //create a new filename
        map<string, string> variables; 
		variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "weighted.cons";
        string conFile = getOutputFileName("tree",variables);							
        outputNames.push_back(conFile); outputTypes["tree"].push_back(conFile); 
        ofstream outTree;
        util.openOutputFile(conFile, outTree);
        
        if (conTree != NULL) { conTree->print(outTree, "boot"); delete conTree; } outTree.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "getConsensusTrees");
		exit(1);
	}
}
/**************************************************************************************************/
vector<Tree*> UnifracWeightedCommand::buildTrees(vector< vector<double> >& dists, int treeNum, CountTable& myct) {
	try {
        vector<Tree*> trees;
        
        //create a new filename
        map<string, string> variables; 
		variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "weighted.all";
        string outputFile = getOutputFileName("tree",variables);				
        outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
        
        ofstream outAll;
        util.openOutputFile(outputFile, outAll);
        int numGroups = Groups.size();

        for (int i = 0; i < dists.size(); i++) { //dists[0] are the dists for the first subsampled tree.
            
            if (m->getControl_pressed()) { break; }
            
            //make matrix with scores in it
            vector< vector<double> > sims;	sims.resize(numGroups);
            for (int j = 0; j < numGroups; j++) { sims[j].resize(numGroups, 0.0); }
            
            int count = 0;
			for (int r=0; r<numGroups; r++) {
				for (int l = 0; l < r; l++) {
                    double sim = -(dists[i][count]-1.0);
					sims[r][l] = sim;
					sims[l][r] = sim;
					count++;
				}
			}

            //create tree
            Tree* tempTree = new Tree(&myct, sims, Groups);
            tempTree->assembleTree();
            trees.push_back(tempTree);
            tempTree->print(outAll); //print tree
        }
        outAll.close();
        
        if (m->getControl_pressed()) {  for (int i = 0; i < trees.size(); i++) {  delete trees[i]; trees[i] = NULL; } util.mothurRemove(outputFile); }
        
        return trees;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "buildTrees");
		exit(1);
	}
}
/**************************************************************************************************/
int UnifracWeightedCommand::runRandomCalcs(Tree* thisTree, CountTable* ct, vector<double> usersScores, int iter, vector<double>& WScoreSig, vector<string> groupComb) {
	try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getSimpleName(treefile);
        variables["[tag]"] = toString(iter);
        string wFileName = getOutputFileName("weighted", variables); 
        ColumnFile output(wFileName, itersString); ofstream out; util.openOutputFile(wFileName, out); out.close();
        outputNames.push_back(wFileName); outputTypes["weighted"].push_back(wFileName);
        
        //calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
        vector< vector<string> > namesOfGroupCombos;
        int numGroups = Groups.size();
        for (int a=0; a<numGroups; a++) {
            for (int l = 0; l < a; l++) {
                vector<string> groups; groups.push_back(Groups[a]); groups.push_back(Groups[l]);
                namesOfGroupCombos.push_back(groups);
            }
        }
        
        vector<vector<int> > randomTreeNodes;
        for (int f = 0; f < numComp; f++) { randomTreeNodes.push_back(thisTree->getNodes(namesOfGroupCombos[f])); }
        vector<vector<int> > savedRandomTreeNodes = randomTreeNodes;
        
        //get scores for random trees
        vector<vector<double> > rScores; rScores.resize(numComp);
        for (int i = 0; i < iters; i++) {
            if (m->getControl_pressed()) { return 0; }
            randomTreeNodes = savedRandomTreeNodes;
            
            for (int f = 0; f < numComp; f++) {   util.mothurRandomShuffle(randomTreeNodes[f]);   }
            
            vector<double> thisItersRScores = createProcesses(thisTree, ct, namesOfGroupCombos, randomTreeNodes);
            
            for (int f = 0; f < numComp; f++) {   rScores[f].push_back(thisItersRScores[f]);  }
            
            if((i+1) % 100 == 0){	m->mothurOut(toString(i+1)+"\n");		}
        }
        
        //find the signifigance of the score for summary file
        for (int f = 0; f < numComp; f++) {
            //sort random scores
            sort(rScores[f].begin(), rScores[f].end());
            
            //the index of the score higher than yours is returned 
            //so if you have 1000 random trees the index returned is 100 
            //then there are 900 trees with a score greater then you. 
            //giving you a signifigance of 0.900
            int index = findIndex(usersScores[f], f, rScores);    if (index == -1) { m->mothurOut("error in UnifracWeightedCommand\n");  exit(1); } //error code
            
            //the signifigance is the number of trees with the users score or higher 
            WScoreSig.push_back((iters-index)/(float)iters);
        }
        
        set<double>  validScores;  //map contains scores from random
        vector< map<double, double> > rScoreFreq;  //map <weighted score, number of random trees with that score.> -vector entry for each combination.
        vector< map<double, double> > rCumul;  //map <weighted score, cumulative percentage of number of random trees with that score or higher.> -vector entry for each c
        calculateFreqsCumuls(validScores, rScores, rScoreFreq, rCumul);
        
        vector<string> tags; tags.push_back("Score"); tags.push_back("RandFreq"); tags.push_back("RandCumul");
        for(int a = 0; a < numComp; a++) {
            output.setLabelName(groupComb[a], tags);
            //print each line
            for (set<double>::iterator it = validScores.begin(); it != validScores.end(); it++) {
                vector<double> data; data.push_back(*it);  data.push_back(rScoreFreq[a][*it]); data.push_back(rCumul[a][*it]);
                output.updateOutput(data);
            } 
            output.resetFile();
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "runRandomCalcs");
		exit(1);
	}
}
/***********************************************************************/
struct weightedRandomData {
    bool includeRoot;
    int count, numComps, start, num;
    vector<string> Groups, Treenames;
    vector<double> scores;
    vector< vector<string> > namesOfGroupCombos;
    vector<vector<int> > randomizedTreeNodes;
    MothurOut* m;
    Tree* t;
    CountTable* ct;
    Utils util;
    
    weightedRandomData(){}
    weightedRandomData(int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count, bool ir, vector<string> g, vector<vector<int> > randomTreeNodes) {
        m = MothurOut::getInstance();
        num = en;
        start = st;
        namesOfGroupCombos = ngc;
        numComps = namesOfGroupCombos.size();
        randomizedTreeNodes = randomTreeNodes;
        t = tree;
        ct = count;
        includeRoot = ir;
        Groups = g;
        Treenames = t->getTreeNames();
        count = 0;
    }
};
/**************************************************************************************************/
void driverWeightedRandom(weightedRandomData* params) {
    try {
        Weighted weighted(params->includeRoot, params->Groups);
        
        params->count = 0;
        
        Tree* randT = new Tree(params->ct, params->Treenames);

        for (int h = params->start; h < (params->start+params->num); h++) {
            
            if (params->m->getControl_pressed()) { break; }
            
            string groupA = params->namesOfGroupCombos[h][0];
            string groupB = params->namesOfGroupCombos[h][1];
            vector<int> treeNodesFromTheseGroups = params->randomizedTreeNodes[h];
                       
            //copy T[i]'s info.
            randT->getCopy(params->t);
            
            //create a random tree with same topology as T[i], but different labels
            randT->assembleRandomUnifracTree(params->randomizedTreeNodes[h]);
            
            if (params->m->getControl_pressed()) { break; }
            
            //get wscore of random tree
            EstOutput randomData = weighted.getValues(randT, groupA, groupB);
            
            if (params->m->getControl_pressed()) { break; }
            
            //save scores
            params->scores.push_back(randomData[0]);
        }
        
        delete randT;
    }
    catch(exception& e) {
        params->m->errorOut(e, "UnifracWeightedCommand", "driver");
        exit(1);
    }
}
/**************************************************************************************************/
vector<double> UnifracWeightedCommand::createProcesses(Tree* t, CountTable* ct, vector< vector<string> > namesOfGroupCombos, vector<vector<int> >& randomizedTreeNodes) {
	try {
        //breakdown work between processors
        vector<linePair> lines;
        int remainingPairs = namesOfGroupCombos.size();
        if (remainingPairs < processors) { processors = remainingPairs; }
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, numPairs)); //startIndex, numPairs
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<weightedRandomData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            CountTable* copyCount = new CountTable();
            copyCount->copy(ct);
            vector<string> Treenames = t->getTreeNames();
            Tree* copyTree = new Tree(copyCount, Treenames);
            copyTree->getCopy(t);
            
            weightedRandomData* dataBundle = new weightedRandomData(lines[i+1].start, lines[i+1].end, namesOfGroupCombos, copyTree, copyCount, includeRoot, Groups, randomizedTreeNodes);
            data.push_back(dataBundle);
            workerThreads.push_back(new std::thread(driverWeightedRandom, dataBundle));
        }
        
        weightedRandomData* dataBundle = new weightedRandomData(lines[0].start, lines[0].end, namesOfGroupCombos, t, ct, includeRoot, Groups, randomizedTreeNodes);
        driverWeightedRandom(dataBundle);
        vector<double> scores = dataBundle->scores;
        
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            scores.insert(scores.end(), data[i]->scores.begin(), data[i]->scores.end());
            
            delete data[i]->t; delete data[i]->ct; delete data[i]; delete workerThreads[i];
        }
        delete dataBundle;
        return scores;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "createProcesses");
		exit(1);
	}
}
/***********************************************************/
void UnifracWeightedCommand::printWSummaryFile(int treeIndex, vector<double> utreeScores, vector<double> WScoreSig, vector<string> groupComb) {
	try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getSimpleName(treefile);
        variables["[tag]"] = toString(treeIndex);
        sumFile = getOutputFileName("wsummary",variables);
        outputNames.push_back(sumFile);  outputTypes["wsummary"].push_back(sumFile);
        ofstream outSum; util.openOutputFile(sumFile, outSum);
        
		//column headers
		outSum << "Tree#" << '\t' << "Groups" << '\t' << "WScore" << '\t';
		m->mothurOut("Tree#\tGroups\tWScore\t");
		if (random) { outSum << "WSig"; m->mothurOut("WSig"); }
		outSum << endl; m->mothurOutEndLine();
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
        int precisionLength = itersString.length();
        for (int j = 0; j < numComp; j++) {
            if (random) {
                if (WScoreSig[j] > (1/(float)iters)) {
                    outSum << setprecision(6) << treeIndex << '\t' << groupComb[j] << '\t' << utreeScores[j] << '\t' << setprecision(precisionLength) << WScoreSig[j] << endl;
                    cout << setprecision(6) << treeIndex << '\t' << groupComb[j] << '\t' << utreeScores[j] << '\t' << setprecision(precisionLength) << WScoreSig[j] << endl;
                    m->mothurOutJustToLog(toString(treeIndex) +"\t" + groupComb[j] +"\t" + toString(utreeScores[j]) +"\t" +  toString(WScoreSig[j]) + "\n");
                }else{
                    outSum << setprecision(6) << treeIndex << '\t' << groupComb[j] << '\t' << utreeScores[j] << '\t' << setprecision(precisionLength) << "<" << (1/float(iters)) << endl;
                    cout << setprecision(6) << treeIndex << '\t' << groupComb[j] << '\t' << utreeScores[j] << '\t' << setprecision(precisionLength) << "<" << (1/float(iters)) << endl;
                    m->mothurOutJustToLog(toString(treeIndex) +"\t" + groupComb[j] +"\t" + toString(utreeScores[j]) +"\t<" +  toString((1/float(iters))) + "\n");
                }
            }else{
                outSum << setprecision(6) << treeIndex << '\t' << groupComb[j] << '\t' << utreeScores[j] << endl;
                cout << setprecision(6) << treeIndex << '\t' << groupComb[j] << '\t' << utreeScores[j]  << endl;
                m->mothurOutJustToLog(toString(treeIndex) +"\t" + groupComb[j] +"\t" + toString(utreeScores[j]) +"\n");
            }
        }
		outSum.close();
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "printWSummaryFile");
		exit(1);
	}
}
/***********************************************************/
void UnifracWeightedCommand::createPhylipFile(int treeIndex, vector<double> utreeScores) {
	try {
		int count = 0;
        int numGroups = Groups.size();
        
        string phylipFileName;
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getSimpleName(treefile);
        variables["[tag]"] = toString(treeIndex);
        if ((outputForm == "lt") || (outputForm == "square")) {
            variables["[tag2]"] = "weighted.phylip";
            phylipFileName = getOutputFileName("phylip",variables);
            outputNames.push_back(phylipFileName); outputTypes["phylip"].push_back(phylipFileName);
        }else { //column
            variables["[tag2]"] = "weighted.column";
            phylipFileName = getOutputFileName("column",variables);
            outputNames.push_back(phylipFileName); outputTypes["column"].push_back(phylipFileName);
        }
        
        ofstream out; util.openOutputFile(phylipFileName, out);
        
        if ((outputForm == "lt") || (outputForm == "square")) { out << numGroups << endl; }
        
        //make matrix with scores in it
        vector< vector<float> > dists;	dists.resize(numGroups);
        for (int i = 0; i < numGroups; i++) { dists[i].resize(numGroups, 0.0); }
        
        //flip it so you can print it
        for (int r=0; r< numGroups; r++) { for (int l = 0; l < r; l++) { dists[r][l] = utreeScores[count]; dists[l][r] = utreeScores[count]; count++; } }
        
        //output to file
        for (int r=0; r<numGroups; r++) {
            string name = Groups[r];
            if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } } //pad with spaces to make compatible
            
            if (outputForm == "lt")          { out << name; for (int l = 0; l < r; l++)         {	out << '\t' << dists[r][l];     } out << endl;  }
            else if (outputForm == "square") { out << name; for (int l = 0; l < numGroups; l++) {	out << '\t' << dists[r][l];     } out << endl;  }
            else{
                for (int l = 0; l < r; l++) {
                    string otherName = Groups[l];
                    if (otherName.length() < 10) { while (otherName.length() < 10) {  otherName += " ";  } }  //pad with spaces to make compatible
                    
                    out  << name << '\t' << otherName  << '\t' << dists[r][l] << endl;  
                }
            }
        }
        out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "createPhylipFile");
		exit(1);
	}
}
/***********************************************************/
int UnifracWeightedCommand::findIndex(float score, int index, vector< vector<double> >& rScores) {
	try{
        int results = rScores[index].size();
        
		for (int i = 0; i < rScores[index].size(); i++) { if (rScores[index][i] >= score)	{	results = i; 	break; } }
        
		return results;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "findIndex");
		exit(1);
	}
}
/***********************************************************/
void UnifracWeightedCommand::calculateFreqsCumuls(set<double>& validScores, vector< vector<double> > rScores, vector< map<double, double> >& rScoreFreq, vector< map<double, double> >& rCumul) {
	try {
		//clear out old tree values
		rScoreFreq.clear(); rScoreFreq.resize(numComp); rCumul.clear(); rCumul.resize(numComp); validScores.clear();
	
		//calculate frequency
		for (int f = 0; f < numComp; f++) {
			for (int i = 0; i < rScores[f].size(); i++) { //looks like 0,0,1,1,1,2,4,7...  you want to make a map that say rScoreFreq[0] = 2, rScoreFreq[1] = 3...
				validScores.insert(rScores[f][i]);
				map<double,double>::iterator it = rScoreFreq[f].find(rScores[f][i]);
                
				if (it != rScoreFreq[f].end())  { rScoreFreq[f][rScores[f][i]]++;   }
				else                            { rScoreFreq[f][rScores[f][i]] = 1; }
			}
		}
		
		//calculate rcumul
		for(int a = 0; a < numComp; a++) {
			float rcumul = 1.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (set<double>::iterator it = validScores.begin(); it != validScores.end(); it++) {
				//make rscoreFreq map and rCumul
				map<double,double>::iterator it2 = rScoreFreq[a].find(*it);
				rCumul[a][*it] = rcumul;
				//get percentage of random trees with that info
				if (it2 != rScoreFreq[a].end()) {  rScoreFreq[a][*it] /= iters; rcumul-= it2->second;  }
				else { rScoreFreq[a][*it] = 0.0000; } //no random trees with that score
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "calculateFreqsCumuls");
		exit(1);
	}
}
/***********************************************************/
