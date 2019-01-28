/*
 *  unifracunweightedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unifracunweightedcommand.h"
#include "treereader.h"
#include "subsample.h"
#include "consensus.h"

//**********************************************************************************************************************
vector<string> UnifracUnweightedCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none","unweighted-uwsummary",false,true,true); parameters.push_back(ptree);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter prandom("random", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prandom);
		CommandParameter pdistance("distance", "Multiple", "column-lt-square-phylip", "column", "", "", "","phylip-column",false,false); parameters.push_back(pdistance);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","",false,false); parameters.push_back(psubsample);
        CommandParameter pconsensus("consensus", "Boolean", "", "F", "", "", "","tree",false,false); parameters.push_back(pconsensus);
        CommandParameter proot("root", "Boolean", "F", "", "", "", "","",false,false); parameters.push_back(proot);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string UnifracUnweightedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The unifrac.unweighted command parameters are tree, group, name, count, groups, iters, distance, processors, root and random.  tree parameter is required unless you have valid current tree file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n";
		helpString += "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n";
		helpString += "The distance parameter allows you to create a distance file from the results. The default is false. You may set distance to lt, square or column.\n";
		helpString += "The random parameter allows you to shut off the comparison to random trees. The default is false, meaning compare don't your trees with randomly generated trees.\n";
		helpString += "The root parameter allows you to include the entire root in your calculations. The default is false, meaning stop at the root for this comparision instead of the root of the entire tree.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The unifrac.unweighted command should be in the following format: unifrac.unweighted(groups=yourGroups, iters=yourIters).\n";
        helpString += "The subsample parameter allows you to enter the size pergroup of the sample or you can set subsample=T and mothur will use the size of your smallest group. The subsample parameter may only be used with a group file.\n";
        helpString += "The consensus parameter allows you to indicate you would like trees built from distance matrices created with the results of the subsampling, as well as a consensus tree built from these trees. Default=F.\n";
		helpString += "Example unifrac.unweighted(groups=A-B-C, iters=500).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and iters is 1000.\n";
		helpString += "The unifrac.unweighted command output two files: .unweighted and .uwsummary their descriptions are in the manual.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string UnifracUnweightedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        if (type == "unweighted")            {  pattern = "[filename],unweighted-[filename],[tag],unweighted";   }
        else if (type == "uwsummary")        {  pattern = "[filename],uwsummary";   }
        else if (type == "phylip")           {  pattern = "[filename],[tag],[tag2],dist";   }
        else if (type == "column")           {  pattern = "[filename],[tag],[tag2],dist";   }
        else if (type == "tree")             {  pattern = "[filename],[tag],[tag2],tre";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "UnifracUnweightedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
UnifracUnweightedCommand::UnifracUnweightedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["unweighted"] = tempOutNames;
		outputTypes["uwsummary"] = tempOutNames;
		outputTypes["phylip"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
        outputTypes["tree"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "UnifracUnweightedCommand");
		exit(1);
	}
}
/***********************************************************/
UnifracUnweightedCommand::UnifracUnweightedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["unweighted"] = tempOutNames;
			outputTypes["uwsummary"] = tempOutNames;
			outputTypes["phylip"] = tempOutNames;
			outputTypes["column"] = tempOutNames;
            outputTypes["tree"] = tempOutNames;
			
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
            //check for required parameters
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { abort = true; }
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
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
			
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(treefile);	}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
				
			itersString = validParameter.valid(parameters, "iters");				if (itersString == "not found") { itersString = "1000"; }
			util.mothurConvert(itersString, iters); 
			
			string temp = validParameter.valid(parameters, "distance");
			if (temp == "not found") { phylip = false; outputForm = ""; }
			else{
                if (temp=="phylip") { temp = "lt"; }
				if ((temp == "lt") || (temp == "column") || (temp == "square")) {  phylip = true;  outputForm = temp; }
				else { m->mothurOut("Options for distance are: lt, square, or column. Using lt."); m->mothurOutEndLine(); phylip = true; outputForm = "lt"; }
			}
			
			temp = validParameter.valid(parameters, "random");					if (temp == "not found") { temp = "f"; }
			random = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "root");					if (temp == "not found") { temp = "F"; }
			includeRoot = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			current->setProcessors(temp);
			util.mothurConvert(temp, processors); 
			
            temp = validParameter.valid(parameters, "subsample");		if (temp == "not found") { temp = "F"; }
			if (util.isNumeric1(temp)) { util.mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (util.isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
			
            if (!subsample) { subsampleIters = 0;   }
            else { subsampleIters = iters;          }
            
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

			if (!random) {  iters = 0;  } //turn off random calcs
			
			if (countfile=="") {
                if (namefile == "") {
                    vector<string> files; files.push_back(treefile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                } 
            }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "UnifracUnweightedCommand");
		exit(1);
	}
}

/***********************************************************/
int UnifracUnweightedCommand::execute() {
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

        TreeReader* reader = NULL;
        if (countfile == "") { reader = new TreeReader(treefile, groupfile, namefile); }
        else { reader = new TreeReader(treefile, countfile); }
        
        vector<Tree*> T; T = reader->getTrees();   //user trees
        CountTable* ct; ct = T[0]->getCountTable();
        if ((Groups.size() == 0) || (Groups.size() < 2)) {  Groups = ct->getNamesOfGroups();  } //must have at least 2 groups to compare
        delete reader;
        
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(treefile));
		sumFile = getOutputFileName("uwsummary",variables);
		outputNames.push_back(sumFile); outputTypes["uwsummary"].push_back(sumFile);
		util.openOutputFile(sumFile, outSum);

		long start = time(NULL);
        
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
        if ((Groups.size() == 0) || (Groups.size() < 2)) {  Groups = ct->getNamesOfGroups();  } //must have at least 2 groups to compare
		Unweighted unweighted(includeRoot, Groups);
        util.getCombos(groupComb, Groups, numComp);
        numGroups = Groups.size();
        
		if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
        
		if (numComp < processors) { processors = numComp;  m->mothurOut("Reducing processors to " + toString(numComp) + ".\n"); }
        
        if (consensus && (numComp < 2)) { m->mothurOut("consensus can only be used with numComparisions greater than 1, setting consensus=f.\n"); consensus=false; }
		
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" <<'\t';
		m->mothurOut("Tree#\tGroups\tUWScore\t");
		if (random) { outSum << "UWSig"; m->mothurOut("UWSig"); }
		outSum << endl; m->mothurOutEndLine();
	 
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			if (m->getControl_pressed()) { break; }
			
            counter = 0; rscoreFreq.resize(numComp);   rCumul.resize(numComp);   utreeScores.resize(numComp);   UWScoreSig.resize(numComp);
            
            vector<double> userData; userData.resize(numComp,0);  //weighted score info for user tree. data[0] = weightedscore AB, data[1] = weightedscore AC...

			userData = unweighted.getValues(T[i], processors, outputDir);  //userData[0] = unweightedscore
		
			if (m->getControl_pressed()) { break; }
			
			//output scores for each combination
			for(int k = 0; k < numComp; k++) {
				//saves users score
				utreeScores[k].push_back(userData[k]);
				
				//add users score to validscores
				validScores[userData[k]] = userData[k];
                
                if (!random) { UWScoreSig[k].push_back(0.0);	}
			}
            
            if (random) {  runRandomCalcs(T[i], userData);  }
			
			if (m->getControl_pressed()) { break; }
            
            int startSubsample = time(NULL);
            
            //subsample loop
            SubSample sample;
            vector< vector<double> > calcDistsTotals;  //each iter, each groupCombos dists. this will be used to make .dist files
            for (int thisIter = 0; thisIter < subsampleIters; thisIter++) { //subsampleIters=0, if subsample=f.
                if (m->getControl_pressed()) { break; }
                
                //copy to preserve old one - would do this in subsample but memory cleanup becomes messy.
                CountTable* newCt = new CountTable();
                 
                //uses method of setting groups to doNotIncludeMe
                int sampleTime = 0;
                if (m->getDebug()) { sampleTime = time(NULL); }
                
                Tree* subSampleTree = sample.getSample(T[i], ct, newCt, subsampleSize, Groups);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: iter " + toString(thisIter) + " took " + toString(time(NULL) - sampleTime) + " seconds to sample tree.\n"); }
                
                //call new weighted function
                vector<double> iterData; iterData.resize(numComp,0);
                Unweighted thisUnweighted(includeRoot, Groups);
                iterData = thisUnweighted.getValues(subSampleTree, processors, outputDir); //userData[0] = weightedscore
        
                //save data to make ave dist, std dist
                calcDistsTotals.push_back(iterData);
                
                delete newCt;
                delete subSampleTree;
                
                if((thisIter+1) % 100 == 0){	m->mothurOutJustToScreen(toString(thisIter+1)+"\n"); 		}
            }
            if (subsample) { m->mothurOut("It took " + toString(time(NULL) - startSubsample) + " secs to run the subsampling."); m->mothurOutEndLine(); }
            
            if (m->getControl_pressed()) { break; }

            if (subsample) {  getAverageSTDMatrices(calcDistsTotals, i); }
            if (consensus) {  getConsensusTrees(calcDistsTotals, i);  }
            
            //print output files
			printUWSummaryFile(i);
			if (random)  {	printUnweightedFile(i+1);	}
			if (phylip) {	createPhylipFile(i);		}
			
			rscoreFreq.clear();  rCumul.clear();  validScores.clear();  utreeScores.clear();  UWScoreSig.clear();
		}
		
		outSum.close();
		delete ct; 
		for (int i = 0; i < T.size(); i++) { delete T[i]; }
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  }	return 0; }
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run unifrac.unweighted."); m->mothurOutEndLine();
		
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
		m->errorOut(e, "UnifracUnweightedCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int UnifracUnweightedCommand::getAverageSTDMatrices(vector< vector<double> >& dists, int treeNum) {
	try {
        
        //we need to find the average distance and standard deviation for each groups distance
        //finds sum
        vector<double> averages = util.getAverages(dists);
        
        //find standard deviation
        vector<double> stdDev = util.getStandardDeviation(dists, averages);
        
        //make matrix with scores in it
        vector< vector<double> > avedists;	//avedists.resize(m->getNumGroups());
        for (int i = 0; i < numGroups; i++) {
            vector<double> temp;
            for (int j = 0; j < numGroups; j++) { temp.push_back(0.0); }
            avedists.push_back(temp);
        }
        
        //make matrix with scores in it
        vector< vector<double> > stddists;	//stddists.resize(m->getNumGroups());
        for (int i = 0; i < numGroups; i++) {
            vector<double> temp;
            for (int j = 0; j < numGroups; j++) { temp.push_back(0.0); }
            //stddists[i].resize(m->getNumGroups(), 0.0);
            stddists.push_back(temp);
        }
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: about to fill matrix.\n"); }
        
        //flip it so you can print it
        int count = 0;
        for (int r=0; r<numGroups; r++) {
            for (int l = 0; l < r; l++) {
                avedists[r][l] = averages[count];
                avedists[l][r] = averages[count];
                stddists[r][l] = stdDev[count];
                stddists[l][r] = stdDev[count];
                count++;
            }
        }
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: done filling matrix.\n"); }
        
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "unweighted.ave";
        string aveFileName = getOutputFileName("phylip",variables);
        if (outputForm != "column") { outputNames.push_back(aveFileName); outputTypes["phylip"].push_back(aveFileName);  }
        else { outputNames.push_back(aveFileName); outputTypes["column"].push_back(aveFileName);  }
        ofstream out;
        util.openOutputFile(aveFileName, out);
        
        variables["[tag2]"] = "unweighted.std";
        string stdFileName = getOutputFileName("phylip",variables);
        if (outputForm != "column") { outputNames.push_back(stdFileName); outputTypes["phylip"].push_back(stdFileName); }
        else { outputNames.push_back(stdFileName); outputTypes["column"].push_back(stdFileName); }
        ofstream outStd;
        util.openOutputFile(stdFileName, outStd);
        
        if ((outputForm == "lt") || (outputForm == "square")) {
            //output numSeqs
            out << numGroups << endl;
            outStd << numGroups << endl;
        }
        
        //output to file
        for (int r=0; r<numGroups; r++) {
            //output name
            string name = Groups[r];
            if (name.length() < 10) { //pad with spaces to make compatible
                while (name.length() < 10) {  name += " ";  }
            }
            
            if (outputForm == "lt") {
                out << name;
                outStd << name;
                
                //output distances
                for (int l = 0; l < r; l++) {	out  << '\t' << avedists[r][l];  outStd  << '\t' << stddists[r][l];}
                out << endl;  outStd << endl;
            }else if (outputForm == "square") {
                out << name;
                outStd << name;
                
                //output distances
                for (int l = 0; l < numGroups; l++) {	out  << '\t' << avedists[r][l]; outStd   << '\t' << stddists[r][l]; }
                out << endl; outStd << endl;
            }else{
                //output distances
                for (int l = 0; l < r; l++) {	
                    string otherName = Groups[l];
                    if (otherName.length() < 10) { //pad with spaces to make compatible
                        while (otherName.length() < 10) {  otherName += " ";  }
                    }
                    
                    out  << name << '\t' << otherName  << '\t' << avedists[r][l] << endl;
                    outStd  << name << '\t' << otherName  << '\t' << stddists[r][l] << endl;
                }
            }
        }
        out.close();
        outStd.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getAverageSTDMatrices");
		exit(1);
	}
}

/**************************************************************************************************/
int UnifracUnweightedCommand::getConsensusTrees(vector< vector<double> >& dists, int treeNum) {
	try {
        //create treemap class from groupmap for tree class to use
        CountTable newCt;
        set<string> nameMap;
        map<string, string> groupMap;
        set<string> gps;
        for (int i = 0; i < Groups.size(); i++) {
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
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "unweighted.cons";
        string conFile = getOutputFileName("tree",variables);				
        outputNames.push_back(conFile); outputTypes["tree"].push_back(conFile); 
        ofstream outTree;
        util.openOutputFile(conFile, outTree);
        
        if (conTree != NULL) { conTree->print(outTree, "boot"); delete conTree; }
        outTree.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getConsensusTrees");
		exit(1);
	}
}
/**************************************************************************************************/

vector<Tree*> UnifracUnweightedCommand::buildTrees(vector< vector<double> >& dists, int treeNum, CountTable& myct) {
	try {
        
        vector<Tree*> trees;
        
        //create a new filename
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "unweighted.all";
        string outputFile = getOutputFileName("tree",variables);				
        outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
        
        ofstream outAll;
        util.openOutputFile(outputFile, outAll);
        
        
        for (int i = 0; i < dists.size(); i++) { //dists[0] are the dists for the first subsampled tree.
            
            if (m->getControl_pressed()) { break; }
            
            //make matrix with scores in it
            vector< vector<double> > sims;	sims.resize(Groups.size());
            for (int j = 0; j < Groups.size(); j++) {
                sims[j].resize(Groups.size(), 0.0);
            }
            
            int count = 0;
			for (int r=0; r<Groups.size(); r++) {
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
            
            //print tree
            tempTree->print(outAll);
        }
        
        outAll.close();
        
        if (m->getControl_pressed()) {  for (int i = 0; i < trees.size(); i++) {  delete trees[i]; trees[i] = NULL; } util.mothurRemove(outputFile); }
        
        return trees;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "buildTrees");
		exit(1);
	}
}
/**************************************************************************************************/

int UnifracUnweightedCommand::runRandomCalcs(Tree* thisTree, vector<double> usersScores) {
	try {
        vector<double> randomData; randomData.resize(numComp,0); //weighted score info for random trees. data[0] = weightedscore AB, data[1] = weightedscore AC...
        
        Unweighted unweighted(includeRoot, Groups);
        
        vector< vector<string> > namesOfGroupCombos;
        numGroups = Groups.size();
        for (int a=0; a<numGroups; a++) {
            for (int l = 0; l < a; l++) {
                vector<string> groups; groups.push_back(Groups[a]); groups.push_back(Groups[l]);
                namesOfGroupCombos.push_back(groups);
            }
        }
        vector<vector<int> > randomTreeNodes;
        for (int f = 0; f < numComp; f++) {
            vector<int> randomNodesForThisCombo = thisTree->getNodes(namesOfGroupCombos[f]);
            randomTreeNodes.push_back(randomNodesForThisCombo);
        }
        vector<vector<int> > savedRandomTreeNodes = randomTreeNodes;
        
        //get unweighted scores for random trees - if random is false iters = 0
        for (int j = 0; j < iters; j++) {
            
            randomTreeNodes = savedRandomTreeNodes;
            
            for (int f = 0; f < numComp; f++) { util.mothurRandomShuffle(randomTreeNodes[f]);  } //randomize labels
            
            //we need a different getValues because when we swap the labels we only want to swap those in each pairwise comparison
            randomData = unweighted.getValues(thisTree, randomTreeNodes, processors, outputDir);
            
            if (m->getControl_pressed()) { return 0; }
			
            for(int k = 0; k < numComp; k++) {	
                //add trees unweighted score to map of scores
                map<float,float>::iterator it = rscoreFreq[k].find(randomData[k]);
                if (it != rscoreFreq[k].end()) {//already have that score
                    rscoreFreq[k][randomData[k]]++;
                }else{//first time we have seen this score
                    rscoreFreq[k][randomData[k]] = 1;
                }
				
                //add randoms score to validscores
                validScores[randomData[k]] = randomData[k];
            }
        }
        
        for(int a = 0; a < numComp; a++) {
            float rcumul = 1.0000;
    
            //this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
            for (map<float,float>::iterator it = validScores.begin(); it != validScores.end(); it++) {
                //make rscoreFreq map and rCumul
                map<float,float>::iterator it2 = rscoreFreq[a].find(it->first);
                rCumul[a][it->first] = rcumul;
                //get percentage of random trees with that info
                if (it2 != rscoreFreq[a].end()) {  rscoreFreq[a][it->first] /= iters; rcumul-= it2->second;  }
                else { rscoreFreq[a][it->first] = 0.0000; } //no random trees with that score
            }
            UWScoreSig[a].push_back(rCumul[a][usersScores[a]]);
        }
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "runRandomCalcs");
		exit(1);
	}
}
/***********************************************************/
void UnifracUnweightedCommand::printUnweightedFile(int treeNum) {
	try {
		vector<double> data;
		vector<string> tags;
		
		tags.push_back("Score");
		tags.push_back("RandFreq"); tags.push_back("RandCumul");
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getSimpleName(treefile);
        variables["[tag]"] = toString(treeNum);
        string unFileName = getOutputFileName("unweighted", variables);
        FileOutput* output = new ColumnFile(unFileName, itersString);
        outputNames.push_back(unFileName); outputTypes["unweighted"].push_back(unFileName);
			
		for(int a = 0; a < numComp; a++) {
			output->setLabelName(groupComb[a], tags);
			//print each line
			for (map<float,float>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]);						
				output->updateOutput(data);
				data.clear();
			} 
			output->resetFile();
		}
        delete output;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "printUnweightedFile");
		exit(1);
	}
}

/***********************************************************/
void UnifracUnweightedCommand::printUWSummaryFile(int i) {
	try {
				
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
			
		//print each line

		for(int a = 0; a < numComp; a++) {
			outSum << i+1 << '\t';
			m->mothurOut(toString(i+1) + "\t");
			
			if (random) {
				if (UWScoreSig[a][0] > (1/(float)iters)) {
					outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << UWScoreSig[a][0] << endl;
					cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << UWScoreSig[a][0] << endl; 
					m->mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t" + toString(UWScoreSig[a][0])+ "\n"); 
				}else {
					outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl;
					cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
					m->mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t<" + toString((1/float(iters))) + "\n"); 
				}
			}else{
				outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0]  << endl;
				cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0]  << endl; 
				m->mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0]) + "\n");
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "printUWSummaryFile");
		exit(1);
	}
}
/***********************************************************/
void UnifracUnweightedCommand::createPhylipFile(int i) {
	try {
		string phylipFileName;
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getSimpleName(treefile);
        variables["[tag]"] = toString(i+1);
		if ((outputForm == "lt") || (outputForm == "square")) {
            variables["[tag2]"] = "unweighted.phylip";
			phylipFileName = getOutputFileName("phylip",variables);
			outputNames.push_back(phylipFileName); outputTypes["phylip"].push_back(phylipFileName); 
		}else { //column
            variables["[tag2]"] = "unweighted.column";
			phylipFileName = getOutputFileName("column",variables);
			outputNames.push_back(phylipFileName); outputTypes["column"].push_back(phylipFileName); 
		}
		
		ofstream out;
		util.openOutputFile(phylipFileName, out);
		
        numGroups = Groups.size();
		if ((outputForm == "lt") || (outputForm == "square")) {
			//output numSeqs
			out << numGroups << endl;
		}
		
		//make matrix with scores in it
		vector< vector<float> > dists;	dists.resize(numGroups);
		for (int i = 0; i < numGroups; i++) {
			dists[i].resize(numGroups, 0.0);
		}
		
		//flip it so you can print it
		int count = 0;
		for (int r=0; r<numGroups; r++) {
			for (int l = 0; l < r; l++) {
				dists[r][l] = utreeScores[count][0];
				dists[l][r] = utreeScores[count][0];
				count++;
			}
		}
		
		//output to file
		for (int r=0; r<numGroups; r++) {
			//output name
			string name = Groups[r];
            //pad with spaces to make compatible
			if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } }
			
			if (outputForm == "lt") { //output distances
				out << name;
				for (int l = 0; l < r; l++) {	out  << '\t' << dists[r][l];  }   out << endl;
			}else if (outputForm == "square") {  //output distances
				out << name;
				for (int l = 0; l < numGroups; l++) {	out  << '\t' << dists[r][l];  }  out << endl;
			}else{
				//output distances
				for (int l = 0; l < r; l++) {	
					string otherName = Groups[l];
                    //pad with spaces to make compatible
					if (otherName.length() < 10) { while (otherName.length() < 10) {  otherName += " ";  } }
					
					out  << name << '\t' << otherName << '\t' << dists[r][l] << endl;
				}
			}
		}
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "createPhylipFile");
		exit(1);
	}
}
/***********************************************************/




