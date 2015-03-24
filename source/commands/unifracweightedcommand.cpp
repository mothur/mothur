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
        CommandParameter pconsensus("consensus", "Boolean", "", "F", "", "", "","tree",false,false); parameters.push_back(pconsensus);
        CommandParameter prandom("random", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prandom);
		CommandParameter pdistance("distance", "Multiple", "column-lt-square-phylip", "column", "", "", "","phylip-column",false,false); parameters.push_back(pdistance);
		CommandParameter proot("root", "Boolean", "F", "", "", "", "","",false,false); parameters.push_back(proot);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
        helpString += "The consensus parameter allows you to indicate you would like trees built from distance matrices created with the results, as well as a consensus tree built from these trees. Default=F.\n";
        helpString += "The unifrac.weighted command should be in the following format: unifrac.weighted(groups=yourGroups, iters=yourIters).\n";
		helpString += "Example unifrac.weighted(groups=A-B-C, iters=500).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and iters is 1000.\n";
		helpString += "The unifrac.weighted command output two files: .weighted and .wsummary their descriptions are in the manual.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
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
        else if (type == "wsummary")        {  pattern = "[filename],wsummary";   }
        else if (type == "phylip")           {  pattern = "[filename],[tag],[tag2],dist";   }
        else if (type == "column")           {  pattern = "[filename],[tag],[tag2],dist";   }
        else if (type == "tree")             {  pattern = "[filename],[tag],[tag2],tre";   }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "UnifracWeightedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
UnifracWeightedCommand::UnifracWeightedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["weighted"] = tempOutNames;
		outputTypes["wsummary"] = tempOutNames;
		outputTypes["phylip"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
        outputTypes["tree"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "UnifracWeightedCommand");
		exit(1);
	}
}

/***********************************************************/
UnifracWeightedCommand::UnifracWeightedCommand(string option) {
	try {
		abort = false; calledHelp = false;   
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters=parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["weighted"] = tempOutNames;
			outputTypes["wsummary"] = tempOutNames;
			outputTypes["phylip"] = tempOutNames;
			outputTypes["column"] = tempOutNames;
            outputTypes["tree"] = tempOutNames;
			
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
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
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { treefile = ""; abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = m->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { m->setTreeFile(treefile); }	
			
			//check for required parameters
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { m->setCountTableFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(treefile);	}
			
																	
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
				
			itersString = validParameter.validFile(parameters, "iters", false);			if (itersString == "not found") { itersString = "1000"; }
			m->mothurConvert(itersString, iters); 
			
			string temp = validParameter.validFile(parameters, "distance", false);			
			if (temp == "not found") { phylip = false; outputForm = ""; }
			else{
                if (temp=="phylip") { temp = "lt"; }
				if ((temp == "lt") || (temp == "column") || (temp == "square")) {  phylip = true;  outputForm = temp; }
				else { m->mothurOut("Options for distance are: lt, square, or column. Using lt."); m->mothurOutEndLine(); phylip = true; outputForm = "lt"; }
			}
			
			temp = validParameter.validFile(parameters, "random", false);				if (temp == "not found") { temp = "F"; }
			random = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "root", false);					if (temp == "not found") { temp = "F"; }
			includeRoot = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
            temp = validParameter.validFile(parameters, "subsample", false);		if (temp == "not found") { temp = "F"; }
			if (m->isNumeric1(temp)) { m->mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (m->isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
			
            if (!subsample) { subsampleIters = 0;   }
            else { subsampleIters = iters;          }
            
            temp = validParameter.validFile(parameters, "consensus", false);					if (temp == "not found") { temp = "F"; }
			consensus = m->isTrue(temp);
            
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
            
			if (countfile=="") {
                if (namefile == "") {
                    vector<string> files; files.push_back(treefile);
                    parser.getNameFile(files);
                } 
            }
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
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		m->setTreeFile(treefile);
		
        TreeReader* reader;
        if (countfile == "") { reader = new TreeReader(treefile, groupfile, namefile); }
        else { reader = new TreeReader(treefile, countfile); }
        T = reader->getTrees();
        ct = T[0]->getCountTable();
        delete reader;
        
        if (m->control_pressed) {  delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } return 0; }
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getSimpleName(treefile);
		sumFile = getOutputFileName("wsummary",variables);
		m->openOutputFile(sumFile, outSum);
		outputNames.push_back(sumFile);  outputTypes["wsummary"].push_back(sumFile);
		
        SharedUtil util;
		string s; //to make work with setgroups
		Groups = m->getGroups();
		vector<string> nameGroups = ct->getNamesOfGroups();
        if (nameGroups.size() < 2) { m->mothurOut("[ERROR]: You cannot run unifrac.weighted with less than 2 groups, aborting.\n"); delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } return 0; }
 		util.setGroups(Groups, nameGroups, s, numGroups, "weighted");	//sets the groups the user wants to analyze
		m->setGroups(Groups);
		
        if (m->control_pressed) {  delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } return 0; }
        
		Weighted weighted(includeRoot);
			
		int start = time(NULL);
            
        //set or check size
        if (subsample) {
            //user has not set size, set size = smallest samples size
            if (subsampleSize == -1) { 
                vector<string> temp; temp.push_back(Groups[0]);
                subsampleSize = ct->getGroupCount(Groups[0]); //num in first group
                for (int i = 1; i < Groups.size(); i++) {
                    int thisSize = ct->getGroupCount(Groups[i]);
                    if (thisSize < subsampleSize) {	subsampleSize = thisSize;	}
                }
                m->mothurOut("\nSetting subsample size to " + toString(subsampleSize) + ".\n\n");
            }else { //eliminate any too small groups
                vector<string> newGroups = Groups;
                Groups.clear();
                for (int i = 0; i < newGroups.size(); i++) {
                    int thisSize = ct->getGroupCount(newGroups[i]);
                    
                    if (thisSize >= subsampleSize) {    Groups.push_back(newGroups[i]);	}
                    else {   m->mothurOut("You have selected a size that is larger than "+newGroups[i]+" number of sequences, removing "+newGroups[i]+".\n"); }
                } 
                m->setGroups(Groups);
            }
        }
        
        //here in case some groups are removed by subsample
        util.getCombos(groupComb, Groups, numComp);
        
        if (numComp < processors) { processors = numComp; }
        
        if (consensus && (numComp < 2)) { m->mothurOut("consensus can only be used with numComparisions greater than 1, setting consensus=f.\n"); consensus=false; }
        
        //get weighted scores for users trees
        for (int i = 0; i < T.size(); i++) {
            
            if (m->control_pressed) { delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
            
            counter = 0;
            rScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
            uScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
            
            vector<double> userData; userData.resize(numComp,0);  //weighted score info for user tree. data[0] = weightedscore AB, data[1] = weightedscore AC...
            vector<double> randomData; randomData.resize(numComp,0); //weighted score info for random trees. data[0] = weightedscore AB, data[1] = weightedscore AC...
            
            if (random) {  
                variables["[filename]"] = outputDir + m->getSimpleName(treefile);
                variables["[tag]"] = toString(i+1);
                string wFileName = getOutputFileName("weighted", variables);
                output = new ColumnFile(wFileName, itersString);
				outputNames.push_back(wFileName); outputTypes["weighted"].push_back(wFileName);
            } 
            
            userData = weighted.getValues(T[i], processors, outputDir); //userData[0] = weightedscore
            if (m->control_pressed) { delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } if (random) { delete output; } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
            
            //save users score
            for (int s=0; s<numComp; s++) {
                //add users score to vector of user scores
                uScores[s].push_back(userData[s]);
                //save users tree score for summary file
                utreeScores.push_back(userData[s]);
            }
            
            if (random) {  runRandomCalcs(T[i], userData); }
            
            //clear data
            rScores.clear();
            uScores.clear();
            validScores.clear();
            
            //subsample loop
            vector< vector<double> > calcDistsTotals;  //each iter, each groupCombos dists. this will be used to make .dist files
            for (int thisIter = 0; thisIter < subsampleIters; thisIter++) { //subsampleIters=0, if subsample=f.
                if (m->control_pressed) { break; }
                
                //copy to preserve old one - would do this in subsample but memory cleanup becomes messy.
                CountTable* newCt = new CountTable();
                
                int sampleTime = 0;
                if (m->debug) { sampleTime = time(NULL); }
                
                //uses method of setting groups to doNotIncludeMe
                SubSample sample;
                Tree* subSampleTree = sample.getSample(T[i], ct, newCt, subsampleSize);
                
                if (m->debug) { m->mothurOut("[DEBUG]: iter " + toString(thisIter) + " took " + toString(time(NULL) - sampleTime) + " seconds to sample tree.\n"); }
                
                //call new weighted function
                vector<double> iterData; iterData.resize(numComp,0);
                Weighted thisWeighted(includeRoot);
                iterData = thisWeighted.getValues(subSampleTree, processors, outputDir); //userData[0] = weightedscore
                
                //save data to make ave dist, std dist
                calcDistsTotals.push_back(iterData);
                
                delete newCt;
                delete subSampleTree;
                
                if((thisIter+1) % 100 == 0){	m->mothurOutJustToScreen(toString(thisIter+1)+"\n"); 	}
            }
            
            if (m->control_pressed) { delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } if (random) { delete output; } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
            
            if (subsample) {  getAverageSTDMatrices(calcDistsTotals, i); }
            if (consensus) {  getConsensusTrees(calcDistsTotals, i);  }
        }
        
		
		if (m->control_pressed) { delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0;  }
		
        if (phylip) {	createPhylipFile();		}
    
		printWSummaryFile();
		
		//clear out users groups
		m->clearGroups();
		delete ct; 
		for (int i = 0; i < T.size(); i++) { delete T[i]; }
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run unifrac.weighted."); m->mothurOutEndLine();
		
		//set phylip file as new current phylipfile
		string current = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setPhylipFile(current); }
		}
		
		//set column file as new current columnfile
		itTypes = outputTypes.find("column");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setColumnFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
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
        //we need to find the average distance and standard deviation for each groups distance
        
        //finds sum
        vector<double> averages = m->getAverages(dists);        
        
        //find standard deviation
        vector<double> stdDev = m->getStandardDeviation(dists, averages);
        
        //make matrix with scores in it
        vector< vector<double> > avedists;	//avedists.resize(m->getNumGroups());
        for (int i = 0; i < m->getNumGroups(); i++) {
            vector<double> temp;
            for (int j = 0; j < m->getNumGroups(); j++) { temp.push_back(0.0); }
            avedists.push_back(temp);
        }
        
        //make matrix with scores in it
        vector< vector<double> > stddists;	//stddists.resize(m->getNumGroups());
        for (int i = 0; i < m->getNumGroups(); i++) {
            vector<double> temp;
            for (int j = 0; j < m->getNumGroups(); j++) { temp.push_back(0.0); }
            //stddists[i].resize(m->getNumGroups(), 0.0);
            stddists.push_back(temp);
        }

        
        //flip it so you can print it
        int count = 0;
        for (int r=0; r<m->getNumGroups(); r++) { 
            for (int l = 0; l < r; l++) {
                avedists[r][l] = averages[count];
                avedists[l][r] = averages[count];
                stddists[r][l] = stdDev[count];
                stddists[l][r] = stdDev[count];
                count++;
            }
        }
        
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getSimpleName(treefile);
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "weighted.ave";
        string aveFileName = getOutputFileName("phylip",variables);
        if (outputForm != "column") { outputNames.push_back(aveFileName); outputTypes["phylip"].push_back(aveFileName);  }
        else { outputNames.push_back(aveFileName); outputTypes["column"].push_back(aveFileName);  }
        ofstream out;
        m->openOutputFile(aveFileName, out);
        
        variables["[tag2]"] = "weighted.std";
        string stdFileName = getOutputFileName("phylip",variables);
        if (outputForm != "column") { outputNames.push_back(stdFileName); outputTypes["phylip"].push_back(stdFileName); }
        else { outputNames.push_back(stdFileName); outputTypes["column"].push_back(stdFileName); }        
        ofstream outStd;
        m->openOutputFile(stdFileName, outStd);
        
        if ((outputForm == "lt") || (outputForm == "square")) {
            //output numSeqs
            out << m->getNumGroups() << endl;
            outStd << m->getNumGroups() << endl;
        }
        
        //output to file
        for (int r=0; r<m->getNumGroups(); r++) { 
            //output name
            string name = (m->getGroups())[r];
            if (name.length() < 10) { //pad with spaces to make compatible
                while (name.length() < 10) {  name += " ";  }
            }
            
            if (outputForm == "lt") {
                out << name << '\t';
                outStd << name << '\t';
                
                //output distances
                for (int l = 0; l < r; l++) {	out  << avedists[r][l] << '\t';  outStd  << stddists[r][l] << '\t';}
                out << endl;  outStd << endl;
            }else if (outputForm == "square") {
                out << name << '\t';
                outStd << name << '\t';
                
                //output distances
                for (int l = 0; l < m->getNumGroups(); l++) {	out  << avedists[r][l] << '\t'; outStd  << stddists[r][l] << '\t'; }
                out << endl; outStd << endl;
            }else{
                //output distances
                for (int l = 0; l < r; l++) {	
                    string otherName = (m->getGroups())[l];
                    if (otherName.length() < 10) { //pad with spaces to make compatible
                        while (otherName.length() < 10) {  otherName += " ";  }
                    }
                    
                    out  << name << '\t' << otherName << avedists[r][l] << endl;  
                    outStd  << name << '\t' << otherName << stddists[r][l] << endl; 
                }
            }
        }
        out.close();
        outStd.close();
        
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
        
        //used in tree constructor 
        m->runParse = false;
        
        ///create treemap class from groupmap for tree class to use
        CountTable newCt;
        set<string> nameMap;
        map<string, string> groupMap;
        set<string> gps;
        for (int i = 0; i < m->getGroups().size(); i++) { 
            nameMap.insert(m->getGroups()[i]); 
            gps.insert(m->getGroups()[i]); 
            groupMap[m->getGroups()[i]] = m->getGroups()[i];
        }
        newCt.createTable(nameMap, groupMap, gps);
        
        //clear  old tree names if any
        m->Treenames.clear();
        
        //fills globaldatas tree names
        m->Treenames = m->getGroups();
        
        vector<Tree*> newTrees = buildTrees(dists, treeNum, newCt); //also creates .all.tre file containing the trees created
        
        if (m->control_pressed) { return 0; }
        
        Consensus con;
        Tree* conTree = con.getTree(newTrees);
        
        //create a new filename
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "weighted.cons";
        string conFile = getOutputFileName("tree",variables);							
        outputNames.push_back(conFile); outputTypes["tree"].push_back(conFile); 
        ofstream outTree;
        m->openOutputFile(conFile, outTree);
        
        if (conTree != NULL) { conTree->print(outTree, "boot"); delete conTree; }
        outTree.close();
        
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
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(treefile));
        variables["[tag]"] = toString(treeNum+1);
        variables["[tag2]"] = "weighted.all";
        string outputFile = getOutputFileName("tree",variables);				
        outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile); 
        
        ofstream outAll;
        m->openOutputFile(outputFile, outAll);
        

        for (int i = 0; i < dists.size(); i++) { //dists[0] are the dists for the first subsampled tree.
            
            if (m->control_pressed) { break; }
            
            //make matrix with scores in it
            vector< vector<double> > sims;	sims.resize(m->getNumGroups());
            for (int j = 0; j < m->getNumGroups(); j++) {
                sims[j].resize(m->getNumGroups(), 0.0);
            }
            
            int count = 0;
			for (int r=0; r<m->getNumGroups(); r++) { 
				for (int l = 0; l < r; l++) {
                    double sim = -(dists[i][count]-1.0);
					sims[r][l] = sim;
					sims[l][r] = sim;
					count++;
				}
			}

            //create tree
            Tree* tempTree = new Tree(&myct, sims);
            tempTree->assembleTree();
            
            trees.push_back(tempTree);
            
            //print tree
            tempTree->print(outAll);
        }
        
        outAll.close();
        
        if (m->control_pressed) {  for (int i = 0; i < trees.size(); i++) {  delete trees[i]; trees[i] = NULL; } m->mothurRemove(outputFile); }
        
        return trees;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "buildTrees");
		exit(1);
	}
}
/**************************************************************************************************/

int UnifracWeightedCommand::runRandomCalcs(Tree* thisTree, vector<double> usersScores) {
	try {
        
        //calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
        vector< vector<string> > namesOfGroupCombos;
        for (int a=0; a<numGroups; a++) { 
            for (int l = 0; l < a; l++) {	
                vector<string> groups; groups.push_back((m->getGroups())[a]); groups.push_back((m->getGroups())[l]);
                namesOfGroupCombos.push_back(groups);
            }
        }
        
        lines.clear();
        
        //breakdown work between processors
        int remainingPairs = namesOfGroupCombos.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, numPairs)); //startIndex, numPairs
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }

        
        
        //get scores for random trees
        for (int j = 0; j < iters; j++) {
//#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            //if(processors == 1){
              //  driver(thisTree,  namesOfGroupCombos, 0, namesOfGroupCombos.size(),  rScores);
           // }else{
                createProcesses(thisTree,  namesOfGroupCombos, rScores);
           // }
//#else
            //driver(thisTree, namesOfGroupCombos, 0, namesOfGroupCombos.size(), rScores);
//#endif
 
            
            if (m->control_pressed) { delete ct;  for (int i = 0; i < T.size(); i++) { delete T[i]; } delete output; outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
            
        }
        lines.clear();
        
        //find the signifigance of the score for summary file
        for (int f = 0; f < numComp; f++) {
            //sort random scores
            sort(rScores[f].begin(), rScores[f].end());
            
            //the index of the score higher than yours is returned 
            //so if you have 1000 random trees the index returned is 100 
            //then there are 900 trees with a score greater then you. 
            //giving you a signifigance of 0.900
            int index = findIndex(usersScores[f], f);    if (index == -1) { m->mothurOut("error in UnifracWeightedCommand"); m->mothurOutEndLine(); exit(1); } //error code
            
            //the signifigance is the number of trees with the users score or higher 
            WScoreSig.push_back((iters-index)/(float)iters);
        }
        
        //out << "Tree# " << i << endl;
        calculateFreqsCumuls();
        printWeightedFile();
        
        delete output;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "runRandomCalcs");
		exit(1);
	}
}
/**************************************************************************************************/

int UnifracWeightedCommand::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos, vector< vector<double> >& scores) {
	try {
        int process = 1;
		vector<int> processIDS;
		EstOutput results;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driver(t, namesOfGroupCombos, lines[process].start, lines[process].end, scores);
			
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputDir + m->mothurGetpid(process) + ".weightedcommand.results.temp";
				m->openOutputFile(tempFile, out);
				for (int i = lines[process].start; i < (lines[process].start + lines[process].end); i++) { out << scores[i][(scores[i].size()-1)] << '\t';  } out << endl;
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		driver(t, namesOfGroupCombos, lines[0].start, lines[0].end, scores);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<(processors-1);i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//get data created by processes
		for (int i=0;i<(processors-1);i++) { 
	
			ifstream in;
			string s = outputDir + toString(processIDS[i]) + ".weightedcommand.results.temp";
			m->openInputFile(s, in);
			
			double tempScore;
			for (int j = lines[(i+1)].start; j < (lines[(i+1)].start + lines[(i+1)].end); j++) { in >> tempScore; scores[j].push_back(tempScore); }
			in.close();
			m->mothurRemove(s);
		}
#else
        //fill in functions
        vector<weightedRandomData*> pDataArray;
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1];
        vector<CountTable*> cts;
        vector<Tree*> trees;
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
            CountTable* copyCount = new CountTable();
            copyCount->copy(ct);
            Tree* copyTree = new Tree(copyCount);
            copyTree->getCopy(t);
            
            cts.push_back(copyCount);
            trees.push_back(copyTree);
            
            vector< vector<double> > copyScores = rScores;
            
            weightedRandomData* tempweighted = new weightedRandomData(m, lines[i].start, lines[i].end, namesOfGroupCombos, copyTree, copyCount, includeRoot, copyScores);
			pDataArray.push_back(tempweighted);
			processIDS.push_back(i);
            
			hThreadArray[i-1] = CreateThread(NULL, 0, MyWeightedRandomThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
		}
		
		driver(t, namesOfGroupCombos, lines[0].start, lines[0].end, scores);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            for (int j = pDataArray[i]->start; j < (pDataArray[i]->start+pDataArray[i]->num); j++) {
                scores[j].push_back(pDataArray[i]->scores[j][pDataArray[i]->scores[j].size()-1]);
            }
			delete cts[i];
            delete trees[i];
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}

		
#endif	
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/
int UnifracWeightedCommand::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num, vector< vector<double> >& scores) { 
 try {
		Tree* randT = new Tree(ct);
     
        Weighted weighted(includeRoot);
     
		for (int h = start; h < (start+num); h++) {
	
			if (m->control_pressed) { return 0; }
		
			//initialize weighted score
			string groupA = namesOfGroupCombos[h][0]; 
			string groupB = namesOfGroupCombos[h][1];
			
			//copy T[i]'s info.
			randT->getCopy(t);
			 
			//create a random tree with same topology as T[i], but different labels
			randT->assembleRandomUnifracTree(groupA, groupB);
			
			if (m->control_pressed) { delete randT;  return 0;  }

			//get wscore of random tree
			EstOutput randomData = weighted.getValues(randT, groupA, groupB);
		
			if (m->control_pressed) { delete randT;  return 0;  }
										
			//save scores
			scores[h].push_back(randomData[0]);
		}
	
		delete randT;
	
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "driver");
		exit(1);
	}
}
/***********************************************************/
void UnifracWeightedCommand::printWeightedFile() {
	try {
		vector<double> data;
		vector<string> tags;
		tags.push_back("Score"); tags.push_back("RandFreq"); tags.push_back("RandCumul");
		
		for(int a = 0; a < numComp; a++) {
			output->initFile(groupComb[a], tags);
			//print each line
			for (map<double,double>::iterator it = validScores.begin(); it != validScores.end(); it++) {
				data.push_back(it->first);  data.push_back(rScoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				output->output(data);
				data.clear();
			} 
			output->resetFile();
		}
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "printWeightedFile");
		exit(1);
	}
}


/***********************************************************/
void UnifracWeightedCommand::printWSummaryFile() {
	try {
		//column headers
		outSum << "Tree#" << '\t' << "Groups" << '\t' << "WScore" << '\t';
		m->mothurOut("Tree#\tGroups\tWScore\t");
		if (random) { outSum << "WSig"; m->mothurOut("WSig"); }
		outSum << endl; m->mothurOutEndLine();
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		int count = 0;
		for (int i = 0; i < T.size(); i++) { 
			for (int j = 0; j < numComp; j++) {
				if (random) {
					if (WScoreSig[count] > (1/(float)iters)) {
						outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << WScoreSig[count] << endl; 
						cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << WScoreSig[count] << endl; 
						m->mothurOutJustToLog(toString(i+1) +"\t" + groupComb[j] +"\t" + toString(utreeScores[count]) +"\t" +  toString(WScoreSig[count]) + "\n");   
					}else{
						outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
						cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
						m->mothurOutJustToLog(toString(i+1) +"\t" + groupComb[j] +"\t" + toString(utreeScores[count]) +"\t<" +  toString((1/float(iters))) + "\n");  
					}
				}else{
					outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << endl; 
					cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count]  << endl; 
					m->mothurOutJustToLog(toString(i+1) +"\t" + groupComb[j] +"\t" + toString(utreeScores[count]) +"\n"); 
				}
				count++;
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
void UnifracWeightedCommand::createPhylipFile() {
	try {
		int count = 0;
		//for each tree
		for (int i = 0; i < T.size(); i++) { 
		
            string phylipFileName;
			map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getSimpleName(treefile);
            variables["[tag]"] = toString(i+1);
            if ((outputForm == "lt") || (outputForm == "square")) {
                variables["[tag2]"] = "weighted.phylip";
                phylipFileName = getOutputFileName("phylip",variables);
                outputNames.push_back(phylipFileName); outputTypes["phylip"].push_back(phylipFileName); 
            }else { //column
                variables["[tag2]"] = "weighted.column";
                phylipFileName = getOutputFileName("column",variables);
                outputNames.push_back(phylipFileName); outputTypes["column"].push_back(phylipFileName); 
            }

			
			ofstream out;
			m->openOutputFile(phylipFileName, out);
			
			if ((outputForm == "lt") || (outputForm == "square")) {
				//output numSeqs
				out << m->getNumGroups() << endl;
			}

			//make matrix with scores in it
			vector< vector<float> > dists;	dists.resize(m->getNumGroups());
			for (int i = 0; i < m->getNumGroups(); i++) {
				dists[i].resize(m->getNumGroups(), 0.0);
			}
			
			//flip it so you can print it
			for (int r=0; r<m->getNumGroups(); r++) { 
				for (int l = 0; l < r; l++) {
					dists[r][l] = utreeScores[count];
					dists[l][r] = utreeScores[count];
					count++;
				}
			}

			//output to file
			for (int r=0; r<m->getNumGroups(); r++) { 
				//output name
				string name = (m->getGroups())[r];
				if (name.length() < 10) { //pad with spaces to make compatible
					while (name.length() < 10) {  name += " ";  }
				}
				
				if (outputForm == "lt") {
					out << name << '\t';
					
					//output distances
					for (int l = 0; l < r; l++) {	out  << dists[r][l] << '\t';  }
					out << endl;
				}else if (outputForm == "square") {
					out << name << '\t';
					
					//output distances
					for (int l = 0; l < m->getNumGroups(); l++) {	out  << dists[r][l] << '\t';  }
					out << endl;
				}else{
					//output distances
					for (int l = 0; l < r; l++) {	
						string otherName = (m->getGroups())[l];
						if (otherName.length() < 10) { //pad with spaces to make compatible
							while (otherName.length() < 10) {  otherName += " ";  }
						}
						
						out  << name << '\t' << otherName << dists[r][l] << endl;  
					}
				}
			}
			out.close();
		}
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "createPhylipFile");
		exit(1);
	}
}
/***********************************************************/
int UnifracWeightedCommand::findIndex(float score, int index) {
	try{
		for (int i = 0; i < rScores[index].size(); i++) {
			if (rScores[index][i] >= score)	{	return i;	}
		}
		return rScores[index].size();
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "findIndex");
		exit(1);
	}
}

/***********************************************************/

void UnifracWeightedCommand::calculateFreqsCumuls() {
	try {
		//clear out old tree values
		rScoreFreq.clear();
		rScoreFreq.resize(numComp);
		rCumul.clear();
		rCumul.resize(numComp);
		validScores.clear();
	
		//calculate frequency
		for (int f = 0; f < numComp; f++) {
			for (int i = 0; i < rScores[f].size(); i++) { //looks like 0,0,1,1,1,2,4,7...  you want to make a map that say rScoreFreq[0] = 2, rScoreFreq[1] = 3...
				validScores[rScores[f][i]] = rScores[f][i];
				map<double,double>::iterator it = rScoreFreq[f].find(rScores[f][i]);
				if (it != rScoreFreq[f].end()) {
					rScoreFreq[f][rScores[f][i]]++;
				}else{
					rScoreFreq[f][rScores[f][i]] = 1;
				}
			}
		}
		
		//calculate rcumul
		for(int a = 0; a < numComp; a++) {
			float rcumul = 1.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (map<double,double>::iterator it = validScores.begin(); it != validScores.end(); it++) {
				//make rscoreFreq map and rCumul
				map<double,double>::iterator it2 = rScoreFreq[a].find(it->first);
				rCumul[a][it->first] = rcumul;
				//get percentage of random trees with that info
				if (it2 != rScoreFreq[a].end()) {  rScoreFreq[a][it->first] /= iters; rcumul-= it2->second;  }
				else { rScoreFreq[a][it->first] = 0.0000; } //no random trees with that score
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "UnifracWeightedCommand", "calculateFreqsCums");
		exit(1);
	}
}
/***********************************************************/





