/*
 *  parsimonycommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "parsimonycommand.h"


//**********************************************************************************************************************
vector<string> ParsimonyCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none","parsimony-psummary",false,true,true); parameters.push_back(ptree);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter prandom("random", "String", "", "", "", "", "","",false,false); parameters.push_back(prandom);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["parsimony"] = tempOutNames;
        outputTypes["psummary"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParsimonyCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The parsimony command parameters are tree, group, name, count, random, groups, processors and iters.  tree parameter is required unless you have valid current tree file or are using random.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n";
		helpString += "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n";
		helpString += "The parsimony command should be in the following format: parsimony(random=yourOutputFilename, groups=yourGroups, iters=yourIters).\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "Example parsimony(random=out, iters=500).\n";
		helpString += "The default value for random is "" (meaning you want to use the trees in your inputfile, randomtree=out means you just want the random distribution of trees outputted to out.rd_parsimony),\n";
		helpString += "and iters is 1000.  The parsimony command output two files: .parsimony and .psummary their descriptions are in the manual.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParsimonyCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "parsimony") {  pattern = "[filename],parsimony"; } 
        else if (type == "psummary") {  pattern = "[filename],psummary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ParsimonyCommand", "getOutputPattern");
        exit(1);
    }
}
/***********************************************************/
ParsimonyCommand::ParsimonyCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			
			randomtree = validParameter.valid(parameters, "random");		if (randomtree == "not found") { randomtree = ""; }
			
			//are you trying to use parsimony without reading a tree or saying you want random distribution
			if (randomtree == "")  {
				//check for required parameters
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
			}
			
			if (outputdir == ""){  	if (randomtree == "")  { outputdir += util.hasPath(treefile); } }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "";  }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
				
			itersString = validParameter.valid(parameters, "iters");			if (itersString == "not found") { itersString = "1000"; }
			util.mothurConvert(itersString, iters); 
			
			string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
						
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "ParsimonyCommand");
		exit(1);
	}
}
/***********************************************************/
int ParsimonyCommand::execute() {
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        Treenames = util.parseTreeFile(treefile); //extract treenames
		
		//randomtree will tell us if user had their own treefile or if they just want the random distribution
		//user has entered their own tree
		if (randomtree == "") {
			current->setTreeFile(treefile);
			
            TreeReader* reader;
            if (countfile == "") { reader = new TreeReader(treefile, groupfile, namefile); }
            else { reader = new TreeReader(treefile, countfile); }
            T = reader->getTrees();
            ct = T[0]->getCountTable();
            delete reader;
	
			if(outputdir == "") { outputdir += util.hasPath(treefile); }
            map<string, string> variables; 
            variables["[filename]"] = outputdir + util.getSimpleName(treefile) +  ".";
            
			output = new ColumnFile(getOutputFileName("parsimony",variables), itersString);
			outputNames.push_back(getOutputFileName("parsimony",variables));
			outputTypes["parsimony"].push_back(getOutputFileName("parsimony",variables));
				
			sumFile = getOutputFileName("psummary",variables);
			util.openOutputFile(sumFile, outSum);
			outputNames.push_back(sumFile);
			outputTypes["psummary"].push_back(sumFile);
		}else { //user wants random distribution
			getUserInput();
				
			if(outputdir == "") { outputdir += util.hasPath(randomtree); }
			output = new ColumnFile(outputdir+ util.getSimpleName(randomtree), itersString);
			outputNames.push_back(outputdir+ util.getSimpleName(randomtree));
			outputTypes["parsimony"].push_back(outputdir+ util.getSimpleName(randomtree));
		}
			
		//set users groups to analyze
		vector<string> tGroups = ct->getNamesOfGroups();
        //check that groups are valid
        for (int i = 0; i < Groups.size(); i++) {
            if (!util.inUsersGroups(Groups[i], tGroups)) {
                m->mothurOut(Groups[i] + " is not a valid group, and will be disregarded.\n"); 
                // erase the invalid group from userGroups
                Groups.erase(Groups.begin()+i);
                i--;
            }
        }
        if (Groups.size() == 0) { Groups = tGroups;  }
		util.getCombos(groupComb, Groups, numComp);
			
		if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
        
        if (numComp < processors) { m->mothurOut("Reducing processors to " + toString(numComp) + ".\n");  }
		Parsimony pars(Groups);
		counter = 0;
	
		
		
		if (m->getControl_pressed()) { 
			 delete output;
			delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; }
			if (randomtree == "") {  outSum.close();  }
			for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear();
			return 0;
		}
			
		
		//get pscore for users tree
		userData.resize(numComp,0);  //data = AB, AC, BC, ABC.
		randomData.resize(numComp,0);  //data = AB, AC, BC, ABC.
		rscoreFreq.resize(numComp); uscoreFreq.resize(numComp); rCumul.resize(numComp); uCumul.resize(numComp); userTreeScores.resize(numComp); UScoreSig.resize(numComp);
				
		if (randomtree == "") {
			//get pscores for users trees
			for (int i = 0; i < T.size(); i++) {
				userData = pars.getValues(T[i], processors, outputdir);  //data = AB, AC, BC, ABC.
				
				if (m->getControl_pressed()) { 
					 delete output;
					delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; }
					if (randomtree == "") {  outSum.close();  }
					for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear();
					return 0;
				}

				//output scores for each combination
				for(int k = 0; k < numComp; k++) {

					//update uscoreFreq
					map<int,double>::iterator it = uscoreFreq[k].find(userData[k]);
					if (it == uscoreFreq[k].end()) {//new score
						uscoreFreq[k][userData[k]] = 1;
					}else{ uscoreFreq[k][userData[k]]++; }
					
					//add users score to valid scores
					validScores[userData[k]] = userData[k];
					
					//save score for summary file
					userTreeScores[k].push_back(userData[k]);
				}
			}
			
            Utils* stableRandom = new Utils();
			//get pscores for random trees
			for (int j = 0; j < iters; j++) {
								
				//create new tree with same num nodes and leaves as users
				randT = new Tree(ct, Treenames);

				//create random relationships between nodes
				randT->assembleRandomTree(stableRandom);

				//get pscore of random tree
				randomData = pars.getValues(randT, processors, outputdir);
				
				if (m->getControl_pressed()) { 
                      delete output; delete randT; delete stableRandom;
					if (randomtree == "") {  outSum.close();  }
					for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear();
					delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; }
					return 0;
				}
					
				for(int r = 0; r < numComp; r++) {
					//add trees pscore to map of scores
					map<int,double>::iterator it = rscoreFreq[r].find(randomData[r]);
					if (it != rscoreFreq[r].end()) {//already have that score
						rscoreFreq[r][randomData[r]]++;
					}else{//first time we have seen this score
						rscoreFreq[r][randomData[r]] = 1;
					}
			
					//add randoms score to validscores
					validScores[randomData[r]] = randomData[r];
				}
				
								
				delete randT;
			}
            delete stableRandom;
		}else {
            Utils* stableRandom = new Utils();
			//get pscores for random trees
			for (int j = 0; j < iters; j++) {
								
				//create new tree with same num nodes and leaves as users
				randT = new Tree(ct, Treenames);
				//create random relationships between nodes

				randT->assembleRandomTree(stableRandom);
				
				if (m->getControl_pressed()) { 
					 delete output; delete randT; delete ct;  delete stableRandom;
					for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear(); return 0;
				}

				//get pscore of random tree
				randomData = pars.getValues(randT, processors, outputdir);
				
				if (m->getControl_pressed()) { 
					 delete output; delete randT; delete ct;
					for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear(); return 0;
				}
			
				for(int r = 0; r < numComp; r++) {
					//add trees pscore to map of scores
					map<int,double>::iterator it = rscoreFreq[r].find(randomData[r]);
					if (it != rscoreFreq[r].end()) {//already have that score
						rscoreFreq[r][randomData[r]]++;
					}else{//first time we have seen this score
						rscoreFreq[r][randomData[r]] = 1;
					}
			
					//add randoms score to validscores
					validScores[randomData[r]] = randomData[r];
				}
				
				
				delete randT;
			}
            delete stableRandom;
		}

		for(int a = 0; a < numComp; a++) {
			float rcumul = 0.0000;
			float ucumul = 0.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (map<int,double>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				if (randomtree == "") {
					map<int,double>::iterator it2 = uscoreFreq[a].find(it->first);
					//user data has that score 
					if (it2 != uscoreFreq[a].end()) { uscoreFreq[a][it->first] /= T.size(); ucumul+= it2->second;  }
					else { uscoreFreq[a][it->first] = 0.0000; } //no user trees with that score
					//make uCumul map
					uCumul[a][it->first] = ucumul;
				}
			
				//make rscoreFreq map and rCumul
				map<int,double>::iterator it2 = rscoreFreq[a].find(it->first);
				//get percentage of random trees with that info
				if (it2 != rscoreFreq[a].end()) {  rscoreFreq[a][it->first] /= iters; rcumul+= it2->second;  }
				else { rscoreFreq[a][it->first] = 0.0000; } //no random trees with that score
				rCumul[a][it->first] = rcumul;
			}
			
			//find the signifigance of each user trees score when compared to the random trees and save for printing the summary file
			for (int h = 0; h < userTreeScores[a].size(); h++) {
				UScoreSig[a].push_back(rCumul[a][userTreeScores[a][h]]);
			}
		}
		
		if (m->getControl_pressed()) { 
				delete output;
				delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; }
				if (randomtree == "") {  outSum.close();  }
				for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear();
				return 0;
		}
				
		printParsimonyFile();
		if (randomtree == "") { printUSummaryFile(); }
				
        delete output; delete ct; for (int i = 0; i < T.size(); i++) { delete T[i]; }
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } outputTypes.clear(); return 0;}
		
		m->mothurOut("\nOutput File Names: \n");
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "execute");
		exit(1);
	}
}

/***********************************************************/
void ParsimonyCommand::printParsimonyFile() {
	try {
		vector<double> data;
		vector<string> tags;
		
		if (randomtree == "") { tags.push_back("Score"); tags.push_back("UserFreq"); tags.push_back("UserCumul"); tags.push_back("RandFreq"); tags.push_back("RandCumul"); }
		else { tags.push_back("Score"); tags.push_back("RandFreq"); tags.push_back("RandCumul"); }

		for(int a = 0; a < numComp; a++) {
			output->setLabelName(groupComb[a], tags);
			//print each line
			for (map<int,double>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				if (randomtree == "") {
					data.push_back(it->first);  data.push_back(uscoreFreq[a][it->first]); data.push_back(uCumul[a][it->first]); data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				}else{
					data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				}
				output->updateOutput(data);
				data.clear();
			} 
			output->resetFile();
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "printParsimonyFile");
		exit(1);
	}
}
/***********************************************************/
int ParsimonyCommand::printUSummaryFile() {
	try {
		//column headers
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "ParsScore" << '\t' << "ParsSig" <<  endl;
		m->mothurOut("Tree#\tGroups\tParsScore\tParsSig\n");
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		for (int i = 0; i< T.size(); i++) {
			for(int a = 0; a < numComp; a++) {
				if (m->getControl_pressed()) {  outSum.close(); return 0; }
				if (UScoreSig[a][i] > (1/(float)iters)) {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[a]  << '\t' << userTreeScores[a][i] << setprecision(itersString.length()) << '\t' << UScoreSig[a][i] << endl;
					cout << setprecision(6) << i+1 << '\t' << groupComb[a]  << '\t' << userTreeScores[a][i] << setprecision(itersString.length()) << '\t' << UScoreSig[a][i] << endl;
					m->mothurOutJustToLog(toString(i+1) + "\t" + groupComb[a] + "\t" + toString(userTreeScores[a][i]) + "\t" + toString(UScoreSig[a][i])); m->mothurOutEndLine();
				}else {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << userTreeScores[a][i] << setprecision(itersString.length())  << '\t' << "<" << (1/float(iters)) << endl;
					cout << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << userTreeScores[a][i] << setprecision(itersString.length()) << '\t' << "<" << (1/float(iters)) << endl;
					m->mothurOutJustToLog(toString(i+1) + "\t" + groupComb[a] + "\t" + toString(userTreeScores[a][i]) + "\t" + toString((1/float(iters)))); m->mothurOutEndLine();
				}
			}
		}
		
		outSum.close();
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "printUSummaryFile");
		exit(1);
	}
}

/***********************************************************/
void ParsimonyCommand::getUserInput() {
	try {
		//create treemap
		ct = new CountTable();

		m->mothurOut("Please enter the number of groups you would like to analyze: ");
		cin >> numGroups;
		m->mothurOutJustToLog(toString(numGroups)); m->mothurOutEndLine();
				
		int num, count;
		count = 1;
		numEachGroup.resize(numGroups, 0);  
		
        set<string> nameMap;
        map<string, string> groupMap;
        set<string> gps;
                
		for (int i = 1; i <= numGroups; i++) {
			m->mothurOut("Please enter the number of sequences in group " + toString(i) +  ": ");
			cin >> num;
			m->mothurOutJustToLog(toString(num)); m->mothurOutEndLine();
			
            gps.insert(toString(i));
            
			//set tmaps namesOfSeqs
			for (int j = 0; j < num; j++) {
				groupMap[toString(count)] = toString(i);
				nameMap.insert(toString(count));
				count++;
			}
		}
		ct->createTable(nameMap, groupMap, gps);
        
		//clears buffer so next command doesn't have error
		string s;	
		getline(cin, s);
		
		Treenames = ct->getNamesOfSeqs();
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getUserInput");
		exit(1);
	}
}
/***********************************************************/


