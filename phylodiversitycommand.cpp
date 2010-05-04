/*
 *  phylodiversitycommand.cpp
 *  Mothur
 *
 *  Created by westcott on 4/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "phylodiversitycommand.h"
#include "phylodiversity.h"

//**********************************************************************************************************************
PhyloDiversityCommand::PhyloDiversityCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"freq","rarefy","iters","groups","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			if (globaldata->gTree.size() == 0) {//no trees were read
				m->mothurOut("You must execute the read.tree command, before you may execute the phylo.diversity command."); m->mothurOutEndLine(); abort = true;  }

			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "rarefy", false);			if (temp == "not found") { temp = "F"; }
			rarefy = isTrue(temp);
			if (!rarefy) { iters = 1;  }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; Groups = globaldata->gTreemap->namesOfGroups;  globaldata->Groups = Groups;  }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloDiversityCommand", "PhyloDiversityCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

void PhyloDiversityCommand::help(){
	try {


	}
	catch(exception& e) {
		m->errorOut(e, "PhyloDiversityCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

PhyloDiversityCommand::~PhyloDiversityCommand(){}

//**********************************************************************************************************************

int PhyloDiversityCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//incase the user had some mismatches between the tree and group files we don't want group xxx to be analyzed
		for (int i = 0; i < globaldata->Groups.size(); i++) { if (globaldata->Groups[i] == "xxx") { globaldata->Groups.erase(globaldata->Groups.begin()+i);  break; }  }
		 
		vector<string> outputNames;
		
		//diversity calculator
		PhyloDiversity phylo(globaldata->gTreemap);
		
		vector<Tree*> trees = globaldata->gTree;
		
		//for each of the users trees
		for(int i = 0; i < trees.size(); i++) {
		
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} return 0; }
			
			phylo.setTotalGroupBranchLengths(trees[i]);
			
			string outFile = outputDir + getRootName(getSimpleName(globaldata->getTreeFile()))  + toString(i+1) + ".phylo.diversity";
			if (rarefy) { outFile += ".rarefaction"; }
			outputNames.push_back(outFile);
			
			int numLeafNodes = trees[i]->getNumLeaves();
			
			//create a vector containing indexes of leaf nodes, randomize it, select nodes to send to calculator
			vector<int> randomLeaf;
			for (int j = 0; j < numLeafNodes; j++) {  
				if (inUsersGroups(trees[i]->tree[j].getGroup(), globaldata->Groups) == true) { //is this a node from the group the user selected.
					randomLeaf.push_back(j); 
				}
			}
			
			numLeafNodes = randomLeaf.size();  //reset the number of leaf nodes you are using 
			
			//each group, each sampling, if no rarefy iters = 1;
			vector< vector<float> > diversity;
			diversity.resize(globaldata->Groups.size());
			
			//initialize sampling spots
			vector<int> numSampledList;
			for(int k = 0; k < numLeafNodes; k++){  if((k == 0) || (k+1) % freq == 0){  numSampledList.push_back(k); }   }
			if(numLeafNodes % freq != 0){	numSampledList.push_back(numLeafNodes);   }
			
			//initialize diversity
			for (int j = 0; j < diversity.size(); j++) {   diversity[j].resize(numSampledList.size(), 0.0);  }  //			10sampled	20 sampled ...
																												//groupA		0.0			0.0
																											//then for each iter you add to score and then when printing divide by iters to get average
			for (int l = 0; l < iters; l++) {
				random_shuffle(randomLeaf.begin(), randomLeaf.end());
		
				vector<int> leavesSampled;
				EstOutput data;
				int count = 0;
				for(int k = 0; k < numLeafNodes; k++){
						
					if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} return 0; }
					
					leavesSampled.push_back(randomLeaf[k]);
						
					if((k == 0) || (k+1) % freq == 0){ //ready to calc?
						
						data = phylo.getValues(trees[i], leavesSampled);
						
						//datas results are in the same order as globaldatas groups
						for (int h = 0; h < data.size(); h++) {  diversity[h][count] += data[h];  }
						
						count++;
					}
				}
		
				if(numLeafNodes % freq != 0){	
					
					data = phylo.getValues(trees[i], leavesSampled);
					
					//datas results are in the same order as globaldatas groups
					for (int h = 0; h < data.size(); h++) {  diversity[h][count] += data[h];  }
				}
			}
			
			printData(numSampledList, diversity, outFile);

		}
		
	
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} return 0; }

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloDiversityCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

void PhyloDiversityCommand::printData(vector<int>& num, vector< vector<float> >& div, string file){
	try {
		ofstream out;
		openOutputFile(file, out);
		
		out << "numSampled\t";
		for (int i = 0; i < globaldata->Groups.size(); i++) { out << globaldata->Groups[i] << '\t';  }
		out << endl;
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		for (int i = 0; i < num.size(); i++) {	
			if (i == (num.size()-1)) {  out << num[i] << '\t';  }
			else {  out << (num[i]+1) << '\t';  }
			
			for (int j = 0; j < div.size(); j++) {
				float score = div[j][i] / (float)iters;
				out << setprecision(6) << score << '\t';
			}
			out << endl;
		}
		
		out.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloDiversityCommand", "printData");
		exit(1);
	}
}

//**********************************************************************************************************************
