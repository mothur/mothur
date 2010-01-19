/*
 *  otuhierarchycommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "otuhierarchycommand.h"

//**********************************************************************************************************************
OtuHierarchyCommand::OtuHierarchyCommand(string option){
	try {
		abort = false;
		//allow user to run help
		if(option == "help") {  help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"list","label"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			listFile = validParameter.validFile(parameters, "list", true);
			if (listFile == "not found") { mothurOut("list is a required parameter for the otu.hierarchy command."); mothurOutEndLine(); abort = true; }
			else if (listFile == "not open") { abort = true; }	

			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { mothurOut("label is a required parameter for the otu.hierarchy command."); mothurOutEndLine(); abort = true; }
			else { 
				splitAtDash(label, labels);
				if (labels.size() != 2) { mothurOut("You must provide 2 labels."); mothurOutEndLine(); abort = true; }
			}				
		}
		
	}
	catch(exception& e) {
		errorOut(e, "OtuHierarchyCommand", "OtuHierarchyCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

void OtuHierarchyCommand::help(){
	try {
		mothurOut("The otu.hierarchy command is used to see how otus relate at two distances. \n");
		mothurOut("The otu.hierarchy command parameters are list and label.  Both parameters are required. \n");
		mothurOut("The otu.hierarchy command should be in the following format: \n");
		mothurOut("otu.hierarchy(list=yourListFile, label=yourLabels).\n");
		mothurOut("Example otu.hierarchy(list=amazon.fn.list, label=0.01-0.03).\n");
		mothurOut("The otu.hierarchy command outputs a .otu.hierarchy file which is described on the wiki.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListFile).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "OtuHierarchyCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

OtuHierarchyCommand::~OtuHierarchyCommand(){}

//**********************************************************************************************************************

int OtuHierarchyCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//get listvectors that correspond to labels requested, (or use smart distancing to get closest listvector)
		vector<ListVector> lists = getListVectors();
		
		//determine which is little and which is big, putting little first
		if (lists.size() == 2) {
			//if big is first swap them
			if (lists[0].getNumBins() < lists[1].getNumBins()) {
				reverse(lists.begin(), lists.end());
			}
		}else{
			mothurOut("error getting listvectors, unable to read 2 different vectors, check your label inputs."); mothurOutEndLine(); return 0;
		}
		
		//map sequences to bin number in the "little" otu
		map<string, int> littleBins; 
		for (int i = 0; i < lists[0].getNumBins(); i++) {
			string names = lists[0].get(i); 
			
			//parse bin
			while (names.find_first_of(',') != -1) { 
				string name = names.substr(0,names.find_first_of(','));
				names = names.substr(names.find_first_of(',')+1, names.length());
				littleBins[name] = i;  
			}
			
			//get last name
			littleBins[names] = i;
		}
		
		ofstream out;
		string outputFileName = getRootName(listFile) + lists[0].getLabel() + "-" + lists[1].getLabel() + ".otu.hierarchy";
		openOutputFile(outputFileName, out);
		
		//go through each bin in "big" otu and output the bins in "little" otu which created it
		for (int i = 0; i < lists[1].getNumBins(); i++) {
		
			string names = lists[1].get(i);
			
			//output column 1
			out << names << '\t';
			
			map<int, int> bins; //bin numbers in little that are in this bin in big
			map<int, int>::iterator it;
			
			//parse bin
			while (names.find_first_of(',') != -1) { 
				string name = names.substr(0,names.find_first_of(','));
				names = names.substr(names.find_first_of(',')+1, names.length());
				bins[littleBins[name]] = littleBins[name];  
			}
			
			//get last name
			bins[littleBins[names]] = littleBins[names]; 
			
			string col2 = "";
			for (it = bins.begin(); it != bins.end(); it++) {
				col2 += lists[0].get(it->first) + "\t";
			}
			
			//output column 2
			out << col2 << endl;
		}
		
		out.close();
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "OtuHierarchyCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
//returns a vector of listVectors where "little" vector is first
vector<ListVector> OtuHierarchyCommand::getListVectors() {
	try {
		
		int pos; //to use in smart distancing, position of last read in file
		int lastPos;
		vector<ListVector> lists;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		//open file
		ifstream in;
		openInputFile(listFile, in);
		
		//get first list vector in file
		ListVector* list = NULL;
		string lastLabel = "";
		if (!in.eof())	{
			pos = in.tellg();
			lastPos = pos;
			list = new ListVector(in);  
			gobble(in);
			lastLabel = list->getLabel();
		}
		
		while ((list != NULL) && (userLabels.size() != 0)) {
			
			//is this a listvector that we want?
			if(labels.count(list->getLabel()) == 1){
				
				//make copy of listvector
				ListVector temp(*list);
				lists.push_back(temp);
			
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
		
			//you have a label the user want that is smaller than this label and the last label has not already been processed 
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				int savePos = in.tellg();
				
				//get smart distance line
				delete list;
				in.seekg(lastPos);
				if (!in.eof())	{	
					list = new ListVector(in);  
				}else { list = NULL; }
				
				//make copy of listvector
				ListVector temp(*list);
				lists.push_back(temp);
					
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());					
										
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
				in.seekg(savePos);
			}
			
			lastLabel = list->getLabel();
			lastPos = pos;
			
			//get next line
			delete list;
			if (!in.eof())	{	
				pos = in.tellg();
				list = new ListVector(in);  
				gobble(in);
			}else { list = NULL; }
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
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			
			in.seekg(lastPos);
			if (!in.eof())	{	
				list = new ListVector(in); 
				
				//make copy of listvector
				ListVector temp(*list);
				lists.push_back(temp);
				
				delete list;
			}
		}
		
		
		return lists;
	}
	catch(exception& e) {
		errorOut(e, "OtuHierarchyCommand", "getListVectors");
		exit(1);
	}
}

//**********************************************************************************************************************





