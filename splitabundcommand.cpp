/*
 *  splitabundcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitabundcommand.h"

//**********************************************************************************************************************
SplitAbundCommand::SplitAbundCommand(string option)  {
	try {
		abort = false;
		wroteRareList = false;
		wroteAbundList = false;
		allLines = 1;
			
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"list","name","group","label","accnos","cutoff","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
						
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}

			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }	
			
			//check for required parameters
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else {  
				groupMap = new GroupMap(groupfile);
				
				int error = groupMap->readMap();
				if (error == 1) { abort = true; }
			}
			
			//do you have all files needed
			if ((listfile == "") && (namefile == "")) { m->mothurOut("You must either a listfile or a namefile for the split.abund command. "); m->mothurOutEndLine(); abort = true;  }
			if ((listfile != "") && (namefile != "")) { m->mothurOut("You must either a listfile or a namefile for the split.abund command, but NOT BOTH. "); m->mothurOutEndLine(); abort = true;  }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = "";  allLines = 1; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			string temp = validParameter.validFile(parameters, "accnos", false);		if (temp == "not found") { temp = "F"; }
			accnos = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);				if (temp == "not found") { temp = "0"; }
			convert(temp, cutoff); 

			if (cutoff == 0) {  m->mothurOut("You must provide a cutoff to qualify what is abundant for the split.abund command. "); m->mothurOutEndLine(); abort = true;  }

		}

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "SplitAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void SplitAbundCommand::help(){
	try {
		m->mothurOut("The split.abund command reads a list or a names file splits the sequences into rare and abundant groups.. \n");
		m->mothurOut("The split.abund command parameters are list, name, cutoff, group, label and accnos.\n");
		m->mothurOut("The list or name parameter is required, and you must provide a cutoff value.\n");
		m->mothurOut("The cutoff parameter is used to qualify what is abundant and rare.\n");
		m->mothurOut("The group parameter allows you to parse a group file into rare and abundant groups.\n");
		m->mothurOut("The label parameter is used to read specific labels in your listfile you want to use.\n");
		m->mothurOut("The accnos parameter allows you to output a .rare.accnos and .abund.accnos files to use with the get.seqs and remove.seqs commands.\n");
		m->mothurOut("The split.abund command should be used in the following format: split.abund(list=yourListFile, group=yourGroupFile, label=yourLabels, cutoff=yourCutoff).\n");
		m->mothurOut("Example: split.abundt(list=abrecovery.fn.list, group=abrecovery.groups, label=0.03, cutoff=2).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
SplitAbundCommand::~SplitAbundCommand(){ 
	if (groupfile != "") {  delete groupMap;  } 
}
//**********************************************************************************************************************
int SplitAbundCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		if (namefile != "") {  split();  }
		else {
		
			//remove old files so you can append later....
			string fileroot = outputDir + getRootName(getSimpleName(listfile));
			remove((fileroot + "rare.list").c_str());
			remove((fileroot + "abund.list").c_str());
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;	
			
			input = new InputData(listfile, "list");
			list = input->getListVector();
			string lastLabel = list->getLabel();
			
			if (m->control_pressed) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->control_pressed) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
				
				if(allLines == 1 || labels.count(list->getLabel()) == 1){
						
						m->mothurOut(list->getLabel()); m->mothurOutEndLine();
						split(list);
											
						processedLabels.insert(list->getLabel());
						userLabels.erase(list->getLabel());
				}
				
				if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
						string saveLabel = list->getLabel();
						
						delete list;
						list = input->getListVector(lastLabel); //get new list vector to process
						
						m->mothurOut(list->getLabel()); m->mothurOutEndLine();
						split(list);
						
						processedLabels.insert(list->getLabel());
						userLabels.erase(list->getLabel());
						
						//restore real lastlabel to save below
						list->setLabel(saveLabel);
				}
				
			
				lastLabel = list->getLabel();
					
				delete list;
				list = input->getListVector(); //get new list vector to process
			}
			
			if (m->control_pressed) { delete input;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
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
			
			if (m->control_pressed) { delete input;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			//run last label if you need to
			if (needToRun == true)  {
				if (list != NULL) {	delete list;	}
				list = input->getListVector(lastLabel); //get new list vector to process
				
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				split(list);		
				
				delete list;
			}
			
			delete input;
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }	return 0;	}
			
			if (wroteAbundList) {  outputNames.push_back(fileroot + "abund.list");		}
			if (wroteRareList)	{  outputNames.push_back(fileroot + "rare.list");		}
		}
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "execute");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::split(ListVector* thisList) {
	try {
	
		SAbundVector* sabund = new SAbundVector();
		*sabund = thisList->getSAbundVector();
		
		//find out how many bins are rare and how many are abundant so you can process the list vector one bin at a time
		// and don't have to store the bins until you are done with the whole vector, this save alot of space.
		int numRareBins = 0;
		for (int i = 0; i <= sabund->getMaxRank(); i++) {
			if (i > cutoff) { break; }
			numRareBins += sabund->get(i);
		}
		int numAbundBins = thisList->getNumBins() - numRareBins;
		delete sabund;
		
		//setup output files
		ofstream outListAbund;
		ofstream outListRare;
		ofstream outGroupRare;
		ofstream outGroupAbund;
		ofstream outAccnosRare;
		ofstream outAccnosAbund;
		
		string fileroot = outputDir + getRootName(getSimpleName(listfile));
		if (numRareBins > 0) {
			wroteRareList = true;
			string listRareName = fileroot + "rare.list";
			openOutputFileAppend(listRareName, outListRare);
			outListRare << thisList->getLabel() << '\t' << numRareBins << '\t';
			
			if (accnos) {
				string accnosName = fileroot + thisList->getLabel() + ".rare.accnos";
				openOutputFile(accnosName, outAccnosRare);
				outputNames.push_back(accnosName);
			}
			
			if (groupfile != "") {
				string groupFileName = outputDir + getRootName(getSimpleName(groupfile)) + thisList->getLabel() + ".rare.group";
				openOutputFile(groupFileName, outGroupRare);
				outputNames.push_back(groupFileName);
			}
		}
		
		if (numAbundBins > 0) {
			wroteAbundList = true;
			string listAbundName = fileroot + "abund.list";
			openOutputFileAppend(listAbundName, outListAbund);
			outListAbund << thisList->getLabel() << '\t' << numAbundBins << '\t';
			
			if (accnos) {
				string accnosName = fileroot + thisList->getLabel() + ".abund.accnos";
				openOutputFile(accnosName, outAccnosAbund);
				outputNames.push_back(accnosName);
			}
			
			if (groupfile != "") {
				string groupFileName = outputDir + getRootName(getSimpleName(groupfile)) + thisList->getLabel() + ".abund.group";
				openOutputFile(groupFileName, outGroupAbund);
				outputNames.push_back(groupFileName);
			}
		}
		
		for (int i = 0; i < thisList->getNumBins(); i++) {
			if (m->control_pressed) { break; }
			
			string bin = list->get(i); 
			
			int size = getNumNames(bin);
			
			if (size <= cutoff) {  outListRare << bin << '\t';  }
			else				{  outListAbund << bin << '\t'; }
			
			if ((groupfile != "") || (accnos)) { //you need to parse the bin...
				vector<string> names;
				splitAtComma(bin, names);  //parses bin into individual sequence names
				
				//parse bin into list of sequences in each group
				for (int j = 0; j < names.size(); j++) {
					
					//write to accnos file
					if (accnos) {
						if (size <= cutoff) {  outAccnosRare << names[j] << endl;  }
						else				{  outAccnosAbund << names[j] << endl; }
					}
					
					//write to groupfile
					if (groupfile != "") {
						string group = groupMap->getGroup(names[j]);
					
						if (group == "not found") {  //error in groupfile so close and remove output file and disregard groupfile
							m->mothurOut(names[j] + " is not in your groupfile. disregarding groupfile."); m->mothurOutEndLine(); 
							delete groupMap; 
							if (numAbundBins > 0) { 
								outGroupAbund.close();
								remove((outputDir + getRootName(getSimpleName(groupfile)) + thisList->getLabel() + ".abund.group").c_str()); 
							}
							if (numRareBins > 0) { 
								outGroupRare.close();
								remove((outputDir + getRootName(getSimpleName(groupfile)) + thisList->getLabel() + ".rare.group").c_str());  
							}
							groupfile = "";
						}else {
							if (size <= cutoff) {  outGroupRare << names[j] << '\t' << group << endl;  }
							else				{  outGroupAbund << names[j] << '\t' << group << endl; }
						}
					}
					
				}//end for names
			}//end if parse
		}//end for list
		
		
		//close files
		if (numRareBins > 0) {	
			outListRare << endl;
			outListRare.close();
			if (accnos) {	outAccnosRare.close();	}
			if (groupfile != "") {  outGroupRare.close();	}
		}
		
		if (numAbundBins > 0) {	
			outListAbund << endl;
			outListAbund.close();
			if (accnos) {	outAccnosAbund.close();	}
			if (groupfile != "") {  outGroupAbund.close();	}
		}

		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "split");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::split() { //namefile
	try {
		//setup output files
		ofstream outNameAbund;
		ofstream outNameRare;
		ofstream outGroupRare;
		ofstream outGroupAbund;
		ofstream outAccnosRare;
		ofstream outAccnosAbund;
		
		bool wroteNameAbund = false;
		bool wroteNameRare = false;
		bool wroteGroupRare = false;
		bool wroteGroupAbund = false;
		bool wroteAccnosRare = false;
		bool wroteAccnosAbund = false;
		
		//prepare output files
		string fileroot = outputDir + getRootName(getSimpleName(namefile));
			
		string nameRareName = fileroot + "rare.names";
		openOutputFile(nameRareName, outNameRare);
		string nameAbundName = fileroot + "abund.names";
		openOutputFile(nameAbundName, outNameAbund);
			
		if (accnos) {
			string accnosName = fileroot + "rare.accnos";
			openOutputFile(accnosName, outAccnosRare);
			
			accnosName = fileroot + "abund.accnos";
			openOutputFile(accnosName, outAccnosAbund);
		}
			
		if (groupfile != "") {
			string groupFileName = outputDir + getRootName(getSimpleName(groupfile))  + ".rare.group";
			openOutputFile(groupFileName, outGroupRare);
			
			groupFileName = outputDir + getRootName(getSimpleName(groupfile))  + ".abund.group";
			openOutputFile(groupFileName, outGroupAbund);
		}
		
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol >> secondCol; gobble(in);
			
			int size = getNumNames(secondCol);
				
			if (size <= cutoff) {  outNameRare << firstCol << '\t' << secondCol << endl;  wroteNameRare = true;  }
			else				{  outNameAbund << firstCol << '\t' << secondCol << endl; wroteNameAbund = true; }

			
			if ((groupfile != "") || (accnos)) { //you need to parse the bin...
				vector<string> names;
				splitAtComma(secondCol, names);  //parses bin into individual sequence names
				
				//parse bin into list of sequences in each group
				for (int j = 0; j < names.size(); j++) {
					
					//write to accnos file
					if (accnos) {
						if (size <= cutoff) {  outAccnosRare << names[j] << endl;  wroteAccnosRare = true; }
						else				{  outAccnosAbund << names[j] << endl; wroteAccnosAbund = true; }
					}
					
					//write to groupfile
					if (groupfile != "") {
						string group = groupMap->getGroup(names[j]);
					
						if (group == "not found") {  //error in groupfile so close and remove output file and disregard groupfile
							m->mothurOut(names[j] + " is not in your groupfile. disregarding groupfile."); m->mothurOutEndLine(); 
							delete groupMap; 
						
							outGroupAbund.close();
							remove((outputDir + getRootName(getSimpleName(groupfile))  + ".abund.group").c_str()); 
							outGroupRare.close();
							remove((outputDir + getRootName(getSimpleName(groupfile)) + ".rare.group").c_str());  
							
							groupfile = "";
							wroteGroupRare = false;
							wroteGroupAbund = false;
						}else {
							if (size <= cutoff) {  outGroupRare << names[j] << '\t' << group << endl;  wroteGroupRare = true; }
							else				{  outGroupAbund << names[j] << '\t' << group << endl; wroteGroupAbund = true; }
						}
					}
					
				}//end for names
			}//end if parse
		}//end while
		
		
		//close files
		in.close();
		outNameRare.close();
		outNameAbund.close();
		if (!wroteNameRare) { remove((fileroot + "rare.names").c_str());  }
		else { outputNames.push_back((fileroot + "rare.names"));  }
		if (!wroteNameAbund) { remove((fileroot + "abund.names").c_str());  }
		else { outputNames.push_back((fileroot + "abund.names"));  }
		
		if (groupfile != "") {  
			outGroupRare.close();	 outGroupAbund.close();
			if (!wroteGroupRare) { remove((outputDir + getRootName(getSimpleName(groupfile))  + ".rare.group").c_str());  }
			else { outputNames.push_back((outputDir + getRootName(getSimpleName(groupfile))  + ".rare.group"));  }
			if (!wroteGroupAbund) { remove((outputDir + getRootName(getSimpleName(groupfile))  + ".abund.group").c_str());  }
			else { outputNames.push_back((outputDir + getRootName(getSimpleName(groupfile))  + ".abund.group"));  }
		}
	
		if (accnos) {	
			outAccnosAbund.close();	outAccnosRare.close();
			if (!wroteAccnosRare) { remove((fileroot + "rare.accnos").c_str());  }
			else { outputNames.push_back((fileroot + "rare.accnos"));  }
			if (!wroteAccnosAbund) { remove((fileroot + "abund.accnos").c_str());  }
			else { outputNames.push_back((fileroot + "abund.accnos"));  }
		}
						
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "split");
		exit(1);
	}
}

/**********************************************************************************************************************/


