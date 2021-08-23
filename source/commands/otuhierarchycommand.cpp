/*
 *  otuhierarchycommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "otuhierarchycommand.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> OtuHierarchyCommand::setParameters(){	
	try {
		CommandParameter poutput("output", "Multiple", "name-otulabel", "name", "", "", "","",false,false); parameters.push_back(poutput);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","otuheirarchy",false,true,true); parameters.push_back(plist);
        CommandParameter palist("asvlist", "InputTypes", "", "", "none", "none", "none","otuheirarchy",false,true,true); parameters.push_back(palist);
        CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","constaxonomy",false,true,true); parameters.push_back(ptaxonomy);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; asv = false;
        
        vector<string> tempOutNames;
        outputTypes["otuheirarchy"] = tempOutNames;
        outputTypes["asvconstaxonomy"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string OtuHierarchyCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The otu.hierarchy command is used to see how otus relate at two distances, or to see how ASVs relate to OTUs. \n";
		helpString += "The otu.hierarchy command parameters are list, asvlist, count, taxonomy, label and output.  list and label parameters are required for relating OTUs at different distances. asvlist, list, taxonomy and count are required for ASV to OTU relation. \n";
		helpString += "The output parameter allows you to output the names of the sequence in the OTUs or the OTU labels. Options are name and otulabel, default is name. \n";
		helpString += "The otu.hierarchy command should be in the following format: \n";
        helpString += "otu.hierarchy(list=yourListFile, asvlist=yourAsvListFile, taxonomy=yourTaxonomyFile, count=yourCountFile, label=yourLabels).\n";
		helpString += "otu.hierarchy(list=yourListFile, label=yourLabels).\n";
		helpString += "Example otu.hierarchy(list=amazon.fn.list, label=0.01-0.03).\n";
		helpString += "The otu.hierarchy command outputs a .otu.hierarchy file which is described on the wiki.\n";
		
        getCommonQuestions();
        
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string OtuHierarchyCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
       // string issue = "...template is not aligned, aborting. What do I do?"; issues.push_back(issue);
        //string ianswer = "\tMothur requires the reference file to be aligned to generate aligned sequences. You can download mothur's aligned silva references here, https://mothur.org/wiki/Silva_reference_files. For ITS sequences, see 'how to' below.\n"; ianswers.push_back(ianswer);
        
        //issue = "...xxx of your sequences generated alignments that eliminated too many bases... What does this mean?"; issues.push_back(issue);
        //ianswer = "\tBy default, mothur will align the reverse compliment of your sequences when the alignment process removes more than 50% of the bases indicating the read may be flipped. This process assembles the best possible alignment, and downstream analysis will remove any poor quality reads remaining.\n"; ianswers.push_back(ianswer);
        
        
        string howto = "How do I find the OTUs and taxonomies my ASVs are clustered in?"; howtos.push_back(howto);
        string hanswer = "\tYou can use the otu.hierarchy command to create a *.cons.taxonomy file. The first column is the ASVLabel, the second column is the abundance of the ASV, and the third column is the ASVs taxonomy with the OTULabel appended.\n\nmothur > otu.hierarchy(list=final.opti_mcc.list, asvlist=final.asv.list, taxonomy=final.taxonomy, count=final.count_table)\n"; hanswers.push_back(hanswer);
        
       // howto = "How do I create a custom reference for the region I am studying?"; howtos.push_back(howto);
       // hanswer = "\tYou can tailor your reference using this method: http://blog.mothur.org/2016/07/07/Customization-for-your-region/.\n"; hanswers.push_back(hanswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);

        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "OtuHierarchyCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
string OtuHierarchyCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "otuheirarchy") {  pattern = "[filename],[distance1],[tag],[distance2],otu.hierarchy"; }
        if (type == "asvconstaxonomy") {  pattern = "[filename],[tag],asv.cons.taxonomy"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "OtuHierarchyCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
OtuHierarchyCommand::OtuHierarchyCommand(string option) : Command() {
	try {
		if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
            OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			listFile = validParameter.validFile(parameters, "list");
			if (listFile == "not found") { 
				listFile = current->getListFile(); 
				if (listFile != "") {  m->mothurOut("Using " + listFile + " as input file for the list parameter.\n");  }
				else { 
					m->mothurOut("No valid current list file. You must provide a list file.\n");  
					abort = true;
				}
			}else if (listFile == "not open") { abort = true; }	
			else { current->setListFile(listFile); }
            
            asvlistFile = validParameter.validFile(parameters, "asvlist");
            if (asvlistFile == "not found") { asvlistFile = ""; }
            else if (asvlistFile == "not open") { asvlistFile = ""; abort = true; }
            else { asv = true; }
			
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }
            else { current->setCountFile(countfile); }
            
            taxfile = validParameter.validFile(parameters, "taxonomy");
            if (taxfile == "not found") { taxfile = "";  }
            else if (taxfile == "not open") { taxfile = "";  abort = true; }
            else { current->setTaxonomyFile(taxfile); }
            
            if (outputdir == ""){	 outputdir += util.hasPath(listFile);  }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") {
                if (!asv) { m->mothurOut("[ERROR]: label is a required parameter for the otu.hierarchy command, please correct.\n");  abort = true; }
                else {  m->mothurOut("\nNo label provided, I will use the first label in the list file.\n"); }
            
            
            }else {
                util.splitAtDash(label, mylabels);
                if (!asv) {    if (mylabels.size() != 2) { m->mothurOut("You must provide 2 labels.\n");  abort = true;  }  }
			}	
			
			output = validParameter.valid(parameters, "output");			if (output == "not found") { output = "name"; }
			if ((output != "name") && (output != "otulabel")) { m->mothurOut("output options are name and otulabel. I will use name.\n");  output = "name"; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "OtuHierarchyCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

int OtuHierarchyCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        if (asv)    { processASV();             }
        else        { processHierarchy();       }
        
        if (m->getControl_pressed()) { outputTypes.clear();   for (int j = 0; j < outputNames.size(); j++) {    util.mothurRemove(outputNames[j]);    }  return 0; }
		
        m->mothurOut("\nOutput File Names:\n");
        for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i]); m->mothurOutEndLine();    }
        m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void OtuHierarchyCommand::processASV() {
    try {
        set<string> labels;
        if (mylabels.size() != 0) { labels.insert(*mylabels.begin()); }
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
    
        //read otu list file
        InputData inputOTU(listFile, "list", nullVector);
        ListVector* list = util.getNextList(inputOTU, false, userLabels, processedLabels, lastLabel);
        string otuListLable = list->getLabel();
      
        //read taxonomy file
        map<string, string> taxMap; map<string, string>::iterator itTax;
        util.readTax(taxfile, taxMap, true);
        
        //append OTU label to taxonomy
        for (int i = 0; i < list->getNumBins(); i++) {
            
            if (m->getControl_pressed()) {  return;  }
            
            string binnames = list->get(i);
            string otuLabel = list->getOTUName(i);
            
            //parse names in bin
            vector<string> names; util.splitAtComma(binnames, names);
            
            for (int j = 0; j < names.size(); j++) {
                
                itTax = taxMap.find(names[j]);
                
                if (itTax != taxMap.end()) {
                    itTax->second += otuLabel + ";";
                    
                }else{ m->mothurOut("\n[ERROR]: " + names[j] + " is missing from your taxonomy file, please correct.\n"); m->setControl_pressed(true); }
            }
        }
        delete list;
        
        //add redundant counts
        CountTable ct; bool hasCount = false;
        if (countfile != "") { ct.readTable(countfile, true, false); hasCount = true; }
        
        if (m->getControl_pressed()) {  return;  }
        
        //read asvlist file
        labels.clear(); processedLabels.clear(); lastLabel = "";
        userLabels = labels;
        
        InputData input(asvlistFile, "list", nullVector);
        ListVector* asvlist = util.getNextList(input, false, userLabels, processedLabels, lastLabel);
        string asvLabel = asvlist->getLabel();
        
        if (m->getControl_pressed()) {  return;  }
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(listFile));
        variables["[tag]"] = asvLabel + "-" + otuListLable;
        string outputFileName = getOutputFileName("asvconstaxonomy",variables);
        outputNames.push_back(outputFileName); outputTypes["asvconstaxonomy"].push_back(outputFileName);
        ofstream out; util.openOutputFile(outputFileName, out);
        
        out << "ASVLabel\tASV_Abundance\tTaxonomy_OTULabel\n";
        
        for (int i = 0; i < asvlist->getNumBins(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            string binnames = asvlist->get(i);
            string asvOtuLabel = asvlist->getOTUName(i);
            
            //parse names in bin
            vector<string> names; util.splitAtComma(binnames, names);
            
            for (int j = 0; j < names.size(); j++) {
                
                itTax = taxMap.find(names[j]);
                
                int abund = 1;
                if (itTax != taxMap.end()) {
                    if (hasCount) { abund = ct.getNumSeqs(names[j]); }
                    
                    out << asvOtuLabel << '\t' << abund << '\t' << itTax->second << endl;
                    
                }else{ m->mothurOut("\n[ERROR]: " + names[j] + " is missing from your taxonomy file, please correct.\n"); m->setControl_pressed(true); }
            }
        }
        out.close();
        
        delete asvlist;
    }
    catch(exception& e) {
        m->errorOut(e, "OtuHierarchyCommand", "processASV");
        exit(1);
    }
}
//**********************************************************************************************************************
void OtuHierarchyCommand::processHierarchy() {
    try {
        //get listvectors that correspond to labels requested, (or use smart distancing to get closest listvector)
        vector< vector<string> > lists = getListVectors();
        
        if (m->getControl_pressed()) { return; }
        
        //determine which is little and which is big, putting little first
        if (lists.size() == 4) {
            //if big is first swap them
            if (lists[0].size() < lists[2].size()) {
                vector< vector<string> > tempLists;
                tempLists.push_back(lists[2]);
                tempLists.push_back(lists[3]);
                tempLists.push_back(lists[0]);
                tempLists.push_back(lists[1]);
                lists = tempLists;
                string tempLabel = list2Label;
                list2Label = list1Label;
                list1Label = tempLabel;
            }
        }else{  m->mothurOut("[ERROR]: error getting listvectors, unable to read 2 different vectors, check your label inputs.\n");  return; }
        
        //map sequences to bin number in the "little" otu
        map<string, int> littleBins;
        vector<string> binLabels0 = lists[0];
        for (int i = 0; i < lists[0].size(); i++) {
        
            if (m->getControl_pressed()) {  return; }
            string bin = lists[1][i];
            vector<string> names; util.splitAtComma(bin, names);
            for (int j = 0; j < names.size(); j++) { littleBins[names[j]] = i; }
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(listFile));
        variables["[distance1]"] = list1Label;
        variables["[tag]"] = "-";
        variables["[distance2]"] = list2Label;
        string outputFileName = getOutputFileName("otuheirarchy",variables);
        outputNames.push_back(outputFileName); outputTypes["otuheirarchy"].push_back(outputFileName);
        ofstream out; util.openOutputFile(outputFileName, out);
        
        //go through each bin in "big" otu and output the bins in "little" otu which created it
        vector<string> binLabels1 = lists[2];
        for (int i = 0; i < lists[2].size(); i++) {
        
            if (m->getControl_pressed()) {  break; }
            
            string binnames = lists[3][i];
            vector<string> names; util.splitAtComma(binnames, names);
            
            //output column 1
            if (output == "name")    {   out << binnames << '\t';    }
            else                    {    out << binLabels1[i] << '\t';        }
            
            map<int, int> bins; //bin numbers in little that are in this bin in big
            map<int, int>::iterator it;
            
            //parse bin
            for (int j = 0; j < names.size(); j++) { bins[littleBins[names[j]]] = littleBins[names[j]];   }
            
            string col2 = "";
            for (it = bins.begin(); it != bins.end(); it++) {
                if (output == "name")    {   col2 += lists[1][it->first] + "\t";    }
                else                    {    col2 += binLabels0[it->first] + "\t";        }
            }
            
            //output column 2
            out << col2 << endl;
        }
        
        out.close();
        
    }
    catch(exception& e) {
        m->errorOut(e, "OtuHierarchyCommand", "processHierarchy");
        exit(1);
    }
}
//**********************************************************************************************************************
//returns a vector of listVectors where "little" vector is first
vector< vector<string> > OtuHierarchyCommand::getListVectors() { //return value [0] -> otulabelsFirstLabel [1] -> binsFirstLabel [2] -> otulabelsSecondLabel [3] -> binsSecondLabel
	try {
		vector< vector<string> > lists;
        
        int count = 0;
        for (set<string>::iterator it = mylabels.begin(); it != mylabels.end(); it++) {
            string realLabel;
            vector< vector<string> > thisList = getListVector(*it, realLabel);
            
            if (m->getControl_pressed()) {  return lists; }
            
            for (int i = 0; i < thisList.size(); i++) { lists.push_back(thisList[i]); }
            
            if (count == 0) {  list1Label = realLabel; count++; }
            else {  list2Label = realLabel; }
        }
        
        return lists;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "getListVectors");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<string> > OtuHierarchyCommand::getListVector(string label, string& realLabel){ //return value [0] -> otulabels [1] -> bins
	try {
        vector< vector<string> > myList;
        
        InputData input(listFile, "list", nullVector);
        set<string> labels; labels.insert(label);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        ListVector* list = util.getNextList(input, false, userLabels, processedLabels, lastLabel);
               
        if (list != NULL) {
                   
            //at this point the list vector has the right distance
            vector<string> bins, listlabels;
            for (int i = 0; i < list->getNumBins(); i++) {
                if (m->getControl_pressed()) {  return myList;  }
                bins.push_back(list->get(i));
                listlabels.push_back(list->getOTUName(i));
            }
            myList.push_back(listlabels);
            myList.push_back(bins);
            realLabel = list->getLabel();
            
            delete list;
        }
        
		return myList;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "getListVector");
		exit(1);
	}
}

//**********************************************************************************************************************





