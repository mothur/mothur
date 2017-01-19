/*
 *  summarytaxcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/23/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "summarytaxcommand.h"
#include "phylosummary.h"

//**********************************************************************************************************************
vector<string> SummaryTaxCommand::setParameters(){	
	try {
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(ptaxonomy);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter prelabund("relabund", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prelabund);
        CommandParameter poutput("output", "Multiple", "simple-detail", "detail", "", "", "","",false,false, true); parameters.push_back(poutput);
        CommandParameter pthreshold("threshold", "Number", "", "0", "", "", "","",false,true); parameters.push_back(pthreshold);
        CommandParameter pprintlevel("printlevel", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pprintlevel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryTaxCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The summary.tax command reads a taxonomy file and an optional name file, and summarizes the taxonomy information.\n";
		helpString += "The summary.tax command parameters are taxonomy, count, group, name and relabund. taxonomy is required, unless you have a valid current taxonomy file.\n";
		helpString += "The name parameter allows you to enter a name file associated with your taxonomy file. \n";
		helpString += "The group parameter allows you add a group file so you can have the summary totals broken up by group.\n";
        helpString += "The count parameter allows you add a count file so you can have the summary totals broken up by group.\n";
        helpString += "The threshold parameter allows you to specify a cutoff for the taxonomy file that is being inputted. Once the classification falls below the threshold the mothur will refer to it as unclassified when calculating the concensus.  This feature is similar to adjusting the cutoff in classify.seqs. Default=0.\n";
        helpString += "The output parameter allows you to specify format of your summary file. Options are simple and detail. The default is detail.\n";
        helpString += "The printlevel parameter allows you to specify taxlevel of your summary file to print to. Options are 1 to the maz level in the file.  The default is -1, meaning max level.  If you select a level greater than the level your sequences classify to, mothur will print to the level your max level. \n";
        helpString += "The relabund parameter allows you to indicate you want the summary file values to be relative abundances rather than raw abundances. Default=F. \n";
		helpString += "The summary.tax command should be in the following format: \n";
		helpString += "summary.tax(taxonomy=yourTaxonomyFile) \n";
		helpString += "Note: No spaces between parameter labels (i.e. taxonomy), '=' and parameters (i.e.yourTaxonomyFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryTaxCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],tax.summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SummaryTaxCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SummaryTaxCommand::SummaryTaxCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "SummaryTaxCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SummaryTaxCommand::SummaryTaxCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}

			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//check for required parameters
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not open") { abort = true; }
			else if (taxfile == "not found") { 				
				taxfile = m->getTaxonomyFile(); 
				if (taxfile != "") { m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current taxonomy file and the taxonomy parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setTaxonomyFile(taxfile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            
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
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(taxfile); //if user entered a file with a path then preserve it	
			}
            
            string temp = validParameter.validFile(parameters, "relabund", false);		if (temp == "not found"){	temp = "false";			}
			relabund = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "printlevel", false);		if (temp == "not found"){	temp = "-1";		}
            m->mothurConvert(temp, printlevel);

            
            output = validParameter.validFile(parameters, "output", false);		if(output == "not found"){	output = "detail"; }
            if ((output != "simple") && (output != "detail")) { m->mothurOut(output + " is not a valid output form. Options are simple and detail. I will use detail."); m->mothurOutEndLine(); output = "detail"; }
			
            temp = validParameter.validFile(parameters, "threshold", false);			if (temp == "not found") { temp = "0"; }
            m->mothurConvert(temp, threshold);
            
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(taxfile);
                    parser.getNameFile(files);
                }
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "SummaryTaxCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SummaryTaxCommand::execute(){
	try{
        
        int maxLevel = findMaxLevel(taxfile);
        
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		int start = time(NULL);
		
        GroupMap* groupMap = NULL;
        CountTable* ct = NULL;
        if (groupfile != "") {
            groupMap = new GroupMap(groupfile);
            groupMap->readMap();
        }else if (countfile != "") {
            ct = new CountTable();
            ct->readTable(countfile, true, false);
        }
		
        PhyloSummary* taxaSum;
        if (countfile != "") { taxaSum = new PhyloSummary(ct, relabund, printlevel);
        }else { taxaSum = new PhyloSummary(groupMap, relabund, printlevel);  }
        
		if (m->control_pressed) { if (groupMap != NULL) { delete groupMap; } if (ct != NULL) { delete ct; } delete taxaSum; return 0; }
		
		int numSeqs = 0;
        map<string, vector<string> > nameMap;
        map<string, vector<string> >::iterator itNames;
        if (namefile != "") { m->readNames(namefile, nameMap); }
		
        ifstream in;
        m->openInputFile(taxfile, in);
        
        string name, taxon;
        while(!in.eof()){
            
            if (m->control_pressed) { break; }
            
            in >> name; m->gobble(in);
            taxon = m->getline(in); m->gobble(in);
            
            string newTax = m->addUnclassifieds(taxon, maxLevel, true);
            
            if (threshold != 0) {  newTax = processTaxMap(newTax);  }
            
            //add sequence to summary, countfile info included from Phylosummary constructor
            if (namefile != "") {
                itNames = nameMap.find(name);
                
                if (itNames == nameMap.end()) {
                    m->mothurOut(name + " is not in your name file please correct."); m->mothurOutEndLine(); exit(1);
                }else{
                    for (int i = 0; i < itNames->second.size(); i++) {
                        taxaSum->addSeqToTree(itNames->second[i], newTax);  //add it as many times as there are identical seqs
                    }
                    itNames->second.clear();
                    nameMap.erase(itNames->first);
                }
            }else {
                taxaSum->addSeqToTree(name, newTax);
            }
            
        }
        in.close();
        
		
		if (m->control_pressed) {  if (groupMap != NULL) { delete groupMap; } if (ct != NULL) { delete ct; } delete taxaSum; return 0; }
		
		//print summary file
		ofstream outTaxTree;
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(taxfile));
		string summaryFile = getOutputFileName("summary",variables);
		m->openOutputFile(summaryFile, outTaxTree);
		taxaSum->print(outTaxTree, output);
		outTaxTree.close();
		
		delete taxaSum;
        if (groupMap != NULL) { delete groupMap; } if (ct != NULL) { delete ct; }
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(summaryFile); m->mothurOutEndLine();	outputNames.push_back(summaryFile); outputTypes["summary"].push_back(summaryFile);
		m->mothurOutEndLine();
					
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
string SummaryTaxCommand::processTaxMap(string tax) {
    try{
        
        string newTax = tax;
        
        vector<string> taxons;
        int taxLength = tax.length();
        string taxon = "";
        int spot = 0;
        
        for(int i=0;i<taxLength;i++){
            
            
            if(tax[i] == ';'){
                
                int openParen = taxon.find_last_of('(');
                int closeParen = taxon.find_last_of(')');
                
                string newtaxon, confidence;
                if ((openParen != string::npos) && (closeParen != string::npos)) {
                    string confidenceScore = taxon.substr(openParen+1, (closeParen-(openParen+1)));
                    if (m->isNumeric1(confidenceScore)) {  //its a confidence
                        newtaxon = taxon.substr(0, openParen); //rip off confidence
                        confidence = taxon.substr((openParen+1), (closeParen-openParen-1));
                    }else { //its part of the taxon
                        newtaxon = taxon;
                        confidence = "0";
                    }
                }else{
                    newtaxon = taxon;
                    confidence = "-1";
                }
                float con = 0;
                
                m->mothurConvert(confidence, con);
                
                if (con == -1) { i += taxLength; } //not a confidence score, no confidence scores on this taxonomy
                else if ( con < threshold)  { spot = i; break; } //below threshold, set all to unclassified
                else {} //acceptable, move on
                taxons.push_back(taxon);
                
                taxon = "";
            }
            else{
                taxon += tax[i];
            }
            
        }
        
        if (spot != 0) {
            newTax = "";
            for (int i = 0; i < taxons.size(); i++) {  newTax += taxons[i] + ";";  }
            //for (int i = spot; i < taxLength; i++) {
                //if(tax[i] == ';'){   newTax += "unclassified;"; }
                m->removeConfidences(newTax);
            //}
        }else { m->removeConfidences(tax); newTax = tax; } //leave tax alone
        
        return newTax;
    }
    catch(exception& e) {
        m->errorOut(e, "SummaryTaxCommand", "processTaxMap");
        exit(1);
    }
}
/**************************************************************************************************/
int SummaryTaxCommand::findMaxLevel(string file) {
    try{
        GroupMap* groupMap = NULL;
        PhyloSummary taxaSum(groupMap, false, -1);
        
        taxaSum.summarize(file);
       
        return taxaSum.getMaxLevel();
    }
    catch(exception& e) {
        m->errorOut(e, "SummaryTaxCommand", "findMaxLevel");
        exit(1);
    }
}
/**************************************************************************************/


