//
//  makelefse.cpp
//  Mothur
//
//  Created by SarahsWork on 6/3/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "makelefsecommand.h"
#include "designmap.h"

//**********************************************************************************************************************
vector<string> MakeLefseCommand::setParameters(){
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none","lefse",false,false,true); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none","lefse",false,false,true); parameters.push_back(prelabund);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","",false,false,false); parameters.push_back(pconstaxonomy);
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,false, true); parameters.push_back(pdesign);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
        CommandParameter pscale("scale", "Multiple", "totalgroup-totalotu-averagegroup-averageotu", "totalgroup", "", "", "","",false,false); parameters.push_back(pscale);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["lefse"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeLefseCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The make.lefse command allows you to create a lefse formatted input file from mothur's output files.\n";
		helpString += "The make.lefse command parameters are: shared, relabund, constaxonomy, design, scale, groups and label.  The shared or relabund are required.\n";
		helpString += "The shared parameter is used to input your shared file, http://www.wiki.mothur.org/wiki/Shared_file.\n";
        helpString += "The relabund parameter is used to input your relabund file, http://www.wiki.mothur.org/wiki/Relabund_file.\n";
        helpString += "The design parameter is used to input your design file, http://www.wiki.mothur.org/wiki/Design_File.\n";
        helpString += "The constaxonomy parameter is used to input your taxonomy file. http://www.wiki.mothur.org/wiki/Constaxonomy_file. The contaxonomy file is the taxonomy file outputted by classify.otu(list=yourListfile, taxonomy=yourTaxonomyFile). Be SURE that the you are the constaxonomy file distance matches the shared file distance.  ie, for *.0.03.cons.taxonomy set label=0.03. Mothur is smart enough to handle shared files that have been subsampled. \n";
        helpString += "The scale parameter allows you to select what scale you would like to use to convert your shared file abundances to relative abundances. Choices are totalgroup, totalotu, averagegroup, averageotu, default is totalgroup.\n";
		helpString += "The label parameter allows you to select what distance level you would like used, if none is given the first distance is used.\n";
		helpString += "The make.lefse command should be in the following format: make.lefse(shared=yourSharedFile)\n";
		helpString += "make.lefse(shared=final.an.shared)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeLefseCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "lefse") {  pattern = "[filename],[distance],lefse"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeLefseCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeLefseCommand::MakeLefseCommand(string option)  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { abort = true; }
			else if (designfile == "not found") { designfile = ""; }
			else { current->setDesignFile(designfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else {  current->setRelAbundFile(relabundfile); }
            
            constaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
			if (constaxonomyfile == "not open") { constaxonomyfile = ""; abort = true; }
			else if (constaxonomyfile == "not found") {  constaxonomyfile = "";  }
			
			label = validParameter.valid(parameters, "label");
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile.\n");  label=""; }
            
            string groups = validParameter.valid(parameters, "groups");
			if (groups == "not found") { groups = "";  }
			else { util.splitAtDash(groups, Groups); if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } } }


            if ((relabundfile == "") && (sharedfile == "")) {
				//is there are current file available for either of these?
				//give priority to shared, then relabund
				sharedfile = current->getSharedFile();
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else {
					relabundfile = current->getRelAbundFile();
					if (relabundfile != "") {   m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter.\n");  }
					else {
						m->mothurOut("No valid current files. You must provide a shared or relabund.\n"); 
						abort = true;
					}
				}
			}
			
			if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("[ERROR]: You may not use both a shared and relabund file.\n");  abort = true;  }
			 
            scale = validParameter.valid(parameters, "scale");				if (scale == "not found") { scale = "totalgroup"; }
			
			if ((scale != "totalgroup") && (scale != "totalotu") && (scale != "averagegroup") && (scale != "averageotu")) {
				m->mothurOut(scale + " is not a valid scaling option for the get.relabund command. Choices are totalgroup, totalotu, averagegroup, averageotu.\n");  abort = true;
			}
        }
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "MakeLefseCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MakeLefseCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
      
        map<int, consTax2> consTax;
        if (constaxonomyfile != "") {  util.readConsTax(constaxonomyfile, consTax);  }
        
        if (m->getControl_pressed()) { return 0; }
        
        if (sharedfile != "") {
            inputFile = sharedfile;
            SharedRAbundFloatVectors* lookup = getSharedRelabund();
            runRelabund(consTax, lookup);
        }else {
            inputFile = relabundfile;
            SharedRAbundFloatVectors* lookup = getRelabund();
            runRelabund(consTax, lookup);
        }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeLefseCommand::runRelabund(map<int, consTax2>& consTax, SharedRAbundFloatVectors*& lookup){
	try {
        if (outputdir == "") { outputdir = util.hasPath(inputFile); }
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFile));
        variables["[distance]"] = lookup->getLabel();
		string outputFile = getOutputFileName("lefse",variables);
		outputNames.push_back(outputFile); outputTypes["lefse"].push_back(outputFile);
        
        ofstream out;
        util.openOutputFile(outputFile, out);

        DesignMap* designMap = NULL;
        vector<string> namesOfGroups = lookup->getNamesGroups();
        if (designfile != "") {
            designMap = new DesignMap(designfile);
            vector<string> categories = designMap->getNamesOfCategories();
            
            if (categories.size() > 3) {  m->mothurOut("\n[NOTE]: LEfSe input files allow for a class, subclass and subject.  More than 3 categories can cause formatting errors.\n\n"); }
            
            for (int j = 0; j < categories.size(); j++) {
                out << categories[j];
                for (int i = 0; i < namesOfGroups.size()-1; i++) {
                    if (m->getControl_pressed()) { out.close(); delete designMap; return 0; }
                    string value = designMap->get(namesOfGroups[i], categories[j]);
                    if (value == "not found") {
                        m->mothurOut("[ERROR]: " + namesOfGroups[i] + " is not in your design file, please correct.\n"); m->setControl_pressed(true);
                    }else { out  << '\t' << value; }
                }
                string value = designMap->get(namesOfGroups[namesOfGroups.size()-1], categories[j]);
                if (value == "not found") {
                    m->mothurOut("[ERROR]: " + namesOfGroups[namesOfGroups.size()-1] + " is not in your design file, please correct.\n"); m->setControl_pressed(true);
                }else { out << '\t' << value; }
                out << endl;
            }
        }
        
        out << "group";
        for (int i = 0; i < namesOfGroups.size(); i++) {  out  << '\t' << namesOfGroups[i]; }
        out << endl;
        
        for (int i = 0; i < lookup->getNumBins(); i++) { //process each otu
            if (m->getControl_pressed()) { break; }
            string nameOfOtu = lookup->getOTUName(i);
            
            if (constaxonomyfile != "") { //try to find the otuName in consTax to replace with consensus taxonomy
                int simpleLabel;
                util.mothurConvert(util.getSimpleLabel(nameOfOtu), simpleLabel);
                map<int, consTax2>::iterator it = consTax.find(simpleLabel);
                if (it != consTax.end()) {
                    nameOfOtu = it->second.taxonomy;
                    
                    //add sanity check abundances here??
                    string fixedName = "";
                    //remove confidences and change ; to |
                    util.removeConfidences(nameOfOtu);
                    for (int j = 0; j < nameOfOtu.length(); j++) {
                        if (nameOfOtu[j] == ';') { fixedName += '|'; }
                        else { fixedName += nameOfOtu[j]; }
                    }
                    nameOfOtu = fixedName + lookup->getOTUName(i) + "|";
                }else {
                    m->mothurOut("[ERROR]: can't find " + nameOfOtu + " in constaxonomy file. Do the distances match, did you forget to use the label parameter?\n"); m->setControl_pressed(true);
                }
                
            }
            //print name
            out << nameOfOtu;
            
            //print out relabunds for each otu
            vector<float> abunds = lookup->getOTU(i);
            for (int j = 0; j < abunds.size(); j++) { out  << '\t' << abunds[j]; }
            out << endl;
        }
        out.close();

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundFloatVectors* MakeLefseCommand::getSharedRelabund(){
	try {
		InputData input(sharedfile, "sharedfile", Groups);
		SharedRAbundVectors* templookup = input.getSharedRAbundVectors();
        Groups = templookup->getNamesGroups();
		string lastLabel = templookup->getLabel();
        
		if (label == "") { label = lastLabel;  }
		else {
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> labels; labels.insert(label);
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((templookup != NULL) && (userLabels.size() != 0)) {
                if (m->getControl_pressed()) {  delete templookup; return NULL;  }
                
                if(labels.count(templookup->getLabel()) == 1){
                    processedLabels.insert(templookup->getLabel());
                    userLabels.erase(templookup->getLabel());
                    break;
                }
                
                if ((util.anyLabelsToProcess(templookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = templookup->getLabel();
                    
                    delete templookup;
                    templookup = input.getSharedRAbundVectors(lastLabel);
                    
                    processedLabels.insert(templookup->getLabel());
                    userLabels.erase(templookup->getLabel());
                    
                    //restore real lastlabel to save below
                    templookup->setLabels(saveLabel);
                    break;
                }
                
                lastLabel = templookup->getLabel();
                
                //get next line to process
                //prevent memory leak
                delete templookup;
                templookup = input.getSharedRAbundVectors();
            }
            
            
            if (m->getControl_pressed()) { delete templookup; return NULL;  }
            
            //output error messages about any remaining user labels
            set<string>::iterator it;
            bool needToRun = false;
            for (it = userLabels.begin(); it != userLabels.end(); it++) {
                m->mothurOut("Your file does not include the label " + *it);
                if (processedLabels.count(lastLabel) != 1) {
                    m->mothurOut(". I will use " + lastLabel + ".\n"); 
                    needToRun = true;
                }else {
                    m->mothurOut(". Please refer to " + lastLabel + ".\n"); 
                }
            }
            
            //run last label if you need to
            if (needToRun )  {
                delete templookup;
                templookup = input.getSharedRAbundVectors(lastLabel);
            }
		}
        
        vector<SharedRAbundVector*> data = templookup->getSharedRAbundVectors();
        vector<string> otuNames = templookup->getOTUNames();
        delete templookup;
        SharedRAbundFloatVectors* lookup = new SharedRAbundFloatVectors();
        
        //convert to relabund
        for (int i = 0; i < data.size(); i++) {
            SharedRAbundFloatVector* rel = new SharedRAbundFloatVector();
            rel->setGroup(data[i]->getGroup());
            rel->setLabel(data[i]->getLabel());
            
			for (int j = 0; j < data[i]->getNumBins(); j++) {
                
				if (m->getControl_pressed()) { for (int k = 0; k < data.size(); k++) {  delete data[k];  } return lookup; }
                
				int abund = data[i]->get(j);
				float relabund = 0.0;
				
				if (scale == "totalgroup") {
					relabund = abund / (float) data[i]->getNumSeqs();
				}else if (scale == "totalotu") {
					//calc the total in this otu
					int totalOtu = 0;
					for (int l = 0; l < data.size(); l++) {  totalOtu += data[l]->get(j); }
					relabund = abund / (float) totalOtu;
				}else if (scale == "averagegroup") {
					relabund = abund / (float) (data[i]->getNumSeqs() / (float) data[i]->getNumBins());
				}else if (scale == "averageotu") {
					//calc the total in this otu
					int totalOtu = 0;
					for (int l = 0; l < data.size(); l++) {  totalOtu += data[l]->get(j); }
					float averageOtu = totalOtu / (float) data.size();
					
					relabund = abund / (float) averageOtu;
				}else{ m->mothurOut(scale + " is not a valid scaling option.\n");  m->setControl_pressed(true);  }
				
				rel->push_back(relabund);
			}
            
            lookup->push_back(rel);
            
        }
        for (int k = 0; k < data.size(); k++) {  delete data[k];  } data.clear();
        lookup->setOTUNames(otuNames);
        lookup->eliminateZeroOTUS();
        
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "getSharedRelabund");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundFloatVectors* MakeLefseCommand::getRelabund(){
	try {
		InputData input(relabundfile, "relabund", Groups);
		SharedRAbundFloatVectors* lookupFloat = input.getSharedRAbundFloatVectors();
		string lastLabel = lookupFloat->getLabel();
        Groups = lookupFloat->getNamesGroups();
		
		if (label == "") { label = lastLabel; return lookupFloat; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookupFloat != NULL) && (userLabels.size() != 0)) {
			
			if (m->getControl_pressed()) {  return lookupFloat;  }
			
			if(labels.count(lookupFloat->getLabel()) == 1){
				processedLabels.insert(lookupFloat->getLabel());
				userLabels.erase(lookupFloat->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(lookupFloat->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookupFloat->getLabel();
				
                delete lookupFloat;
				lookupFloat = input.getSharedRAbundFloatVectors(lastLabel);
				
				processedLabels.insert(lookupFloat->getLabel());
				userLabels.erase(lookupFloat->getLabel());
				
				//restore real lastlabel to save below
				lookupFloat->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookupFloat->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookupFloat;
			lookupFloat = input.getSharedRAbundFloatVectors();
		}
		
		
		if (m->getControl_pressed()) { return lookupFloat;  }
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it);
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + ".\n"); 
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + ".\n"); 
			}
		}
		
		//run last label if you need to
		if (needToRun )  {
			delete lookupFloat;
			lookupFloat = input.getSharedRAbundFloatVectors(lastLabel);
		}
		
		return lookupFloat;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLefseCommand", "getRelabund");
        exit(1);
	}
}

//**********************************************************************************************************************

