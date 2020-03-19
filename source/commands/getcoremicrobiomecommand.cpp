//
//  GetCoreMicroBiomeCommand.cpp
//  Mothur
//
//  Created by John Westcott on 5/8/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "getcoremicrobiomecommand.h"
#include "getrelabundcommand.h"


//**********************************************************************************************************************
vector<string> GetCoreMicroBiomeCommand::setParameters(){	
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none","coremicrobiom",false,false, true); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none","coremicrobiom",false,false, true); parameters.push_back(prelabund);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter poutput("output", "Multiple", "fraction-count", "fraction", "", "", "","",false,false); parameters.push_back(poutput);
        CommandParameter pabund("abundance", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pabund);
		CommandParameter psamples("samples", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(psamples);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCoreMicroBiomeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetCoreMicroBiomeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.coremicrobiome determines the fraction of OTUs that are found in varying numbers of samples for different minimum relative abundances.\n";
		helpString += "The get.coremicrobiome parameters are: shared, relabund, groups, label, output, abundance and samples. Shared or relabund is required.\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like analyzed.\n";
        helpString += "The output parameter is used to specify whether you would like the fraction of OTU's or OTU count outputted. Options are fraction or count. Default=fraction.\n";
		helpString += "The abundance parameter allows you to specify an abundance you would like the OTU names outputted for. Values 1 to 100, will be treated as the percentage.  For example relabund=0.01 can be set with abundance=1 or abundance=0.01.  For abundance values < 1 percent, abundance=0.001 will specify OTUs with relative abundance of 0.001.\n";
        helpString += "The samples parameter allows you to specify the minimum number of samples you would like the OTU names outputted for. Must be an interger between 1 and number of samples in your file.\n";
		helpString += "The new command should be in the following format: get.coremicrobiome(shared=yourSharedFile)\n";
		helpString += "get.coremicrobiom(shared=final.an.shared, abund=30)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCoreMicroBiomeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetCoreMicroBiomeCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "coremicrobiome") {  pattern = "[filename],[tag],core.microbiome"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetCoreMicroBiomeCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetCoreMicroBiomeCommand::GetCoreMicroBiomeCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["coremicrobiome"] = tempOutNames; 
	}
	catch(exception& e) {
		m->errorOut(e, "GetCoreMicroBiomeCommand", "GetCoreMicroBiomeCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetCoreMicroBiomeCommand::GetCoreMicroBiomeCommand(string option)  {
	try {
		abort = false; calledHelp = false;   allLines = true;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
            outputTypes["coremicrobiome"] = tempOutNames; 

			//check for parameters
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; format = "sharedfile"; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; format = "relabund"; current->setRelAbundFile(relabundfile); }
            
            if ((relabundfile == "") && (sharedfile == "")) { 
				//is there are current file available for either of these?
				//give priority to shared, then relabund
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") {  inputFileName = sharedfile; format="sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") {  inputFileName = relabundfile; format="relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a shared or relabund."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = util.hasPath(inputFileName); //if user entered a file with a path then preserve it	
			}
            
            string groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { util.splitAtDash(groups, Groups); if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } } }
            
            string label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
            
            output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "fraction"; }
						
			if ((output != "fraction") && (output != "count")) { m->mothurOut(output + " is not a valid output form. Options are fraction and count. I will use fraction."); m->mothurOutEndLine(); output = "fraction"; }
            
            string temp = validParameter.valid(parameters, "abundance");	if (temp == "not found"){	temp = "-1";	}
			util.mothurConvert(temp, abund);
            
            if (!util.isEqual(abund, -1)) {
                if ((abund < 0) || (abund > 100)) { m->mothurOut(toString(abund) + " is not a valid number for abund. Must be between 0 and 100.\n"); }
                if (abund < 1) { //convert
                    string temp = toString(abund); string factorString = "1";
                    bool found = false;
                    for (int i = 0; i < temp.length(); i++) {
                        if (temp[i] == '.') { found = true; }
                        else {
                            if (found) { factorString += "0"; }
                        }
                    }
                    util.mothurConvert(factorString, factor);
                }else {
                    factor = 100;
                    abund /= 100;
                }
            }else {
                factor = 100;
            }
            
            temp = validParameter.valid(parameters, "samples");	if (temp == "not found"){	temp = "-1";	}
			util.mothurConvert(temp, samples);

		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetCoreMicroBiomeCommand", "GetCoreMicroBiomeCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetCoreMicroBiomeCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (format == "sharedfile") { //convert to relabund
            string options = "shared=" + sharedfile;
            if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: get.relabund(" + options + ")\n");
            Command* relabundCommand = new GetRelAbundCommand(options);
            
            relabundCommand->execute();
            map<string, vector<string> > filenames = relabundCommand->getOutputFiles();
            relabundfile = filenames["relabund"][0];
            inputFileName = relabundfile; format="relabund";
            
            delete relabundCommand;
            current->setMothurCalling(false);
            m->mothurOut("/******************************************/\n");
        }
        
        InputData input(inputFileName, format, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundFloatVectors* lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        
        if (samples != -1) {
            if ((samples < 1) || (samples > lookup->size())) { m->mothurOut(toString(samples) + " is not a valid number for samples. Must be an integer between 1 and the number of samples in your file. Your file contains " + toString(lookup->size()) + " samples, so I will use that.\n"); samples = lookup->size(); }
        }
            
        while (lookup != NULL) {
                
            if (m->getControl_pressed()) { delete lookup; break; }
                
            createTable(lookup); delete lookup;
                
            lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
        }

        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); }  return 0; }
        
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "GetCoreMicroBiomeCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetCoreMicroBiomeCommand::createTable(SharedRAbundFloatVectors*& lookup){
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[tag]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("coremicrobiome", variables);
        outputNames.push_back(outputFileName);  outputTypes["coremicrobiome"].push_back(outputFileName);
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
        int numSamples = lookup->size();
        int numOtus = lookup->getNumBins();
        
        //table is 100 by numsamples
        //question we are answering is: what fraction of OTUs in a study have a relative abundance at or above %X
        //in at least %Y samples. x goes from 0 to 100, y from 1 to numSamples
        vector< vector<double> > table; table.resize(factor+1);
        for (int i = 0; i < table.size(); i++) { table[i].resize(numSamples, 0.0); }
        
        map<int, vector<string> > otuNames;
        if (!(util.isEqual(abund, -1)) && (samples == -1)) { //fill with all samples
            for (int i = 0; i < numSamples; i++) {
                vector<string> temp;
                otuNames[i+1] = temp;
            }
        }else if ((util.isEqual(abund, -1)) && (samples != -1)) { //fill with all relabund
            for (int i = 0; i < factor+1; i++) {
                vector<string> temp;
                otuNames[i] = temp;
            }
        }else if (!(util.isEqual(abund, -1)) && (samples != -1)) { //only one line is wanted
            vector<string> temp;
            int thisAbund = abund*factor;
            otuNames[thisAbund] = temp;
        }
        
        vector<string> sampleNames = lookup->getNamesGroups();
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < numOtus; i++) {
            
            if (m->getControl_pressed()) { break; }
            
            //count number of samples in this otu with a relabund >= spot in count
            vector<int> counts; counts.resize(factor+1, 0);
            
            for (int j = 0; j < sampleNames.size(); j++) {
                double relabund = lookup->get(i, sampleNames[j]);
                int wholeRelabund = (int) (floor(relabund*factor));
                for (int k = 0; k < wholeRelabund+1; k++) { counts[k]++; }
            }
            
            //add this otus info to table
            for (int j = 0; j < table.size(); j++) {
                for (int k = 0; k < counts[j]; k++) { table[j][k]++; }
                
                if ((util.isEqual(abund, -1)) && (samples != -1)) { //we want all OTUs with this number of samples
                    if (counts[j] >= samples) { otuNames[j].push_back(currentLabels[i]); }
                }else if (!(util.isEqual(abund, -1)) && (samples == -1)) { //we want all OTUs with this relabund
                    if (j == (int)(abund*factor)) {
                        for (int k = 0; k < counts[j]; k++) {  otuNames[k+1].push_back(currentLabels[i]); }
                    }
                }else if (!(util.isEqual(abund, -1)) && (samples != -1)) { //we want only OTUs with this relabund for this number of samples
                    if ((j == (int)(abund*factor)) && (counts[j] >= samples)) {
                        otuNames[j].push_back(currentLabels[i]);
                    }
                }
            }
        }
		
        //format output
		if (output == "fraction") { out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint); }
        out << "NumSamples\t";
        
        //convert table counts to percents
        int precisionLength = (toString(factor)).length();
        for (int i = 0; i < table.size(); i++) {
            out << "Relabund-" << setprecision(precisionLength-1)<< (float)(i/(float)factor) << "\t";
            if (m->getControl_pressed()) { break; }
            for (int j = 0; j < table[i].size(); j++) {  if (output == "fraction") { table[i][j] /= (double) numOtus; } }
        }
        out << endl;
        
        for (int i = 0; i < numSamples; i++) {
            if (m->getControl_pressed()) { break; }
            out << i+1;
            for (int j = 0; j < table.size(); j++) {  out << setprecision(6) << '\t' << table[j][i]; }
            out << endl;
        }

        out.close();
        
        if (m->getControl_pressed()) { return 0; }
        
        if ((samples != -1) || (!util.isEqual(abund, -1)))  {
            string outputFileName2 = outputDir + util.getRootName(util.getSimpleName(inputFileName)) + lookup->getLabel() + ".core.microbiomelist";
            outputNames.push_back(outputFileName2);  outputTypes["coremicrobiome"].push_back(outputFileName2);
            ofstream out2;
            util.openOutputFile(outputFileName2, out2);
            
            if ((util.isEqual(abund, -1)) && (samples != -1)) { //we want all OTUs with this number of samples
                out2 << "Relabund\tOTUList_for_samples=" << samples << "\n";
            }else if (!(util.isEqual(abund, -1)) && (samples == -1)) { //we want all OTUs with this relabund
                out2 << "Samples\tOTUList_for_abund=" << abund*factor << "\n";
            }else if (!(util.isEqual(abund, -1)) && (samples != -1)) { //we want only OTUs with this relabund for this number of samples
                out2 << "Relabund\tOTUList_for_samples=" << samples << "\n";
            }

            for (map<int, vector<string> >::iterator it = otuNames.begin(); it != otuNames.end(); it++) {
                if (m->getControl_pressed()) { break; }
                
                vector<string> temp = it->second;
                string list = util.makeList(temp);
                
                if (!(util.isEqual(abund, -1)) && (samples == -1)) { //fill with all samples
                    out2 << it->first << '\t' << list << endl;
                }else  { //fill with relabund
                    out2 << fixed << showpoint << setprecision(precisionLength-1) << (it->first/(float)(factor)) << '\t' << list << endl;
                }
            }
            
            out2.close();
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetCoreMicroBiomeCommand", "createTable");
		exit(1);
	}
}

//**********************************************************************************************************************


