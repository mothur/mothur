//
//  makeclrcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/20/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "makeclrcommand.hpp"

//**********************************************************************************************************************
vector<string> MakeCLRCommand::setParameters(){
    try {
        CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none","clr",false,false,true); parameters.push_back(pshared);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
        CommandParameter pzero("zero", "Number", "", "0.1", "", "", "","",false,false); parameters.push_back(pzero);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["clr"] = tempOutNames;
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {    myArray.push_back(parameters[i].name);        }
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeCLRCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string MakeCLRCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The make.clr command parameters are shared, groups, zero and label. The shared file is required, unless you have a valid current file.\n";
        helpString += "The groups parameter allows you to specify which of the groups in your sharedfile you would like included. The group names are separated by dashes.\n";
        helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
        helpString += "The zero parameter allows you to set an value for zero OTUs. Default is 0.1.\n";
        helpString += "The make.clr command should be in the following format: make.clr(shared=yourSharedFile).\n";
        helpString += "Example make.clr(shared=final.opti_mcc.shared, zero=0.25).\n";
        helpString += "The default value for groups is all the groups in your sharedfile, and all labels in your inputfile will be used.\n";
        helpString += "The make.clr command outputs a .clr file.\n";
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeCLRCommand", "MakeCLRCommand");
        exit(1);
    }
}
//**********************************************************************************************************************
string MakeCLRCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "clr") {  pattern = "[filename],[distance],clr-[filename],clr"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeCLRCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

MakeCLRCommand::MakeCLRCommand(string option) {
    try {

        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") { sharedfile = ""; abort = true; }
            else if (sharedfile == "not found") {
                sharedfile = current->getSharedFile();
                if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n"); }
                else { m->mothurOut("[ERROR]: No valid current shared file. You must provide a shared file, quitting.\n"); abort = true; }
            }else {  current->setSharedFile(sharedfile); }
            
            if (outputdir == ""){    outputdir = util.hasPath(sharedfile);        }
            
            //check for optional parameter and set defaults
            // ...at some point should added some additional type checking...
            label = validParameter.valid(parameters, "label");
            if (label == "not found") { label = ""; }
            else {
                if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
                else { allLines = true;  }
            }
            
            groups = validParameter.valid(parameters, "groups");
            if (groups == "not found") { groups = "";  }
            else {
                util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
            }
            
            string temp = validParameter.valid(parameters, "zero"); if (temp == "not found") { temp = "0.1";  }
            util.mothurConvert(temp, zeroReplacementValue);
        }

    }
    catch(exception& e) {
        m->errorOut(e, "MakeCLRCommand", "MakeCLRCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int MakeCLRCommand::execute(){
    try {
    
        if (abort) { if (calledHelp) { return 0; }  return 2;    }
        
        InputData input(sharedfile, "sharedfile", Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
        string outputFileName = getOutputFileName("clr",variables);
        bool printHeaders = true;
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName); outputTypes["clr"].push_back(outputFileName);
        
        while (lookup != NULL) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            process(lookup, out, printHeaders); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
        out.close();
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {    util.mothurRemove(outputNames[i]);    } outputTypes.clear(); return 0;}
        
        string currentName = "";
        itTypes = outputTypes.find("clr");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCLRFile(currentName);  }
        }
        
        m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i] +"\n");     } m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeCLRCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
// compute geometric mean through formula
// antilog(((log(1) + log(2) + . . . + log(n))/n)


//x <- c(10, 5, 3, 1,0)
//> x[x==0] <- 0.1
//> log2(x / prod(x)^(1/4))
//[1]  2.3452054  1.3452054  0.6082399 -0.9767226 -4.2986507
void MakeCLRCommand::process(SharedRAbundVectors*& thisLookUp, ofstream& out, bool& printHeaders){
    try {
        vector<string> lookupGroups = thisLookUp->getNamesGroups();
        
        vector<SharedRAbundFloatVector*> lookup = thisLookUp->getSharedRAbundFloatVectors();
        vector<string> otuNames = thisLookUp->getOTUNames();
        
        if (printHeaders) {  out << "label\tGroup\tnumOtus\t" << util.getStringFromVector(otuNames, "\t") << endl; printHeaders = false; }
        
        for (int i = 0; i < lookup.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            vector<float> abunds = lookup[i]->get();
            
            double geoMean = util.geometricMean(abunds, zeroReplacementValue);
            
            for (int j = 0; j < abunds.size(); j++) { lookup[i]->set(j, log2(abunds[j]/geoMean)); }
            
            lookup[i]->print(out);
            
            delete lookup[i];
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeCLRCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************

