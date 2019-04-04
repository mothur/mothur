//
//  mergeotuscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 12/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "mergeotuscommand.hpp"


//**********************************************************************************************************************
vector<string> MergeOTUsCommand::setParameters(){
    try {
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","constaxonomy",false,true, true); parameters.push_back(pconstaxonomy);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","shared",false,true,true); parameters.push_back(pshared);
        CommandParameter prelabund("relabund", "InputTypes", "", "", "none", "none", "none","relabund",false,true,true); parameters.push_back(prelabund);
        CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","list",false,true,true); parameters.push_back(plist);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter ptaxlevel("taxlevel", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(ptaxlevel);
        
        
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string MergeOTUsCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The merge.otus command parameters are shared, list, relabund, constaxonomy, taxlevel and label.  constaxonomy is a required, unless you have a valid current file.\n";
        helpString += "The taxlevel parameter allows you to specify the taxonomy level you would like to use when merging. Default=maxlevel.\n";
        helpString += "Example merge.otus(shared=yourSharedFile, constaxonomy=yourConsTaxonomyFile).\n";
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string MergeOTUsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared")           {  pattern = "[filename],merge,[extension]"; }
        else if (type == "list")        {  pattern = "[filename],merge,[extension]"; }
        else if (type == "relabund")    {  pattern = "[filename],merge,[extension]"; }
        else if (type == "constaxonomy") {  pattern = "[filename],[label],merge,cons.taxonomy"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MergeOTUsCommand::MergeOTUsCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["relabund"] = tempOutNames;
        outputTypes["constaxonomy"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "MergeOTUsCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

MergeOTUsCommand::MergeOTUsCommand(string option)  {
    try {
        abort = false; calledHelp = false;
        allLines = true;
        
        //allow user to run help
        if(option == "help") {  help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters  = parser.getParameters();
            map<string,string>::iterator it;
            
            ValidParameters validParameter;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
            }
            
            vector<string> tempOutNames;
            outputTypes["shared"] = tempOutNames;
            outputTypes["list"] = tempOutNames;
            outputTypes["relabund"] = tempOutNames;
            outputTypes["constaxonomy"] = tempOutNames;
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){	inputDir = "";		}
            else {
                string path;
                it = parameters.find("shared");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["shared"] = inputDir + it->second;		}
                }
                
                it = parameters.find("constaxonomy");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
                }
                
                it = parameters.find("relabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("list");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["list"] = inputDir + it->second;		}
                }
            }
            
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not found") { sharedfile = ""; }
            else if (sharedfile == "not open") { sharedfile = ""; abort = true; }
            else { current->setSharedFile(sharedfile); }
            
            listfile = validParameter.validFile(parameters, "list");
            if (listfile == "not found") { listfile = ""; }
            else if (listfile == "not open") { listfile = ""; abort = true; }
            else { current->setListFile(listfile); }
            
            relabundfile = validParameter.validFile(parameters, "relabund");
            if (relabundfile == "not found") { relabundfile = ""; }
            else if (relabundfile == "not open") { relabundfile = ""; abort = true; }
            else { current->setRelAbundFile(relabundfile); }
            
            constaxfile = validParameter.validFile(parameters, "constaxonomy"); //required
            if (constaxfile == "not found") {
                constaxfile = current->getConsTaxonomyFile();
                if (constaxfile != "") { m->mothurOut("Using " + constaxfile + " as input file for the constaxonomy parameter.\n");  }
                else { 	m->mothurOut("[ERROR]: You have no current constaxonomy file and the constaxonomy parameter is required.\n"); abort = true; }
            }
            else if (constaxfile == "not open") { constaxfile = ""; abort = true; }
            else { current->setConsTaxonomyFile(constaxfile); }
            
            if ((relabundfile == "") && (listfile == "") && (sharedfile == "")) { //no files to merge provided, look for currents
                //is there are current file available for any of these?
                //give priority to shared, then list, then relabund
                //if there is a current shared file, use it
                sharedfile = current->getSharedFile();
                if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
                else {
                    listfile = current->getListFile();
                    if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
                    else {
                        relabundfile = current->getRelAbundFile();
                        if (relabundfile != "") { m->mothurOut("Using " + relabundfile + " as input file for the rabund parameter.\n");  }
                        else {
                            m->mothurOut("[ERROR]: No valid current files. You must provide a list, relabund or shared file.\n");  abort = true;
                        }
                    }
                }
            }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){
                outputDir = "";
                outputDir += util.hasPath(constaxfile); //if user entered a file with a path then preserve it
            }
            
            //check for optional parameter and set defaults
            // ...at some point should added some additional type checking...
            label = validParameter.valid(parameters, "label");
            if (label == "not found") { label = ""; }
            else {
                if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
                else { allLines = true;  }
            }
            
            string temp = validParameter.valid(parameters, "taxlevel");		if (temp == "not found")  { temp = "-1"; }
            util.mothurConvert(temp, taxLevelCutoff); //-1 means use max level
        }
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "MergeOTUsCommand");
        exit(1);
    }
}

//**********************************************************************************************************************

MergeOTUsCommand::~MergeOTUsCommand(){}

//**********************************************************************************************************************

int MergeOTUsCommand::execute(){
    try {
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        //read in consensus taxonomy
        vector<consTax> conTaxes = util.readConsTax(constaxfile);
        
        PhyloTree tree;
        //add consensus taxonomy to phylo tree
        for (size_t i = 0; i < conTaxes.size(); i++) {
            tree.addSeqToTree(conTaxes[i].name, conTaxes[i].taxonomy);
        }
        //build tree
        tree.assignHeirarchyIDs(0);
        
        //get max level of phylotree
        int maxlevel = tree.getMaxLevel();
        
        //is the taxlevel parameter valid - not greater than the max level in the file. If greater reduce to maxlevel.
        if (taxLevelCutoff == -1) { //default, use maxlevel
            taxLevelCutoff = maxlevel;
        }else if ( taxLevelCutoff > maxlevel) { //invalid taxlevel, use maxlevel
            m->mothurOut("[WARNING]: The taxlevel selected is larger than the maxlevel in your constaxonomy file, disregarding. Using the max level of " + toString(maxlevel) + " for the taxlevel parameter.\n");
            taxLevelCutoff = maxlevel;
        }
        
        for (size_t i = 0; i < conTaxes.size(); i++) {
            string otuTax = conTaxes[i].taxonomy;
            if (taxLevelCutoff != maxlevel) { otuTax = util.trimTax(otuTax, taxLevelCutoff); }
            otuLabel2ConsTax[conTaxes[i].name] = otuTax;
            if (listfile != "") {   otuLabel2ConsSize[conTaxes[i].name] = conTaxes[i].abundance; }
        }
        conTaxes.clear();
        
        //extract tree nodes at taxlevel
        vector<TaxNode> thisLevelsNodes = tree.getNodes(taxLevelCutoff);
        
        //merge otus at each node at this level
        if (listfile != "")             {    mergeListOTUs(thisLevelsNodes);       }
        else if (sharedfile != "")      {    mergeSharedOTUs(thisLevelsNodes);     }
        else if (relabundfile != "")    {    mergeRelabundOTUs(thisLevelsNodes);   }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); }  return 0; }
        
        string currentName = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }
        
        itTypes = outputTypes.find("shared");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
        }
        
        itTypes = outputTypes.find("relabund");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRelAbundFile(currentName); }
        }
        
        //set constaxonomy file as new current constaxonomyfile
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
        }
        
        m->mothurOut("\nOutput File Names:\n");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	} m->mothurOutEndLine();

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeOTUsCommand::mergeSharedOTUs(vector<TaxNode>& nodes){
    try {
        string numNodes = toString(nodes.size());
        
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[extension]"] = util.getExtension(sharedfile);
        string outputFileName = getOutputFileName("shared", variables);
        outputTypes["shared"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        InputData input(sharedfile, "sharedfile", Groups);
        SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
        string lastLabel = lookup->getLabel();
        Groups = lookup->getNamesGroups();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        bool printHeaders = true;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->getControl_pressed()) {  out.close(); delete lookup;  return 0; }
            
            if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
                
                m->mothurOut(lookup->getLabel()+"\t"+numNodes+"\n");
                process(lookup, out, printHeaders, nodes);
                
                processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
            }
            
            if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup->getLabel();
                
                delete lookup;
                lookup = input.getSharedRAbundVectors(lastLabel);
                m->mothurOut(lookup->getLabel()+"\t"+numNodes+"\n");
                
                process(lookup, out, printHeaders, nodes);
                
                processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
                
                //restore real lastlabel to save below
                lookup->setLabels(saveLabel);
            }
            
            lastLabel = lookup->getLabel();
            //prevent memory leak
            delete lookup;
            
            if (m->getControl_pressed()) {  out.close();  return 0; }
            
            //get next line to process
            lookup = input.getSharedRAbundVectors();
        }
        
        if (m->getControl_pressed()) { out.close(); return 0; }
        
        //output error messages about any remaining user labels
        bool needToRun = false;
        for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
            m->mothurOut("Your file does not include the label " + *it);
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
            else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
        }
        
        //run last label if you need to
        if (needToRun )  {
            delete lookup;
            lookup = input.getSharedRAbundVectors(lastLabel);
            
            m->mothurOut(lookup->getLabel()+"\t"+numNodes+"\n");
            process(lookup, out, printHeaders, nodes);
            
            delete lookup;
        }
        
        out.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "mergeSharedOTUs");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeOTUsCommand::process(SharedRAbundVectors*& thisLookUp, ofstream& out, bool& printHeaders, vector<TaxNode>& nodes){
    try {
        vector<string> groups = thisLookUp->getNamesGroups();
        //create SharedRAbundVectors for the merged groups. Fill with blank rabundFloatVectors
        SharedRAbundVectors* merged; merged = new SharedRAbundVectors();
        for (int i = 0; i < groups.size(); i++) {
            SharedRAbundVector* myLookup = new SharedRAbundVector();
            myLookup->setLabel(thisLookUp->getLabel());
            myLookup->setGroup(groups[i]);
            merged->push_back(myLookup);
        }
        
        //translate otuNames to bin numbers
        map<string, int> otuName2BinNumber;
        map<string, int>::iterator it;
        for (int j = 0; j < thisLookUp->getNumBins(); j++) {
            if (m->getControl_pressed()) { break; }
            
            otuName2BinNumber[thisLookUp->getOTUName(j)] = j;
        }
        
        if (m->getControl_pressed()) { delete merged; return 0; }
        
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(constaxfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxfile));
        variables["[label]"] = thisLookUp->getLabel();
        string outputFileName = getOutputFileName("constaxonomy", variables);
        outputTypes["constaxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        ofstream outCons;
        util.openOutputFile(outputFileName, outCons);
        
        outCons << "OTU\tSize\tTaxonomy\n";
    
        for (int i = 0; i < nodes.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            vector<string> otuNames = nodes[i].accessions; //otus to merge
            string newOTUName = otuNames[0];
            vector<int> abunds; abunds.resize(groups.size(), 0);
            
            for (int j = 0; j < otuNames.size(); j++) {
                it = otuName2BinNumber.find(otuNames[j]); //do we have this otu in the shared file
                
                if (it != otuName2BinNumber.end()) { //we found it
                    vector<int> thisOtusAbunds = thisLookUp->getOTU(it->second);
                    
                    for (int k = 0; k < thisOtusAbunds.size(); k++) {  abunds[k] += thisOtusAbunds[k];  } //add this otus abunds to merged otu abunds
                    
                }else { m->mothurOut("[ERROR]: missing otu " + otuNames[j] + " from shared file, cannot continue.\n"); m->setControl_pressed(true); break;  }
            }
            
            merged->push_back(abunds, newOTUName);
            
            //merge consensus taxonomy results
            int sumOtu = util.sum(abunds);
            
            outCons << newOTUName << '\t' << sumOtu << '\t' << otuLabel2ConsTax[newOTUName] << endl;
        }
        
        if (m->getControl_pressed()) { delete merged; return 0; }
        
        merged->eliminateZeroOTUS(); // remove any zero OTUs created by median option.
        
        //print new file
        merged->print(out, printHeaders);
        delete merged;
        
        outCons.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeOTUsCommand::mergeListOTUs(vector<TaxNode>& nodes){
    try {
        string numNodes = toString(nodes.size());
        
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        string outputFileName = getOutputFileName("list", variables);
        outputTypes["list"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        InputData input(listfile, "list", Groups);
        ListVector* list = input.getListVector();
        string lastLabel = list->getLabel();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        bool printHeaders = true;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->getControl_pressed()) {  out.close(); delete list;  return 0; }
            
            if(allLines == 1 || labels.count(list->getLabel()) == 1){
                
                m->mothurOut(list->getLabel()+"\t"+numNodes+"\n");
                process(list, out, printHeaders, nodes);
                
                processedLabels.insert(list->getLabel()); userLabels.erase(list->getLabel());
            }
            
            if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = list->getLabel();
                
                delete list;
                list = input.getListVector(lastLabel);
                m->mothurOut(list->getLabel()+"\t"+numNodes+"\n");
                
                process(list, out, printHeaders, nodes);
                
                processedLabels.insert(list->getLabel()); userLabels.erase(list->getLabel());
                
                //restore real lastlabel to save below
                list->setLabel(saveLabel);
            }
            
            lastLabel = list->getLabel();
            //prevent memory leak
            delete list;
            
            if (m->getControl_pressed()) {  out.close();  return 0; }
            
            //get next line to process
            list = input.getListVector();
        }
        
        if (m->getControl_pressed()) { out.close(); return 0; }
        
        //output error messages about any remaining user labels
        bool needToRun = false;
        for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
            m->mothurOut("Your file does not include the label " + *it);
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
            else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
        }
        
        //run last label if you need to
        if (needToRun )  {
            delete list;
            list = input.getListVector(lastLabel);
            
            m->mothurOut(list->getLabel()+"\t"+numNodes+"\n");
            process(list, out, printHeaders, nodes);
            
            delete list;
        }
        
        out.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "mergeListOTUs");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeOTUsCommand::process(ListVector*& list, ofstream& out, bool& printHeaders, vector<TaxNode>& nodes){
    try {
        
        ListVector* merged; merged = new ListVector();
        merged->setLabel(list->getLabel());
        
        //translate otuNames to bin numbers
        map<string, int> otuName2BinNumber;
        map<string, int>::iterator it;
        for (int j = 0; j < list->getNumBins(); j++) {
            if (m->getControl_pressed()) { break; }
            
            otuName2BinNumber[list->getOTUName(j)] = j;
        }
        
        if (m->getControl_pressed()) { delete merged; return 0; }
        
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(constaxfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxfile));
        variables["[label]"] = list->getLabel();
        string outputFileName = getOutputFileName("constaxonomy", variables);
        outputTypes["constaxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        ofstream outCons;
        util.openOutputFile(outputFileName, outCons);
        
        outCons << "OTU\tSize\tTaxonomy\n";
        
        for (int i = 0; i < nodes.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            vector<string> otuNames = nodes[i].accessions; //otus to merge
            string newOTUName = otuNames[0];
            string mergedNames = "";
            it = otuName2BinNumber.find(newOTUName); //do we have this otu in the shared file
            
            int sizeOtu = 0;
            if (it != otuName2BinNumber.end()) { //we found it
                mergedNames = list->get(it->second);
                sizeOtu += otuLabel2ConsSize[newOTUName];
            }else { m->mothurOut("[ERROR]: missing otu " + newOTUName + " from list file, cannot continue.\n"); m->setControl_pressed(true); }

            for (int j = 1; j < otuNames.size(); j++) {
                it = otuName2BinNumber.find(otuNames[j]); //do we have this otu in the shared file
                
                if (it != otuName2BinNumber.end()) { //we found it
                    string bin = list->get(it->second);
                    
                    mergedNames += "," + bin;
                    sizeOtu += otuLabel2ConsSize[otuNames[j]];
                    
                }else { m->mothurOut("[ERROR]: missing otu " + otuNames[j] + " from list file, cannot continue.\n"); m->setControl_pressed(true); break;  }
            }
            
            int sumOtu = util.getNumNames(mergedNames);
            merged->push_back(mergedNames, sumOtu, newOTUName);
            
            outCons << newOTUName << '\t' << sizeOtu << '\t' << otuLabel2ConsTax[newOTUName] << endl;
        }
        
        if (m->getControl_pressed()) { delete merged; return 0; }
        
        //print new file
        merged->print(out, printHeaders);
        delete merged;
        
        outCons.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeOTUsCommand::mergeRelabundOTUs(vector<TaxNode>& nodes){
    try {
        string numNodes = toString(nodes.size());
        
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(relabundfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(relabundfile));
        variables["[extension]"] = util.getExtension(relabundfile);
        string outputFileName = getOutputFileName("relabund", variables);
        outputTypes["relabund"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        InputData input(relabundfile, "relabund", Groups);
        SharedRAbundFloatVectors* lookup = input.getSharedRAbundFloatVectors();
        string lastLabel = lookup->getLabel();
        Groups = lookup->getNamesGroups();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        bool printHeaders = true;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->getControl_pressed()) {  out.close(); delete lookup;  return 0; }
            
            if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
                
                m->mothurOut(lookup->getLabel()+"\t"+numNodes+"\n");
                process(lookup, out, printHeaders, nodes);
                
                processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
            }
            
            if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup->getLabel();
                
                delete lookup;
                lookup = input.getSharedRAbundFloatVectors(lastLabel);
                m->mothurOut(lookup->getLabel()+"\t"+numNodes+"\n");
                
                process(lookup, out, printHeaders, nodes);
                
                processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
                
                //restore real lastlabel to save below
                lookup->setLabels(saveLabel);
            }
            
            lastLabel = lookup->getLabel();
            //prevent memory leak
            delete lookup;
            
            if (m->getControl_pressed()) {  out.close();  return 0; }
            
            //get next line to process
            lookup = input.getSharedRAbundFloatVectors();
        }
        
        if (m->getControl_pressed()) { out.close(); return 0; }
        
        //output error messages about any remaining user labels
        bool needToRun = false;
        for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
            m->mothurOut("Your file does not include the label " + *it);
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
            else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
        }
        
        //run last label if you need to
        if (needToRun )  {
            delete lookup;
            lookup = input.getSharedRAbundFloatVectors(lastLabel);
            
            m->mothurOut(lookup->getLabel()+"\t"+numNodes+"\n");
            process(lookup, out, printHeaders, nodes);
            
            delete lookup;
        }
        
        out.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "mergeRelabundOTUs");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeOTUsCommand::process(SharedRAbundFloatVectors*& thisLookUp, ofstream& out, bool& printHeaders, vector<TaxNode>& nodes){
    try {
        vector<string> groups = thisLookUp->getNamesGroups();
        //create SharedRAbundVectors for the merged groups. Fill with blank rabundFloatVectors
        SharedRAbundFloatVectors* merged; merged = new SharedRAbundFloatVectors();
        for (int i = 0; i < groups.size(); i++) {
            SharedRAbundFloatVector* myLookup = new SharedRAbundFloatVector();
            myLookup->setLabel(thisLookUp->getLabel());
            myLookup->setGroup(groups[i]);
            merged->push_back(myLookup);
        }
        
        //translate otuNames to bin numbers
        map<string, int> otuName2BinNumber;
        map<string, int>::iterator it;
        for (int j = 0; j < thisLookUp->getNumBins(); j++) {
            if (m->getControl_pressed()) { break; }
            
            otuName2BinNumber[thisLookUp->getOTUName(j)] = j;
        }
        
        if (m->getControl_pressed()) { delete merged; return 0; }
        
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(constaxfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxfile));
        variables["[label]"] = thisLookUp->getLabel();
        string outputFileName = getOutputFileName("constaxonomy", variables);
        outputTypes["constaxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        ofstream outCons;
        util.openOutputFile(outputFileName, outCons);
        
        outCons << "OTU\tSize\tTaxonomy\n";
        
        for (int i = 0; i < nodes.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            vector<string> otuNames = nodes[i].accessions; //otus to merge
            string newOTUName = otuNames[0];
            vector<float> abunds; abunds.resize(groups.size(), 0);
            
            for (int j = 0; j < otuNames.size(); j++) {
                it = otuName2BinNumber.find(otuNames[j]); //do we have this otu in the shared file
                
                if (it != otuName2BinNumber.end()) { //we found it
                    vector<float> thisOtusAbunds = thisLookUp->getOTU(it->second);
                    
                    for (int k = 0; k < thisOtusAbunds.size(); k++) {  abunds[k] += thisOtusAbunds[k];  } //add this otus abunds to merged otu abunds
                    
                }else { m->mothurOut("[ERROR]: missing otu " + otuNames[j] + " from relabund file, cannot continue.\n"); m->setControl_pressed(true); break;  }
            }
            
            merged->push_back(abunds, newOTUName);
            
            //merge consensus taxonomy results
            float sumOtu = util.sum(abunds);
            
            outCons << newOTUName << '\t' << sumOtu << '\t' << otuLabel2ConsTax[newOTUName] << endl;
        }
        
        if (m->getControl_pressed()) { delete merged; return 0; }
        
        merged->eliminateZeroOTUS(); // remove any zero OTUs created by median option.
        
        //print new file
        merged->print(out, printHeaders);
        delete merged;
        
        outCons.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************


