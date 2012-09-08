//
//  classifysharedcommand.cpp
//  Mothur
//
//  Created by Abu Zaher Md. Faridee on 8/13/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "classifysharedcommand.h"
#include "fs/randomforest.hpp"
#include "fs/decisiontree.hpp"
#include "fs/treenode.hpp"


  // TODO: finish implementation
ClassifySharedCommand::ClassifySharedCommand() {
  try {
    abort = true; calledHelp = true;
    setParameters();
    vector<string> tempOutNames;
      // TODO update the outputTypes variable
      // START of segmenet to be updated
    outputTypes["fileType1"] = tempOutNames; //filetypes should be things like: shared, fasta, accnos...
    outputTypes["fileType2"] = tempOutNames;
    outputTypes["FileType3"] = tempOutNames;
      // END of segment to be updated
  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "ClassifySharedCommand");
    exit(1);
  }
}

  // TODO: finish implementation
  // WIP
ClassifySharedCommand::ClassifySharedCommand(string option) {
  try {
      ////////////////////////////////////////////////////////
      /////////////////// start leave alone block ////////////
      ////////////////////////////////////////////////////////
    abort = false; calledHelp = false;
    
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
        if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
      }
      
      
        //if the user changes the input directory command factory will send this info to us in the output parameter
      string inputDir = validParameter.validFile(parameters, "inputdir", false);
      if (inputDir == "not found"){	inputDir = "";		}
      else {
        
          ///////////////////////////////////////////////////////////////
          //////////////// stop leave alone block ///////////////////////
          ///////////////////////////////////////////////////////////////
        
          //edit file types below to include only the types you added as parameters
        
        string path;
        
        it = parameters.find("shared");
          //user has given a shared file
        if(it != parameters.end()){
          path = m->hasPath(it->second);
            //if the user has not given a path then, add inputdir. else leave path alone.
          if (path == "") {	parameters["shared"] = inputDir + it->second;		}
        }
        
        it = parameters.find("design");
          //user has given a design file
        if(it != parameters.end()){
          path = m->hasPath(it->second);
            //if the user has not given a path then, add inputdir. else leave path alone.
          if (path == "") {	parameters["design"] = inputDir + it->second;		}
        }
        
      }
        ///////////////////////////////////////////////////////////////////////////////
        /////////// example of getting filenames and checking dependancies ////////////
        // the validParameter class will make sure file exists, fill with correct    //
        // and name is current is given ///////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
      
      
        ///variables for examples below that you will most likely want to put in the header for
        //use by the other class functions.
      string phylipfile, columnfile, namefile, fastafile;
      string method;
      string splitcriteria;
      string otupersplit;
      vector<string> Estimators, Groups;
        //if allLines is used it should be initialized to 1 above.
      
      
        //check for parameters
        //get shared file, it is required
      sharedfile = validParameter.validFile(parameters, "shared", true);
      if (sharedfile == "not open") { sharedfile = ""; abort = true; }
      else if (sharedfile == "not found") {
          //if there is a current shared file, use it
        sharedfile = m->getSharedFile();
        if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
        else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
      }else { m->setSharedFile(sharedfile); }
      
        //get design file, it is required
      designfile = validParameter.validFile(parameters, "design", true);
      if (designfile == "not open") { sharedfile = ""; abort = true; }
      else if (designfile == "not found") {
          //if there is a current shared file, use it
        designfile = m->getDesignFile();
        if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
        else { 	m->mothurOut("You have no current designfile and the design parameter is required."); m->mothurOutEndLine(); abort = true; }
      }else { m->setDesignFile(designfile); }

      
        //if the user changes the output directory command factory will send this info to us in the output parameter
      outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){
        outputDir = m->hasPath(sharedfile); //if user entered a file with a path then preserve it
      }
      
      
        //////////////////////////////////////////////////////////////////////
        ////////// example of getting other types of parameters //////////////
        //////////////////////////////////////////////////////////////////////
      
        //use only one Mutliple type
//      method = validParameter.validFile(parameters, "method", false);
//      if (method == "not found") { method = "average"; }
//      if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { }
//      else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, and weighted."); m->mothurOutEndLine(); abort = true; }
//      
//        //use more than one multiple type. do not check to make sure the entry is valid.
//      string calc = validParameter.validFile(parameters, "calc", false);
//      if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
//      else {
//        if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
//      }
//      m->splitAtDash(calc, Estimators);
      
        // NEW CODE for OTU per split selection criteria
      otupersplit = validParameter.validFile(parameters, "otupersplit", false);
      if (otupersplit == "not found") { otupersplit = "log2"; }
      if ((otupersplit == "squareroot") || (otupersplit == "log2")) {
        optimumFeatureSubsetSelectionCriteria = otupersplit;
      }
      else { m->mothurOut("Not a valid OTU per split selection method. Valid OTU per split selection methods are 'log2' and 'squareroot'."); m->mothurOutEndLine(); abort = true; }
      
        // splitcriteria
      splitcriteria = validParameter.validFile(parameters, "splitcriteria", false);
      if (splitcriteria == "not found") { splitcriteria = "gainratio"; }
      if ((splitcriteria == "gainratio") || (splitcriteria == "infogain")) {
        treeSplitCriterion = splitcriteria;
      }
      else { m->mothurOut("Not a valid tree splitting criterio. Valid tree splitting criteria are 'gainratio' and 'infogain'."); m->mothurOutEndLine(); abort = true; }
      
        // boolean example
//        //Boolean type - m->isTrue looks for t, true, f or false and is case insensitive
//      string timing = validParameter.validFile(parameters, "timing", false);
//      if (timing == "not found") { timing = "F"; }
//      useTiming = m->isTrue(timing);

        // number example
//        //Number type - mothurConvert makes sure the convert can happen to avoid a crash.
//      string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
//      m->setProcessors(temp);
//      m->mothurConvert(temp, processors);
      
      string temp = validParameter.validFile(parameters, "numtrees", false); if (temp == "not found"){	temp = "100";	}
      m->mothurConvert(temp, numDecisionTrees);

        //Groups must be checked later to make sure they are valid. SharedUtilities has functions of check the validity, just make to so m->setGroups() after the checks.  If you are using these with a shared file no need to check the SharedRAbundVector class will call SharedUtilites for you, kinda nice, huh?
      string groups = validParameter.validFile(parameters, "groups", false);
      if (groups == "not found") { groups = ""; }
      else { m->splitAtDash(groups, Groups); }
      m->setGroups(Groups);
      
        //Commonly used to process list, rabund, sabund, shared and relabund files.  Look at "smart distancing" examples below in the execute function.
      string label = validParameter.validFile(parameters, "label", false);
      if (label == "not found") { label = ""; }
      else {
        if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
        else { allLines = 1;  }
      }
    }
    
  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "ClassifySharedCommand");
    exit(1);
  }
}

  // TODO: finish implementation
  // TODO: add more parameters like standardDeviationThreshold that was introduced in
  // python version
  // WIP 13 Aug
vector<string> ClassifySharedCommand::setParameters() {
  
    // define parameters here
    // we'll input a shared file and a design file name
  CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
  CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none", false, true); parameters.push_back(pdesign);
  
    // user will specify number of trees
  CommandParameter pnumtrees("numtrees", "Number", "", "100", "", "", "", false, false); parameters.push_back(pnumtrees);
  
    // user will specify tree splitting criteria
  CommandParameter psplitcriteria("splitcriteria", "Multiple", "gainratio-infogain", "infogain", "", "", "", false, false); parameters.push_back(psplitcriteria);
  
    // user will specify how much OTU to consider for each split of the total number of OTU
  CommandParameter potupersplitcriteria("otupersplit", "Multiple", "squareroot-log2", "log2", "", "", "", false, false); parameters.push_back(potupersplitcriteria);
  
    //  groups parameter is used to select the samples one would like to include from the shared file. The group parameter is used to provide a group file
  CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
  
    // user can specify to run the algo on only the labels specified from the shared file
  CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
  
    // set input and output folder
  CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
  CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);

  try {
    vector<string> myArray;
    for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
    return myArray;
  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "setParameters");
    exit(1);
  }

}

  // TODO: finish implementation
string ClassifySharedCommand::getOutputFileNameTag(string type, string inputName="") {
  try {
    string tag = "";
    map<string, vector<string> >::iterator it;
    
      //is this a type this command creates
    it = outputTypes.find(type);
    if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
    else {
        // TODO: update this code content, keeping this as placeholder
        // START of segment to be udpated
      if (type == "fileType1") {  tag = "tag1"; }
      else if (type == "fileType2") {  tag = "tag2"; }
      else if (type == "fileType3") {  tag = "tag3"; }
        // END of segment to be updated
      else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
    }
    return tag;
  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "getOutputFileName");
    exit(1);
  }

}

  // TODO: udpate the contents of the helpString
string ClassifySharedCommand::getHelpString() {
  try {
    string helpString = "";
    helpString += "dummy help string\n";
    return helpString;
  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "getHelpString");
    exit(1);
  }
}

  // TODO: finish implementation
int ClassifySharedCommand::execute() {
  try {
    
    if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
    
      // reading and processing a shared file code example
      // Note: As long as you set groups and labels as shown in the constructor, you can use this code without modification other than adding your function call which is passed the lookup vector.
      // The classes used below will handle the checking of groups to make sure they are valid and returning only the groups you selected.  The while loop implements mothur "smart distancing" so as long as you filled label as shown above in the constructor the code below will handle bad labels or labels not included in the sharedfile.
    
      //Reads sharefile, binLabels are stored in m->currentBinLabels, lookup will be filled with groups in m->getGroups() or all groups in file if m->getGroups is empty. If groups are selected, some bins maybe eliminated if they only contained seqs from groups not included. No need to worry about the details of this, SharedRAbundVector takes care of it.  Just make sure to use m->currentBinLabels if you are outputting OTU labels so that if otus are eliminated you still have the correct names.
    
      //    /*
    InputData input(sharedfile, "sharedfile");
    vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
        
      //read design file
    designMap = new GroupMap(designfile);
    designMap->readDesignMap();
    
    processSharedAndDesignData(lookup);
    
    string lastLabel = lookup[0]->getLabel();
    
      //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
    set<string> processedLabels;
    set<string> userLabels = labels;
    
    
      //as long as you are not at the end of the file or done wih the lines you want
    while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
      
      if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
      
      if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){
        
        m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
        
          ///////////////////////////////////////////////////////////////////////////////////
          //// Call your function to process specific distance in sharedfile, ie lookup /////
          ///////////////////////////////////////////////////////////////////////////////////
        
        processedLabels.insert(lookup[0]->getLabel());
        userLabels.erase(lookup[0]->getLabel());
      }
      
      if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
        string saveLabel = lookup[0]->getLabel();
        
        for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
        lookup = input.getSharedRAbundVectors(lastLabel);
        m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
        
          ///////////////////////////////////////////////////////////////////////////////////
          //// Call your function to process specific distance in sharedfile, ie lookup /////
          ///////////////////////////////////////////////////////////////////////////////////
        
        
        processedLabels.insert(lookup[0]->getLabel());
        userLabels.erase(lookup[0]->getLabel());
        
          //restore real lastlabel to save below
        lookup[0]->setLabel(saveLabel);
      }
      
      lastLabel = lookup[0]->getLabel();
        //prevent memory leak
      for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
      
      if (m->control_pressed) { return 0; }
      
        //get next line to process
      lookup = input.getSharedRAbundVectors();
    }
    
    if (m->control_pressed) {  return 0; }
    
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
    
      //run last label if you need to
    if (needToRun == true)  {
      for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }
      lookup = input.getSharedRAbundVectors(lastLabel);
      
      m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
      
        ///////////////////////////////////////////////////////////////////////////////////
        //// Call your function to process specific distance in sharedfile, ie lookup /////
        ///////////////////////////////////////////////////////////////////////////////////
      
      
      for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
      
    }
      //     */
        
      // TODO: re-enable this portion later and write to output files
      //output files created by command
//    m->mothurOutEndLine();
//    m->mothurOut("Output File Names: "); m->mothurOutEndLine();
//    for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
//    m->mothurOutEndLine();
    return 0;
    
  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "execute");
    exit(1);
  }
}

void ClassifySharedCommand::processSharedAndDesignData(vector<SharedRAbundVector*> lookup){  
//    for (int i = 0; i < designMap->getNamesOfGroups().size(); i++) {
//      string groupName = designMap->getNamesOfGroups()[i];
//      cout << groupName << endl;
//    }

//    for (int i = 0; i < designMap->getNumSeqs(); i++) {
//      string sharedGroupName = designMap->getNamesSeqs()[i];
//      string treatmentName = designMap->getGroup(sharedGroupName);
//      cout << sharedGroupName << " : " << treatmentName <<  endl;
//    }
  
  map<string, int> treatmentToIntMap;
  map<int, string> intToTreatmentMap;
  for (int  i = 0; i < designMap->getNumGroups(); i++) {
    string treatmentName = designMap->getNamesOfGroups()[i];
    treatmentToIntMap[treatmentName] = i;
    intToTreatmentMap[i] = treatmentName;
  }

  int numSamples = lookup.size();
  int numFeatures = lookup[0]->size();

  int numRows = numSamples;
  int numColumns = numFeatures + 1;           // extra one space needed for the treatment/outcome
  
  vector< vector<int> > dataSet(numRows, vector<int>(numColumns, 0));
  
  for (int i = 0; i < lookup.size(); i++) {
    vector<individual> data = lookup[i]->getData();
    string sharedGroupName = data[0].group;
    string treatmentName = designMap->getGroup(sharedGroupName);
    
    int j = 0;
    for (; j < data.size(); j++) {
      int otuCount = data[j].abundance;
      dataSet[i][j] = otuCount;
    }
    dataSet[i][j] = treatmentToIntMap[treatmentName];
  }
  
//  cout << dataSet << endl;
  
  RandomForest randomForest(dataSet, numDecisionTrees, treeSplitCriterion);
  randomForest.setMothurOut(m);
  
  randomForest.populateDecisionTrees();
  randomForest.calcForrestErrorRate();
  randomForest.calcForrestVariableImportance();
  
  m->mothurOutEndLine();

}

