//
//  classifysharedcommand.cpp
//  Mothur
//
//  Created by Abu Zaher Md. Faridee on 8/13/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "classifysharedcommand.h"

  // TODO: finish implementation
ClassifySharedCommand::ClassifySharedCommand() {
}

  // TODO: finish implementation
ClassifySharedCommand::ClassifySharedCommand(string option) {
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
  return NULL;
}

string ClassifySharedCommand::getHelpString() {
  return NULL;
}

  // TODO: finish implementation
int ClassifySharedCommand::execute() {
  return 0;
}

