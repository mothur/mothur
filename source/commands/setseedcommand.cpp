//
//  setseedcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/24/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "setseedcommand.h"

//**********************************************************************************************************************
vector<string> SetSeedCommand::setParameters(){
    try {
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "SetSeedCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string SetSeedCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The set.seed command is used to seed random.\n";
        helpString += "The set.seed command parameter is seed, and it is required.\n";
        helpString += "To seed random set seed=yourRandomValue. By default mothur seeds random with the start time.\n";
          helpString += "Example set.seed(seed=12345).\n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "SetSeedCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************

SetSeedCommand::SetSeedCommand(string option)  {
    try {
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string, string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
    
            bool seed = false;
            string temp = validParameter.valid(parameters, "seed");
            if (temp == "not found") { random = 0;  m->mothurOut("[ERROR]: You must provide a seed value or set seed to clear.\n");  abort = true;}
            else if (temp == "clear") {
                random = time(NULL);
                seed = true;
            }else {
                if (util.isInteger(temp)) { util.mothurConvert(temp, random); seed = true; }
                else { m->mothurOut("[ERROR]: Seed must be an integer for the set.dir command.\n"); abort = true; }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SetSeedCommand", "SetSeedCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int SetSeedCommand::execute(){
    try {
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        m->setRandomSeed(random);
        m->mothurOut("Setting random seed to " + toString(random) + ".\n\n");
       
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SetSeedCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************/
