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
        helpString += "The set.dir command parameter is seed, and it is required.\n";
        helpString += "To seed random set seed=yourRandomValue. By default mothur seeds random with the start time.\n";
          helpString += "Example set.seed(seed=12345).\n";
        helpString += "Note: No spaces between parameter labels (i.e. seed), '=' and parameters (i.e.yourSeedValue).\n";
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
        abort = false; calledHelp = false;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string, string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            //check to make sure all parameters are valid for command
            for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) {
                if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
            }
            
            bool seed = false;
            string temp = validParameter.validFile(parameters, "seed", false);
            if (temp == "not found") { random = 0;  m->mothurOut("[ERROR]: You must provide a seed value or set seed to clear."); m->mothurOutEndLine(); abort = true;}
            else if (temp == "clear") {
                random = time(NULL);
                seed = true;
            }else {
                if (m->isInteger(temp)) { m->mothurConvert(temp, random); seed = true; }
                else { m->mothurOut("[ERROR]: Seed must be an integer for the set.dir command."); m->mothurOutEndLine(); abort = true; }
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
        
        if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        srand(time(NULL));
        m->mothurOut("\nRandom number(1-100) with mothur's default of current time: " + toString(rand()%100 + 1) + ".\n");
        srand(random);
        m->mothurOut("Random number(1-100) seeded with " + toString(random) + ": " + toString(rand()%100 + 1) + ".\n\n");
       
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SetSeedCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************/
