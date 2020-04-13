//
//  scriptengine.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "scriptengine.hpp"

/***********************************************************************/
ScriptEngine::ScriptEngine(string tpath, string commandString, map<string, string> ev) : Engine(tpath){
    try {
        //remove quotes
        listOfCommands = commandString.substr(1, (commandString.length()-1));
        setEnvironmentVariables(ev);
    }
    catch(exception& e) {
        m->errorOut(e, "ScriptEngine", "ScriptEngine");
        exit(1);
    }
}
/***********************************************************************/

ScriptEngine::~ScriptEngine(){
    time_t end = time(NULL);
    m->mothurOut("\n\nIt took " + toString(end-start) + " seconds to run " + toString(numCommandsRun) + " commands from your script.\n\n");
}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on mothur
bool ScriptEngine::getInput(){
    try {
            
        string input = "";
        string commandName = "";
        string options = "";
        
        int quitCommandCalled = 0;
    
        while(quitCommandCalled != 1){
            
            input = getNextCommand(listOfCommands);
            
            if (input == "") { input = "quit()"; }
            
            m->appendLogBuffer("\nmothur > " + input + "\n");
        
            if (m->getControl_pressed()) { input = "quit()"; }
                
            //allow user to omit the () on the quit command
            if (input == "quit") { input = "quit()"; }
            if (input == "help") { input = "help()"; }

            CommandOptionParser parser(input);
            commandName = parser.getCommandString();
            options = parser.getOptionString();
                                        
            if (commandName != "") {
                numCommandsRun++;
                m->setExecuting(true);
                m->resetCommandErrors();
                
                //executes valid command
                m->setChangedSeqNames(true);
               
                Command* command = cFactory->getCommand(commandName, options);
                quitCommandCalled = command->execute();
                delete command;
                
                //if we aborted command
                if (quitCommandCalled == 2) {  m->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
                
                if (m->getControl_pressed()) { break;  }
                m->setControl_pressed(false);
                m->setExecuting(false);
                                
            }else {    m->mothurOut("[ERROR]: Invalid.\n"); }
        }
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "ScriptEngine", "getInput");
        exit(1);
    }
}
/***********************************************************************/
string ScriptEngine::getNextCommand(string& commandString) {
    try {
        
        string nextcommand = "";
        int count = 0;
        bool ignoreSemiColons = false;
        
        //go through string until you reach ; or end
        while (count < commandString.length()) {
            
             //you want to ignore any ; until you reach the next '
            if ((commandString[count] == '\'') && (!ignoreSemiColons)) {  ignoreSemiColons = true;  }
            else if ((commandString[count] == '\'') && (ignoreSemiColons)) {  ignoreSemiColons = false;  }
                
            if ((commandString[count] == ';') && (!ignoreSemiColons)) {  break;   }
            else {        nextcommand += commandString[count];    }
            
            count++;
        }
        
        //if you are not at the end
        if (count != commandString.length())  {   commandString = commandString.substr(count+1, commandString.length());  }
        else { commandString = ""; }
                
        
        //get rid of spaces in between commands if any
        if (commandString.length() > 0) {
            while (commandString[0] == ' ') {
                commandString = commandString.substr(1,commandString.length());
                if (commandString.length() == 0) {  break;  }
            }
        }
        
        string type = findType(nextcommand);
        
        if (type == "environment") {
            //set environmental variables
            string key, value; value = nextcommand;
            util.splitAtEquals(key, value);
            
            map<string, string>::iterator it = environmentalVariables.find(key);
            if (it == environmentalVariables.end())     { environmentalVariables[key] = value;  }
            else                                        { it->second = value;                   }
            
            m->mothurOut("Setting environment variable " + key + " to " + value + "\n");
            
            nextcommand = getNextCommand(commandString);
            
        }else { //assume command, look for environmental variables to replace
            
            int evPos = nextcommand.find_first_of('$');
            if (evPos == string::npos) { }//no '$' , nothing to do
            else { replaceVariables(nextcommand); }
        }
        
        double ramUsed, total;
        ramUsed = util.getRAMUsed(); total = util.getTotalRAM();
        m->mothurOut("RAM used: " + toString(ramUsed/(double)GIG) + "Gigabytes . Total Ram: " + toString(total/(double)GIG) + "Gigabytes.\n\n");
        
        return nextcommand;
    }
    catch(exception& e) {
        m->errorOut(e, "ScriptEngine", "getNextCommand");
        exit(1);
    }
}
/***********************************************************************/
