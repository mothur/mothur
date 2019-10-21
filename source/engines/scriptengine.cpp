//
//  scriptengine.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "scriptengine.hpp"

/***********************************************************************/
ScriptEngine::ScriptEngine(string tpath, string commandString) : Engine(tpath){
    try {
        //remove quotes
        listOfCommands = commandString.substr(1, (commandString.length()-1));
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
        
        return nextcommand;
    }
    catch(exception& e) {
        m->errorOut(e, "ScriptEngine", "getNextCommand");
        exit(1);
    }
}
/***********************************************************************/
