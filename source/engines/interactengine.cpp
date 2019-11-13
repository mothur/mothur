//
//  interactengine.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "interactengine.hpp"

/***********************************************************************/

InteractEngine::InteractEngine(string tpath) : Engine(tpath) {
    
    if (m->getLogFileName() == "") {
        time_t ltime = time(NULL); /* calendar time */
        string outputPath = current->getOutputDir();
        //if (outputPath == "") { outputPath = current->getDefaultPath();  }
        string logFileName = outputPath + "mothur." + toString(ltime) + ".logfile";
        m->setLogFileName(logFileName, false);
        m->mothurOut("\n");
    }
}

/***********************************************************************/

InteractEngine::~InteractEngine(){}

/***********************************************************************/
//This function allows the user to input commands one line at a time until they quit.
//If the command is garbage it does nothing.
bool InteractEngine::getInput(){
    try {
        string input = "";
        string commandName = "";
        string options = "";
        int quitCommandCalled = 0;
        
        while(quitCommandCalled != 1){
            
            input = getCommand();
            
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
                
                m->setControl_pressed(false);
                m->setExecuting(false);
                
            }else { m->mothurOut("[ERROR]: Invalid.\n"); }
        }
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "InteractEngine", "getInput");
        exit(1);
    }
}
/***********************************************************************/
string InteractEngine::getCommand()  {
    try {
    
        #if defined NON_WINDOWS
            #ifdef USE_READLINE
                char* nextCommand = NULL;
                nextCommand = readline("\nmothur > ");
                
                if(nextCommand != NULL) {  add_history(nextCommand);  }
                else{ //^D causes null string and we want it to quit mothur
                    nextCommand = strdup("quit");
                }
                
                m->mothurOutJustToLog("\nmothur > " + toString(nextCommand) + "\n");
                string returnCommand = nextCommand;
                free(nextCommand);
                return returnCommand;
            #else
                string nextCommand = "";
                m->mothurOut("\nmothur > ");
                getline(cin, nextCommand);
                m->mothurOut("\n");
                m->mothurOutJustToLog("\nmothur > " + toString(nextCommand) + "\n");
                
                return nextCommand;
            #endif
        #else
                string nextCommand = "";
                
                m->mothurOut("\nmothur > ");
                getline(cin, nextCommand);
                m->mothurOut("\n");
                m->mothurOutJustToLog(toString(nextCommand) + "\n");
                
                return nextCommand;
        #endif
    
                        
    }
    catch(exception& e) {
        m->errorOut(e, "InteractEngine", "getCommand");
        exit(1);
    }
}
/***********************************************************************/
