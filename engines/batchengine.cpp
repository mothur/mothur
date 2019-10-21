//
//  batchengine.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "batchengine.hpp"


/***********************************************************************/
//This function opens the batchfile to be used by BatchEngine::getInput.
BatchEngine::BatchEngine(string tpath, string batchFile) : Engine(tpath) {
    try {
        openedBatch = util.openInputFile(batchFile, inputBatchFile, "no error");
        if (!openedBatch) {
            if (util.checkLocations(batchFile, current->getLocations())) { openedBatch = util.openInputFile(batchFile, inputBatchFile); }
            else {  m->mothurOut("[ERROR]: unable to open " + batchFile + " batch file, please correct.\n");  }
        }
        
        batchFileName = batchFile;
                
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "BatchEngine");
        exit(1);
    }
}

/***********************************************************************/

BatchEngine::~BatchEngine(){
    time_t end = time(NULL);
    m->mothurOut("\n\nIt took " + toString(end-start) + " seconds to run " + toString(numCommandsRun) + " commands from " + batchFileName  + " batch file.\n\n");
}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on Dotur
bool BatchEngine::getInput(){
    try {
        //check if this is a valid batchfile
        if (!openedBatch) { return 1;  }
    
        string input = "";
        string commandName = "";
        string options = "";
        
        int quitCommandCalled = 0;
        int count = 0;
        while(quitCommandCalled != 1){
            
            input = getNextCommand(inputBatchFile);
            count++;
            
            if (input[0] != '#') {
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
                                        
                }else {     m->mothurOut("[ERROR]: Invalid.\n"); }
                
            }
            util.gobble(inputBatchFile);
        }
        
        inputBatchFile.close();
        return 1;
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "getInput");
        exit(1);
    }
}
/***********************************************************************/
string BatchEngine::getNextCommand(ifstream& inputBatchFile) {
    try {
            
        string nextcommand = "";
        
        if (inputBatchFile.eof()) { nextcommand = "quit()"; }
        else { nextcommand = util.getline(inputBatchFile); }
        
        return nextcommand;
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "getNextCommand");
        exit(1);
    }
}

/***********************************************************************/
