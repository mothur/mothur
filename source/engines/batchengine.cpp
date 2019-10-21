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
        bstart = time(NULL);
        numBatches = 0;
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "BatchEngine");
        exit(1);
    }
}

/***********************************************************************/

BatchEngine::~BatchEngine(){
    string batchesOutput = "";
    if (numBatches != 0) {
        batchesOutput = " and " + toString(numBatches) + " batch file";
        if (numBatches > 1) { batchesOutput += "s"; }
    }
    time_t end = time(NULL);
    m->mothurOut("\n\nIt took " + toString(end-bstart) + " seconds to run " + toString(numCommandsRun) + " commands" + batchesOutput + " from " + batchFileName  + " batch file.\n\n");
}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on Dotur
bool BatchEngine::getInput(){
    try {
        //check if this is a valid batchfile
        if (!openedBatch) { return true;  }
    
        int quitCommandCalled = 0;
        while(quitCommandCalled != 1){
            
            string input = getNextCommand(inputBatchFile);
            
            m->appendLogBuffer("\nmothur > " + input + "\n");
                
            if (m->getControl_pressed()) { input = "quit()"; }

            CommandOptionParser parser(input);
            string commandName = parser.getCommandString();
            string options = parser.getOptionString();
                                        
            if (commandName != "") {
                numCommandsRun++;
                m->setExecuting(true);
                m->resetCommandErrors();
                m->setChangedSeqNames(true);
                            
                Command* command = cFactory->getCommand(commandName, options);
                quitCommandCalled = command->execute();
                delete command;
                            
                //if we aborted command
                if (quitCommandCalled == 2) {  m->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
                    
                if (m->getControl_pressed()) { break;  }
                m->setControl_pressed(false); m->setExecuting(false);
                                        
            }else {     m->mothurOut("[ERROR]: Invalid command.\n"); }
        }
        
        inputBatchFile.close();
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "getInput");
        exit(1);
    }
}
/***********************************************************************/
string BatchEngine::getNextCommand(ifstream& inputBatchFile) {
    try {
        string nextcommand = "#"; //force grabbing first command
        
        while (nextcommand[0] == '#') { //skip comments
            if (!inputBatchFile.eof()) {
                
                nextcommand = util.getline(inputBatchFile);
                util.gobble(inputBatchFile);
                
            }else { nextcommand = "quit()"; break; } //end of file, quit
        }
        
        //allow user to omit the () on the help and quit commands
        if (nextcommand == "quit") { nextcommand = "quit()"; }
        if (nextcommand == "help") { nextcommand = "help()"; }
        
        //determine if this is a command or batch file
        //we know commands must include '(' characters for search for that
        int openParen = nextcommand.find_first_of('(');
        if (openParen == string::npos) { //no '(' character -> assume not a command, treat as new batchfile
            m->mothurOut("/******************************************/\n");
            m->mothurOut("/******************************************/\n");
            
            BatchEngine newBatchEngine(path, nextcommand);
            bool bail = false;
            while(!bail)    {    bail = newBatchEngine.getInput();    }
            numBatches++;
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("/******************************************/\n");
            nextcommand = getNextCommand(inputBatchFile);
        }
        
        return nextcommand;
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "getNextCommand");
        exit(1);
    }
}

/***********************************************************************/
