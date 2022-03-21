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
BatchEngine::BatchEngine(string tpath, string batchFile, map<string, string> ev) : Engine(tpath) {
    try {
        batchFile = util.removeQuotes(batchFile);
        ifstream inBatchTest;
        openedBatch = util.openInputFile(batchFile, inBatchTest, "no error");
        if (!openedBatch) {
            if (util.checkLocations(batchFile, current->getLocations())) { openedBatch = util.openInputFile(batchFile, inBatchTest); }
            else {  m->mothurOut("[ERROR]: unable to open " + batchFile + " batch file, please correct.\n");  }
        }
        
        batchFileName = batchFile;
        noBufferNeeded = true;
        
        if (openedBatch) { //check for set.logfile
            string nextcommand = "#"; //force grabbing first command

            while (!inBatchTest.eof()) {
                
                nextcommand = util.getline(inBatchTest); util.gobble(inBatchTest);
                
                if (nextcommand[0] != '#') { //skip comments
                    
                    int pos = nextcommand.find("set.logfile");
                    if (pos != string::npos) { noBufferNeeded = false; break; }
                }
            }
            inBatchTest.close();
            
            openedBatch = util.openInputFile(batchFileName, inputBatchFile, "no error");
        }
        
        if (noBufferNeeded) {
            if (m->getLogFileName() == "") {
                time_t ltime = time(nullptr); /* calendar time */
                string outputPath = current->getOutputDir();
                string logFileName = outputPath + "mothur." + toString(ltime) + ".logfile";
                m->setLogFileName(logFileName, false);
                m->mothurOut("\n");
            }
        }
        
        setEnvironmentVariables(ev); //inherit environmental variables from nested batch files

        bstart = time(nullptr);
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
    time_t end = time(nullptr);
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

            CommandOptionParser parser(input);
            string commandName = parser.getCommandString();
            string options = parser.getOptionString();
            
            m->mothurOut("\nmothur > " + input + "\n");
                
            if (m->getControl_pressed()) { input = "quit()"; }
                                        
            if (commandName != "") {
                numCommandsRun++;
                m->setExecuting(true); m->resetCommandErrors(); m->setChangedSeqNames(true); m->setChangedGroupNames(true);
                            
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
        
        string type = findType(nextcommand);
       
        if (type == "batch") {
            m->mothurOut("/*****************************************************************************/\n");
            
            BatchEngine newBatchEngine(path, nextcommand, environmentalVariables);
            
            if (newBatchEngine.getOpenedBatch()) {
                bool bail = false;
                while(!bail)    {    bail = newBatchEngine.getInput();    }
                numBatches++;
            }
            m->mothurOut("/*****************************************************************************/\n");
            
            nextcommand = getNextCommand(inputBatchFile);
        }else if (type == "environment") {
            //set environmental variables
            string key, value; value = nextcommand;
            util.splitAtEquals(key, value);
            
            map<string, string>::iterator it = environmentalVariables.find(key);
            if (it == environmentalVariables.end())     { environmentalVariables[key] = value;  }
            else                                        { it->second = value;                   }
            
            m->mothurOut("Setting environment variable " + key + " to " + value + "\n");
            
            nextcommand = getNextCommand(inputBatchFile);
            
        }else { //assume command, look for environmental variables to replace
            
            int evPos = nextcommand.find_first_of('$');
            if (evPos == string::npos) { //no '$' , check for mothurhome
                evPos = nextcommand.find("mothurhome");
                if (evPos != string::npos) { replaceVariables(nextcommand); }
            }else { replaceVariables(nextcommand); }
        }
       
        if (m->getDebug()) {
            double ramUsed, total;
            ramUsed = util.getRAMUsed(); total = util.getTotalRAM();
            m->mothurOut("RAM used: " + toString(ramUsed/(double)GIG) + " Gigabytes. Total Ram: " + toString(total/(double)GIG) + " Gigabytes.\n\n");
        }
        
        return nextcommand;
    }
    catch(exception& e) {
        m->errorOut(e, "BatchEngine", "getNextCommand");
        exit(1);
    }
}
/***********************************************************************/

