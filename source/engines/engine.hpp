#ifndef ENGINE_HPP
#define ENGINE_HPP

/*
 *  engine.hpp
 *  
 *
 *  Created by Pat Schloss on 8/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */
 
#include "commandoptionparser.hpp"
#include "command.hpp"
#include "commandfactory.hpp"
#include "mothurout.h"

class Engine {
public:
    Engine(string tpath) {
        try {
            cFactory = CommandFactory::getInstance();
            m = MothurOut::getInstance();
            current = CurrentFile::getInstance();
            
            string temppath = tpath.substr(0, (tpath.find_last_of("othur")-5));
            
            //this will happen if you set the path variable to contain mothur's exe location
            if (temppath == "") { path = util.findProgramPath("mothur"); }
            else { path = temppath; }
            
            current->setProgramPath(util.getFullPathName(path));
            current->setBlastPath(current->getProgramPath());
            
            //if you haven't set your own location
            #ifdef MOTHUR_FILES
            #else
                    //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
                    if (current->getProgramPath() != "") { current->setDefaultPath(current->getProgramPath()); }
            #endif
            
            //if you haven't set your own location
            #ifdef MOTHUR_TOOLS
                current->setBlastPath(current->getToolsPath());
            #endif
            
            start = time(NULL);
            numCommandsRun = 0;
        }
        catch(exception& e) {
            m->errorOut(e, "Engine", "Engine");
            exit(1);
        }
    }
	virtual ~Engine(){}
	virtual bool getInput() = 0;
	virtual string getLogFileName()	{	return m->getLogFileName();  }

	vector<string> getOptions()		{	return options;		}
   
    virtual void replaceVariables(string& nextCommand) {
        for (map<string, string>::iterator it = environmentalVariables.begin(); it != environmentalVariables.end(); it++) {
            int pos = nextCommand.find("$"+it->first);
            while (pos != string::npos) { //allow for multiple uses of a environmental variable in a single command
                nextCommand.replace(pos,it->first.length()+1,it->second); //-1 to grab $char
                pos = nextCommand.find(it->first);
            }
        }
    }
    
    virtual string findType(string nextCommand) {
        string type = "command";
           
        //determine if this is a command or environmental variable
        //we know commands must include '(' characters for search for that
        int openParen = nextCommand.find_first_of('(');
        if (openParen == string::npos) { //no '(' character -> assume not a command, treat as environmental variable if contains an equals
            //if no '=' sign than not an environmental variable
            int equalsSign = nextCommand.find_first_of('=');
            if (equalsSign != string::npos) { //no '=' character -> assume not a environmental variable, treat as new batch
                type = "environment";
            }
        }
            
        return type;
    }
    
    virtual void setEnvironmentVariables(map<string, string> ev) {
        environmentalVariables = ev;
        
        //set HOME path is present in environment variables
        string homeEnvironmentTag = "HOMEPATH";
        string homeEnvironmentValue = "";
        #if defined NON_WINDOWS
             homeEnvironmentTag = "HOME";
        #endif
        
        map<string, string>::iterator it = environmentalVariables.find(homeEnvironmentTag);
        if (it != environmentalVariables.end()) { homeEnvironmentValue = it->second; }
        
        //parse PATH to set search locations for mothur tools
        //set HOME path is present in environment variables
        string pathEnvironmentTag = "PATH";
        string pathEnvironmentValue = "";
        char delim = ';';
        #if defined NON_WINDOWS
             delim = ':';
        #endif
        
        it = environmentalVariables.find(pathEnvironmentTag);
        if (it != environmentalVariables.end()) { pathEnvironmentValue = it->second; }

        vector<string> pathDirs;
        util.splitAtChar(pathEnvironmentValue, pathDirs, delim);

        if (m->getDebug()) {
            m->mothurOut("[DEBUG]: dir's in path:\n");
            for (int i = 0; i < pathDirs.size(); i++) { m->mothurOut("[DEBUG]: " + pathDirs[i] + "\n"); }
        }
        
        current->setPaths(pathDirs);
        current->setHomePath(homeEnvironmentValue);
    }
    
protected:
	vector<string> options;
	CommandFactory* cFactory;
	MothurOut* m;
    CurrentFile* current;
    Utils util;
    time_t start;
    int numCommandsRun;
    string path;
    map<string, string> environmentalVariables;
};


#endif
