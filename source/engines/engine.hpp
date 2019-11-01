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
    
protected:
	vector<string> options;
	CommandFactory* cFactory;
	MothurOut* m;
    CurrentFile* current;
    Utils util;
    time_t start;
    int numCommandsRun;
    string path;
};


#endif
