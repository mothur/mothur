#ifndef Mothur_newcommandtemplate_h
#define Mothur_newcommandtemplate_h

//
//  newcommandtemplate.h
//  Mothur
//
//  Created by westcott on 5/3/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//
//test

//*********Be sure to change ifdef and define to a unique name.**************//

/* This class is designed to provide a template for creating new commands. 
 It includes code snippets to make creating the command classes virtually pure
 functions easier. It includes sample parameter declaration and parameter checking,
 as well as reference to other classes you may find helpful.
 It also includes the code needed to read a sharedfile. It is a work in progress so 
 please add things you may find helpful to yourself or other developers trying to 
 add commands to mothur.
 
*/

#include "command.hpp"

/**************************************************************************************************/

class NewCommand : public Command {
public:
    NewCommand(string);
    ~NewCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "newCommandNameToBeSeenByUser";			}
    string getCommandCategory()		{ return "commandCategory";		} 
    
    string getOutputPattern(string);
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
	string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/newCommandNameToBeSeenByUser"; }
    string getDescription()		{ return "brief description"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort;
    vector<string> outputNames;
};

/**************************************************************************************************/




#endif
