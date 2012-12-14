#ifndef COMMAND_HPP
#define COMMAND_HPP
//test2
/*
 *  command.h
 *  
 *
 *  Created by Pat Schloss on 10/23/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

/*This class is a parent to all the command classes.  */


#include "mothur.h"
#include "optionparser.h"
#include "validparameter.h"
#include "mothurout.h"
#include "commandparameter.h"


class Command {
	
	public:
		Command() {  m = MothurOut::getInstance();   } 
		
		//needed by gui
		virtual string getCommandName() = 0;
		virtual string getCommandCategory() = 0;
		virtual string getHelpString() = 0;
		virtual string getCitation() = 0;
		virtual string getDescription() = 0;
		
		virtual map<string, vector<string> > getOutputFiles() { return outputTypes; }
        string getOutputFileName(string type, map<string, string> variableParts) {  //uses the pattern to create an output filename for a given type and input file name.  
               try {
                    string filename = "";
                    map<string, vector<string> >::iterator it;
                    
                    //is this a type this command creates
                    it = outputTypes.find(type);
                    if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
                    else {
                        
                        string patternTemp = getOutputPattern(type);
                        vector<string> patterns;
                        m->splitAtDash(patternTemp, patterns);
                        
                        //find pattern to use based on number of variables passed in
                        string pattern = "";
                        bool foundPattern = false;
                        for (int i = 0; i < patterns.size(); i++) {
                            int numVariables = 0;
                            for (int j = 0; j < patterns[i].length(); j++) { if (patterns[i][j] == '[') { numVariables++; } }
                            
                            if (numVariables == variableParts.size()) { pattern = patterns[i]; foundPattern = true; break; }
                        }
                        
                        if (!foundPattern) {  m->mothurOut("[ERROR]: Not enough variable pieces for " + type + ".\n"); m->control_pressed = true; }
                        
                        if (pattern != "") {
                            int numVariables = 0;
                            for (int i = 0; i < pattern.length(); i++) { if (pattern[i] == '[') { numVariables++; } }
                            
                            vector<string> pieces;
                            m->splitAtComma(pattern, pieces);
                            
                        
                            for (int i = 0; i < pieces.size(); i++) {
                                if (pieces[i][0] == '[') {
                                    map<string, string>::iterator it = variableParts.find(pieces[i]);
                                    if (it == variableParts.end()) {
                                        m->mothurOut("[ERROR]: Did not provide variable for " + pieces[i] + ".\n"); m->control_pressed = true;
                                    }else {
                                        if (it->second != "") {
                                            if (it->first == "[filename]") { filename += it->second; }
                                            else if (it->first == "[extension]") { 
                                                if (filename.length() > 0) { //rip off last "."
                                                    filename = filename.substr(0, filename.length()-1);
                                                }
                                                filename += it->second + "."; 
                                            }else { filename += it->second + "."; }
                                        }
                                    }
                                }else {
                                    filename += pieces[i] + ".";
                                }
                            }
                            if (filename.length() > 0) { //rip off last "."
                                filename = filename.substr(0, filename.length()-1);
                            }
                        }
                    }
                    return filename;
                }
                catch(exception& e) {
                    m->errorOut(e, "command", "getOutputFileName");
                    exit(1);
                }
        }
        
        virtual string getOutputPattern(string) = 0; //pass in type, returns something like: [filename],align or [filename],[distance],subsample.shared  strings in [] means its a variable.  This is used by the gui to predict output file names.  use variable keywords: [filename], [distance], [group], [extension], [tag]
		virtual vector<string> setParameters() = 0; //to fill parameters
		virtual vector<CommandParameter> getParameters() { return parameters; }
	
		virtual int execute() = 0;
		virtual void help() = 0;
		void citation() { m->mothurOutEndLine(); m->mothurOut(getCitation()); m->mothurOutEndLine(); }
		virtual ~Command() { }
	
	protected:
		MothurOut* m;
		bool calledHelp;
			
		map<string, vector<string> > outputTypes;
		vector<CommandParameter> parameters;
	
		map<string, vector<string> >::iterator itTypes;
};

#endif
