#ifndef GETCOMMANDINFOCOMMAND_H
#define GETCOMMANDINFOCOMMAND_H

/*
 *  getcommandinfo.h
 *  Mothur
 *
 *  Created by westcott on 4/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */



#include "command.hpp"
#include "commandfactory.hpp"

/**********************************************************/

class GetCommandInfoCommand : public Command {
	
public:
	GetCommandInfoCommand(string);
	GetCommandInfoCommand() { abort = true; calledHelp = true; setParameters(); }
	~GetCommandInfoCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.commandinfo";		}
	string getCommandCategory()		{ return "Hidden";				}
	string getHelpString();	
    string getOutputPattern(string) {  return "";  }	
	string getCitation() { return "no citation"; }
	string getDescription()		{ return "get.commandinfo"; }
	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	CommandFactory* commandFactory;
	string output;
	bool abort;
	vector<string> outputNames;
	
	int getInfo(vector<CommandParameter>, vector<string>&, vector<string>&, vector<string>&, vector<string>&, vector<string>&, map<string, string>&);

	
	
};

/**********************************************************/

#endif

