#ifndef COMMANDOPTIONPARSER_HPP
#define COMMANDOPTIONPARSER_HPP

#include "mothur.h"
#include "mothurout.h"

//**********************************************************************************************************************

class CommandOptionParser {
public:
	CommandOptionParser(string);
	string getCommandString();
	string getOptionString();
	
private:
	string commandString, optionString;
	MothurOut* m;
};

//**********************************************************************************************************************

#endif
