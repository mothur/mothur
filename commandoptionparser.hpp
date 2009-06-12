#ifndef COMMANDOPTIONPARSER_HPP
#define COMMANDOPTIONPARSER_HPP

#include "mothur.h"

//**********************************************************************************************************************

class CommandOptionParser {
public:
	CommandOptionParser(string);
	string getCommandString();
	string getOptionString();
	
private:
	string commandString, optionString;
};

//**********************************************************************************************************************

#endif
