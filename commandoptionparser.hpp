#ifndef COMMANDOPTIONPARSER_HPP
#define COMMANDOPTIONPARSER_HPP

#include "mothur.h"

//**********************************************************************************************************************

class CommandOptionParser {
public:
	CommandOptionParser(string);
	string getCommandString();
	
private:
	string commandString;
};

//**********************************************************************************************************************

#endif
