#ifndef LIBSHUFFCOMMAND_H
#define LIBSHUFFCOMMAND_H

/*
 *  libshuffcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"

using namespace std;

class GlobalData;

class LibShuffCommand : public Command {
	
	public:
		LibShuffCommand();	
		~LibShuffCommand();
		int execute();	
	
	private:
		GlobalData* globaldata;
		
};

#endif
