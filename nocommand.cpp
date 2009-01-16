/*
 *  nocommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "nocommand.h"

//**********************************************************************************************************************

NoCommand::NoCommand(){}

//**********************************************************************************************************************

NoCommand::~NoCommand(){}

//**********************************************************************************************************************

int NoCommand::execute(){
	//Could choose to give more help here?fdsah
	cout << "Invalid command - choices are read.phylip(), read.column(), read.list(), read.rabund(), read.sabund(), cluster(), collect(), collect.shared(), parselist(), shared(), rarefaction(), summary(), summary.shared(), quit(), help()." << "\n";
	cout << "For more information on command parameters use the help() command." << "\n";
	return 0;
}

//**********************************************************************************************************************
