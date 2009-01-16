#ifndef READDISTCOLUMNFILECOMMAND_H
#define READDISTCOLUMNFILECOMMAND_H
/*
 *  readdistcolumnfilecommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include "command.hpp"
#include "readmatrix.hpp"

/* The read.column command is used to read a distance matrix file in column format.  
The read.column command parameter options are distfile, namefile, cutoff and precision. 
The read.column command should be in the following format: read.column (distfile=yourDistFile, 
namefile=yourNameFile, cutoff=yourCutoff, precision=yourPrecision). The distfile and namefile parameters are required.  
If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.  */

class NameAssignment;
class GlobalData;


class ReadDistColumnFileCommand : public Command {
public:
	ReadDistColumnFileCommand();
	~ReadDistColumnFileCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	double cutoff;
	int precision;
	ReadMatrix* read;
	string filename, format, method;
	NameAssignment* nameMap;
};

#endif