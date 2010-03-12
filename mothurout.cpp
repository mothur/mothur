/*
 *  m->mothurOut.cpp
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"

/******************************************************/
MothurOut* MothurOut::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new MothurOut();
	}
	return _uniqueInstance;
}
/*********************************************************************************************/
void MothurOut::setFileName(string filename)  {
	try {
		logFileName = filename;
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		openOutputFile(filename, out);
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFileName");
		exit(1);
	}
}
/*********************************************************************************************/
MothurOut::~MothurOut() {
	try {
		_uniqueInstance = 0;
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out.close();
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}

/*********************************************************************************************/
void MothurOut::mothurOut(string output) {
	try {
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		cout << output;
		out << output;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine() {
	try {
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		cout << endl;
		out << endl;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToLog(string output) {
	try {
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out << output;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutJustToLog");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::errorOut(exception& e, string object, string function) {
	mothurOut("Error: ");
	mothurOut(toString(e.what()));
	mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
	mothurOutEndLine();
}
/*********************************************************************************************/





