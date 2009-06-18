/*
 *  readdistcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readdistcommand.h"
#include "readphylip.h"
#include "readcolumn.h"
#include "readmatrix.hpp"

ReadDistCommand::ReadDistCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip", "column", "name", "cutoff", "precision", "group"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			globaldata->newRead();
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  globaldata->setPhylipFile(phylipfile);  globaldata->setFormat("phylip"); 	}
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  globaldata->setColumnFile(columnfile); globaldata->setFormat("column");	}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else {  
				globaldata->setGroupFile(groupfile); 
				//groupMap = new GroupMap(groupfile);
				//groupMap->readMap();
			}

			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else {  globaldata->setNameFile(namefile);	}

			
			//you are doing a list and group shared
			if ((phylipfile != "") && (groupfile != "")) { 
			globaldata->setFormat("matrix"); }
			
			if ((phylipfile == "") && (columnfile == "")) { cout << "When executing a read.dist command you must enter a phylip or a column." << endl; abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { cout << "When executing a read.dist command you must enter ONLY ONE of the following: phylip or column." << endl; abort = true; }
		
			if (columnfile != "") {
				if (namefile == "") {  cout << "You need to provide a namefile if you are going to use the column format." << endl; abort = true; }
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter.validFile(parameters, "precision", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			if (abort == false) {
				distFileName = globaldata->inputFileName;
				format = globaldata->getFormat();	
		
				if (format == "column") { read = new ReadColumnMatrix(distFileName); }	
				else if (format == "phylip") { read = new ReadPhylipMatrix(distFileName); }
				else if (format == "matrix") { 
					groupMap = new GroupMap(groupfile);
					groupMap->readMap();
					if (globaldata->gGroupmap != NULL) { delete globaldata->gGroupmap;  }
					globaldata->gGroupmap = groupMap;
				}
		
				if (format != "matrix" ) {
					read->setCutoff(cutoff);
	
					if(namefile != ""){	
						nameMap = new NameAssignment(namefile);
						nameMap->readMap(1,2);
					}else{
						nameMap = NULL;
					}
				}
			}

		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistCommand class Function ReadDistCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistCommand class function ReadDistCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************

void ReadDistCommand::help(){
	try {
		cout << "The read.dist command parameter options are phylip or column, group, name, cutoff and precision" << "\n";
		cout << "The read.dist command can be used in two ways.  The first is to read a phylip or column and run the cluster command" << "\n";
		cout << "For this use the read.dist command should be in the following format: " << "\n";
		cout << "read.dist(phylip=yourDistFile, name=yourNameFile, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The phylip or column parameter is required, but only one may be used.  If you use a column file the name filename is required. " << "\n";
		cout << "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed." << "\n";
		cout << "The second way to use the read.dist command is to read a phylip or column and a group, so you can use the libshuff command." << "\n";
		cout << "For this use the read.dist command should be in the following format: " << "\n";
		cout << "read.dist(phylip=yourPhylipfile, group=yourGroupFile). The cutoff and precision parameters are not valid with this use.  " << "\n";
		cout << "Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourPhylipfile)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

ReadDistCommand::~ReadDistCommand(){
	if (abort == false) {
		if (format != "matrix") { delete read; delete nameMap; }
	}
}

//**********************************************************************************************************************
int ReadDistCommand::execute(){
	try {
		
		if (abort == true) {	return 0;	}
		
		if (format == "matrix") {
			ifstream in;
			openInputFile(distFileName, in);
			matrix = new FullMatrix(in); //reads the matrix file
			in.close();
			//memory leak prevention
			if (globaldata->gMatrix != NULL) { delete globaldata->gMatrix;  }
			globaldata->gMatrix = matrix; //save matrix for coverage commands
		}else {
			read->read(nameMap);
			//to prevent memory leak

			if (globaldata->gListVector != NULL) {  delete globaldata->gListVector;  }
			globaldata->gListVector = read->getListVector();

			if (globaldata->gSparseMatrix != NULL) { delete globaldata->gSparseMatrix;  }
			globaldata->gSparseMatrix = read->getMatrix();

		}
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadDistCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadDistCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
