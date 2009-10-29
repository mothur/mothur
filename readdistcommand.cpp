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
			
			if ((phylipfile == "") && (columnfile == "")) { mothurOut("When executing a read.dist command you must enter a phylip or a column."); mothurOutEndLine(); abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { mothurOut("When executing a read.dist command you must enter ONLY ONE of the following: phylip or column."); mothurOutEndLine(); abort = true; }
		
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
						nameMap->readMap();
					}else{
						nameMap = NULL;
					}
				}
			}

		}

	}
	catch(exception& e) {
		errorOut(e, "ReadDistCommand", "ReadDistCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ReadDistCommand::help(){
	try {
		mothurOut("The read.dist command parameter options are phylip or column, group, name, cutoff and precision\n");
		mothurOut("The read.dist command can be used in two ways.  The first is to read a phylip or column and run the cluster command\n");
		mothurOut("For this use the read.dist command should be in the following format: \n");
		mothurOut("read.dist(phylip=yourDistFile, name=yourNameFile, cutoff=yourCutoff, precision=yourPrecision) \n");
		mothurOut("The phylip or column parameter is required, but only one may be used.  If you use a column file the name filename is required. \n");
		mothurOut("If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.\n");
		mothurOut("The second way to use the read.dist command is to read a phylip or column and a group, so you can use the libshuff command.\n");
		mothurOut("For this use the read.dist command should be in the following format: \n");
		mothurOut("read.dist(phylip=yourPhylipfile, group=yourGroupFile). The cutoff and precision parameters are not valid with this use.  \n");
		mothurOut("Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourPhylipfile).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "ReadDistCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

ReadDistCommand::~ReadDistCommand(){
	if (abort == false) {
		if (format != "matrix") { 
			delete read; 
			delete nameMap; 
		}
	}
}

//**********************************************************************************************************************
int ReadDistCommand::execute(){
	try {
		
		if (abort == true) {	return 0;	}

		//time_t start = time(NULL);
		size_t numDists = 0;
		
		if (format == "matrix") {
			ifstream in;
			openInputFile(distFileName, in);
			matrix = new FullMatrix(in); //reads the matrix file
			in.close();
			//memory leak prevention
			if (globaldata->gMatrix != NULL) { delete globaldata->gMatrix;  }
			globaldata->gMatrix = matrix; //save matrix for coverage commands
			numDists = matrix->getSizes()[1];
		} else {
			read->read(nameMap);
			//to prevent memory leak

			if (globaldata->gListVector != NULL) {  delete globaldata->gListVector;  }
			globaldata->gListVector = read->getListVector();

			if (globaldata->gSparseMatrix != NULL) { delete globaldata->gSparseMatrix;  }
			globaldata->gSparseMatrix = read->getMatrix();
			numDists = globaldata->gSparseMatrix->getNNodes();
			
      int lines = cutoff / (1.0/precision);
      vector<float> dist_cutoff(lines+1,0);
			for (int i = 0; i <= lines;i++) {	
      	dist_cutoff[i] = (i + 0.5) / precision; 
      } 
      vector<int> dist_count(lines+1,0);
      list<PCell>::iterator currentCell;
      SparseMatrix* smatrix = globaldata->gSparseMatrix;
  		for (currentCell = smatrix->begin(); currentCell != smatrix->end(); currentCell++) {
				for (int i = 0; i <= lines;i++) {	
					if (currentCell->dist < dist_cutoff[i]) {
						dist_count[i]++;
            break;
          }
        }
			}

     // string dist_string = "Dist:";
    //  string count_string = "Count: ";
			//for (int i = 0; i <= lines;i++) {	
      	//dist_string = dist_string.append("\t").append(toString(dist_cutoff[i]));
      //	count_string = count_string.append("\t").append(toString(dist_count[i]));
		//	}
      //mothurOut(dist_string); mothurOutEndLine(); mothurOut(count_string); mothurOutEndLine();
		}
		//mothurOut("It took " + toString(time(NULL) - start) + " secs to read " + toString(numDists) + " distances (cutoff: " + toString(cutoff) + ")"); mothurOutEndLine();
		return 0;
		
	}
	catch(exception& e) {
		errorOut(e, "ReadDistCommand", "execute");
		exit(1);
	}
}
