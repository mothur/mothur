/*
 *  readotu.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readotu.h"

/***********************************************************************/

ReadOTUFile::ReadOTUFile(string pf): philFile(pf){
	
	//openInputFile(philFile, fileHandle);
}

/***********************************************************************/
//This function reads the list, rabund or sabund files to be used by collect and rarefact command.
void ReadOTUFile::read(GlobalData* globaldata){
	try {

		if (globaldata->getOrderFile() == "") {
			//you have two inputs because in the next if statement if you only have one then it moves ahead in the same file.  
			//So when you run the collect or summary commands you miss a line.
			input = new InputData(philFile, globaldata->getFormat()); //format tells you whether philFile is list, rabund, sabund.
			inputList = new InputData(philFile, globaldata->getFormat()); //format tells you whether philFile is list, rabund, sabund.
			inputSabund = new InputData(philFile, globaldata->getFormat()); //format tells you whether philFile is list, rabund, sabund or shared.
			inputRabund = new InputData(philFile, globaldata->getFormat());
		}else {//there is an orderfile
			input = new InputData(philFile, globaldata->getOrderFile(), globaldata->getFormat());
		}
	
		//memory leak prevention
		//if (globaldata->ginput != NULL) { delete globaldata->ginput;  }
		globaldata->ginput = input;	//saving to be used by collector and rarefact commands.
	
		if ((globaldata->getFormat() == "list") || (globaldata->getFormat() == "rabund") || (globaldata->getFormat() == "sabund")) {//you are reading a list, rabund or sabund file for collect, rarefaction or summary.

//cout << input << '\t' << globaldata << endl;
			order = input->getOrderVector();
			//memory leak prevention

			//if (globaldata->gorder != NULL) { delete globaldata->gorder;  }
			globaldata->gorder = order;	//saving to be used by collect and rarefact commands.
			sabund = inputSabund->getSAbundVector(); 
			//if (globaldata->sabund != NULL) { delete globaldata->sabund;  }
			globaldata->sabund = sabund; //saving to be used by summary command.
			delete inputSabund;

			rabund = inputRabund->getRAbundVector(); 
			//if (globaldata->rabund != NULL) { delete globaldata->rabund;  }
			globaldata->rabund = rabund; //saving to be used by heatmap.bin command.
			delete inputRabund;

			list = inputList->getListVector();
			//if (globaldata->gListVector != NULL) { delete globaldata->gListVector;  }
			globaldata->gListVector = list;
			delete inputList;

		}else if (globaldata->getFormat() == "shared") {
			SharedList = input->getSharedListVector(); //you are reading for collect.shared, rarefaction.shared, summary.shared, parselist command, or shared commands.
			//memory leak prevention
			//if (globaldata->gSharedList != NULL) { delete globaldata->gSharedList;  }
			globaldata->gSharedList = SharedList;
		}
	}
	catch(exception& e) {
		errorOut(e, "ReadOTUFile", "read");
		exit(1);
	}
}

/***********************************************************************/

ReadOTUFile::~ReadOTUFile(){
//	delete input;
//	delete order;
}

/***********************************************************************/

