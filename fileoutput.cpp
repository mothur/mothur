/*
 *  fileoutput.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "fileoutput.h"

/***********************************************************************/

ThreeColumnFile::~ThreeColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());
};

/***********************************************************************/

void ThreeColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);

			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << '\t' << label << "\tlci\thci" << endl;
		}
		else{
			openOutputFile(outName, outFile);
			outFile << "numsampled\t" << label << "\tlci\thci" << endl;
		}

		outFile.setf(ios::fixed, ios::floatfield);
		outFile.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ThreeColumnFile class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ThreeColumnFile class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void ThreeColumnFile::output(int nSeqs, vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
		}
		else{
			outFile	<< nSeqs << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ThreeColumnFile class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ThreeColumnFile class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************/

void ThreeColumnFile::resetFile(){
	try {
		if(counter != 0){
			outFile.close();
			inFile.close();
		}
		else{
			outFile.close();
		}
		counter = 1;
		
		remove(inName.c_str());
		renameOk = rename(outName.c_str(), inName.c_str());
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) {	cout << "Unable to rename necessary files." << endl; }

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ThreeColumnFile class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ThreeColumnFile class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
/***********************************************************************/

ColumnFile::~ColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());
};

/***********************************************************************/

void ColumnFile::initFile(string label, vector<string> tags){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);

			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << '\t'; 
			for(int i = 0; i < tags.size(); i++) {
				outFile << label + tags[i] << '\t';
			}
			outFile << endl;
		}
		else{
			openOutputFile(outName, outFile);
			for(int i = 0; i < tags.size(); i++) {
				outFile << label + tags[i] << '\t';
			}
			outFile << endl;
		}

		outFile.setf(ios::fixed, ios::floatfield);
		outFile.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ColumnFile class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ColumnFile class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void ColumnFile::output(vector<double> data){
	try {
	
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);

			outFile << inputBuffer << '\t' << setprecision(6) << data[0] << setprecision(globaldata->getIters().length());
			for (int i = 1; i< data.size(); i++) {
				outFile << '\t' << data[i]; 
			}
			outFile << endl;
		}
		else{
			outFile << setprecision(6) << data[0] << setprecision(globaldata->getIters().length());
			for (int i = 1; i< data.size(); i++) {
				outFile << '\t' << data[i]; 
			}
			outFile << endl;
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ColumnFile class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ColumnFile class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************/

void ColumnFile::resetFile(){
	try {
		if(counter != 0){
			outFile.close();
			inFile.close();
		}
		else{
			outFile.close();
		}
		counter = 1;
		
		remove(inName.c_str());
		renameOk = rename(outName.c_str(), inName.c_str());
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { cout << "Unable to rename necessary files." << endl; }

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ColumnFile class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ColumnFile class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
/***********************************************************************/

SharedThreeColumnFile::~SharedThreeColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());
};

/***********************************************************************/

void SharedThreeColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);

			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << '\t' << label << "\tlci\thci" << endl;
		}
		else{
			openOutputFile(outName, outFile);
			outFile << "numsampled\t" << groupLabel << '\t' << label << "\tlci\thci" << endl;
		}

		outFile.setf(ios::fixed, ios::floatfield);
		outFile.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedThreeColumnFile class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedThreeColumnFile class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void SharedThreeColumnFile::output(int nSeqs, vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
		}
		else{
			outFile	<< numGroup << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
			numGroup++;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedThreeColumnFile class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedThreeColumnFile class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************/

void SharedThreeColumnFile::resetFile(){
	try {
		if(counter != 0){
			outFile.close();
			inFile.close();
		}
		else{
			outFile.close();
		}
		counter = 1;
		remove(inName.c_str());
		renameOk = rename(outName.c_str(), inName.c_str());
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { cout << "Unable to rename necessary files." << endl; }

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedThreeColumnFile class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedThreeColumnFile class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

/***********************************************************************/

OneColumnFile::~OneColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());	
};

/***********************************************************************/

void OneColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);
		
			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << '\t' << label << endl;
		}
		else{
			openOutputFile(outName, outFile);
			outFile << "numsequences\t" << label << endl;
		}
	
		outFile.setf(ios::fixed, ios::floatfield);
		outFile.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OneColumnFile class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OneColumnFile class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

void OneColumnFile::output(int nSeqs, vector<double> data){
	try {	
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << setprecision(4) << '\t'  << data[0] << endl;
		}
		else{	
			outFile	<< nSeqs << setprecision(4) << '\t' << data[0] << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OneColumnFile class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OneColumnFile class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
};

/***********************************************************************/

void OneColumnFile::resetFile(){
	try {
		if(counter != 0){
			outFile.close();
			inFile.close();
		}
		else{
			outFile.close();
		}	
		counter = 1;
		remove(inName.c_str());
		renameOk = rename(outName.c_str(), inName.c_str());
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { cout << "Unable to rename necessary files." << endl; }

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OneColumnFile class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OneColumnFile class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
/***********************************************************************/

SharedOneColumnFile::~SharedOneColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());	
};

/***********************************************************************/

void SharedOneColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);
		
			string inputBuffer;
			getline(inFile, inputBuffer);
		
			outFile	<<  inputBuffer << '\t' << label  << endl;

		}
		else{
			openOutputFile(outName, outFile);
			outFile << "sampled\t" << label << endl;
		
		}
	
		outFile.setf(ios::fixed, ios::floatfield);
		outFile.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OneColumnFile class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OneColumnFile class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

void SharedOneColumnFile::output(int nSeqs, vector<double> data){
	try {	
			string dataOutput;
			float sam;
			sam = data[0];
			dataOutput = "";
			for (int i = 0; i < data.size(); i++) {
				dataOutput = dataOutput + "\t" + toString(data[i]);
			}
			if(counter != 0){		
				string inputBuffer;
				getline(inFile, inputBuffer);

				outFile	<<  inputBuffer << setprecision(2) << '\t' << dataOutput << endl;
			}
			else{	
				outFile	<< nSeqs << setprecision(2) << '\t' << dataOutput << endl;
			}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OneColumnFile class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OneColumnFile class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
};

/***********************************************************************/

void SharedOneColumnFile::resetFile(){
	try {
		if(counter != 0){
			outFile.close();
			inFile.close();
		}
		else{
			outFile.close();
		}	
		counter = 1;
		
		remove(inName.c_str());
		renameOk = rename(outName.c_str(), inName.c_str());
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { cout << "Unable to rename necessary files." << endl; }

		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OneColumnFile class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OneColumnFile class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
