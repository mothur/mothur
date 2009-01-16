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
		rename(outName.c_str(), inName.c_str());
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
		rename(outName.c_str(), inName.c_str());
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
		rename(outName.c_str(), inName.c_str());
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

SharedOneColumnFile::~OneColumnFile(){
	
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
		rename(outName.c_str(), inName.c_str());
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
