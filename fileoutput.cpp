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
}

/***********************************************************************/

void ThreeColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);

			string inputBuffer;
			inputBuffer = getline(inFile);
		
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
		errorOut(e, "ThreeColumnFile", "initFile");
		exit(1);
	}
}

/***********************************************************************/

void ThreeColumnFile::output(int nSeqs, vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			inputBuffer = getline(inFile);
		
			outFile	<<  inputBuffer << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
		}
		else{
			outFile	<< nSeqs << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
		}
	}
	catch(exception& e) {
		errorOut(e, "ThreeColumnFile", "output");
		exit(1);
	}
}

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
		
		//renameFile(outName, inName);
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { 
			mothurOut("Unable to rename " + outName); mothurOutEndLine();
			perror(" : ");
		}	
	}
	catch(exception& e) {
		errorOut(e, "ThreeColumnFile", "resetFile");
		exit(1);
	}
}

/***********************************************************************/
/***********************************************************************/

ColumnFile::~ColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());
}

/***********************************************************************/

void ColumnFile::initFile(string label, vector<string> tags){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);

			string inputBuffer;
			inputBuffer = getline(inFile);
		
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
		errorOut(e, "ColumnFile", "initFile");
		exit(1);
	}
}

/***********************************************************************/

void ColumnFile::output(vector<double> data){
	try {
	
		if(counter != 0){		
			string inputBuffer;
			inputBuffer = getline(inFile);

			outFile << inputBuffer << '\t' << setprecision(6) << data[0] << setprecision(iters.length());
			for (int i = 1; i< data.size(); i++) {
				outFile << '\t' << data[i]; 
			}
			outFile << endl;
		}
		else{
			outFile << setprecision(6) << data[0] << setprecision(iters.length());
			for (int i = 1; i< data.size(); i++) {
				outFile << '\t' << data[i]; 
			}
			outFile << endl;
		}

	}
	catch(exception& e) {
		errorOut(e, "ColumnFile", "output");
		exit(1);
	}
}

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
		
		//renameFile(outName, inName);
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { 
			mothurOut("Unable to rename " + outName); mothurOutEndLine();
			perror(" : ");
		}	
	}
	catch(exception& e) {
		errorOut(e, "ColumnFile", "resetFile");
		exit(1);
	}
}

/***********************************************************************/
/***********************************************************************/

SharedThreeColumnFile::~SharedThreeColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());
}

/***********************************************************************/

void SharedThreeColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);

			string inputBuffer;
			inputBuffer = getline(inFile);
		
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
		errorOut(e, "SharedThreeColumnFile", "initFile");
		exit(1);
	}
}

/***********************************************************************/

void SharedThreeColumnFile::output(int nSeqs, vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			inputBuffer = getline(inFile);
		
			outFile	<<  inputBuffer << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
		}
		else{
			outFile	<< numGroup << setprecision(4) << '\t' << data[0] << '\t' << data[1] << '\t' << data[2] << endl;
			numGroup++;
		}
	}
	catch(exception& e) {
		errorOut(e, "SharedThreeColumnFile", "output");
		exit(1);
	}
}

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
		
		//renameFile(outName, inName);
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { 
			mothurOut("Unable to rename " + outName); mothurOutEndLine();
			perror(" : ");
		}	
	}
	catch(exception& e) {
		errorOut(e, "SharedThreeColumnFile", "resetFile");
		exit(1);
	}
}

/***********************************************************************/

/***********************************************************************/

OneColumnFile::~OneColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());	
}

/***********************************************************************/

void OneColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);
		
			string inputBuffer;
			inputBuffer = getline(inFile);
		
			outFile	<<  inputBuffer << '\t' << label << endl;
		}
		else{
			openOutputFile(outName, outFile);
			outFile << "numsampled\t" << label << endl;
		}
	
		outFile.setf(ios::fixed, ios::floatfield);
		outFile.setf(ios::showpoint);
	}
	catch(exception& e) {
		errorOut(e, "OneColumnFile", "initFile");
		exit(1);
	}
}

/***********************************************************************/

void OneColumnFile::output(int nSeqs, vector<double> data){
	try {	
		if(counter != 0){		
			string inputBuffer;
			inputBuffer = getline(inFile);
		
			outFile	<<  inputBuffer << setprecision(4) << '\t'  << data[0] << endl;
		}
		else{	
			outFile	<< nSeqs << setprecision(4) << '\t' << data[0] << endl;
		}
	}
	catch(exception& e) {
		errorOut(e, "OneColumnFile", "output");
		exit(1);
	}
}

/***********************************************************************/

void OneColumnFile::resetFile(){
	try {
		if(counter != 0){
			outFile.close();
			inFile.close();
		}else{
			outFile.close();
		}	
		counter = 1;
		
		remove(inName.c_str());
		renameOk = rename(outName.c_str(), inName.c_str());
		
		//renameFile(outName, inName);
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { 
			mothurOut("Unable to rename " + outName); mothurOutEndLine();
			perror(" : ");
		}	

	}
	catch(exception& e) {
		errorOut(e, "OneColumnFile", "resetFile");
		exit(1);
	}
}

/***********************************************************************/
/***********************************************************************/

SharedOneColumnFile::~SharedOneColumnFile(){
	
	inFile.close();
	outFile.close();
	remove(outName.c_str());	
}

/***********************************************************************/

void SharedOneColumnFile::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(outName, outFile);
			openInputFile(inName, inFile);
		
			string inputBuffer;
			inputBuffer = getline(inFile);
		
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
		errorOut(e, "SharedOneColumnFile", "initFile");
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
				inputBuffer = getline(inFile);

				outFile	<<  inputBuffer << setprecision(2) << '\t' << dataOutput << endl;
			}
			else{	
				outFile	<< nSeqs << setprecision(2) << '\t' << dataOutput << endl;
			}
	}
	catch(exception& e) {
		errorOut(e, "SharedOneColumnFile", "output");
		exit(1);
	}
}

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
		
		//renameFile(outName, inName);
		
		//checks to make sure user was able to rename and remove successfully
		if ((renameOk != 0)) { 
			mothurOut("Unable to rename " + outName); mothurOutEndLine();
			perror(" : ");
		}	
	}
	catch(exception& e) {
		errorOut(e, "SharedOneColumnFile", "resetFile");
		exit(1);
	}
}

/***********************************************************************/
