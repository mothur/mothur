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
void ThreeColumnFile::setLabelName(string label){
	try {
		if(!firstLabel)     { fileHeader += "\t" + label + "\tlci\thci";             }
        else                { fileHeader = "numsampled\t" + label + "\tlci\thci";    }
    }
	catch(exception& e) {
		m->errorOut(e, "ThreeColumnFile", "setLabelName");
		exit(1);
	}
}
/***********************************************************************/
void ThreeColumnFile::updateOutput(int nSeqs, vector<double> data){
	try {
        
        map<int, int>::iterator it = nseqsToRow.find(nSeqs);
        if (it != nseqsToRow.end()) { //new column in row for next label
            int resultsIndex = it->second;
            results[resultsIndex].push_back(data[0]);
            results[resultsIndex].push_back(data[1]);
            results[resultsIndex].push_back(data[2]);
        }else{
            //new row
            nseqsToRow[nSeqs] = results.size();
            vector<double> theseResults;
            theseResults.push_back(nSeqs);
            theseResults.push_back(data[0]);
            theseResults.push_back(data[1]);
            theseResults.push_back(data[2]);
            results.push_back(theseResults);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "ThreeColumnFile", "updateOutput");
		exit(1);
	}
}
/***********************************************************************/
void ColumnFile::setLabelName(string label, vector<string> tags){
	try {
		if(firstLabel){ fileHeader = ""; }
        
        for(int i = 0; i < tags.size(); i++) { fileHeader += label + tags[i] + '\t'; }
	}
	catch(exception& e) {
		m->errorOut(e, "ColumnFile", "setLabelName");
		exit(1);
	}
}
/***********************************************************************/
void ColumnFile::updateOutput(vector<double> data){
	try {
        vector<double> theseResults;
        for (size_t i = 0; i < data.size(); i++) { theseResults.push_back(data[i]);  }
        results.push_back(theseResults);
	}
	catch(exception& e) {
		m->errorOut(e, "ColumnFile", "updateOutput");
		exit(1);
	}
}
/***********************************************************************/
void FileOutput::printFile(){
	try {
        ofstream outFile;
        util.openOutputFile(filename, outFile);
        
        outFile.setf(ios::fixed, ios::floatfield); outFile.setf(ios::showpoint);
        cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
    
        outFile << fileHeader << endl;
        for (size_t i = 0; i < results.size(); i++) {
            for (size_t j = 0; j < results[i].size(); j++) {
                outFile << setprecision(6) << results[i][j] << '\t';
            }
            outFile << endl;
        }
        outFile << endl; 
        
        outFile.close();
    }
	catch(exception& e) {
		m->errorOut(e, "FileOutput", "printFile");
		exit(1);
	}
}
/***********************************************************************/
void SharedThreeColumnFile::setLabelName(string label){
	try {
        if(!firstLabel)     { fileHeader += "\t" + label + "\tlci\thci";                                }
        else                { fileHeader = "numsampled\t" + groupLabel + "\t" + label + "\tlci\thci";   }
	}
	catch(exception& e) {
		m->errorOut(e, "SharedThreeColumnFile", "setLabelName");
		exit(1);
	}
}
/***********************************************************************/
void SharedThreeColumnFile::updateOutput(int nSeqs, vector<double> data){
	try {
        
        map<int, int>::iterator it = nseqsToRow.find(nSeqs);
        if (it != nseqsToRow.end()) { //new column in row for next label
            int resultsIndex = it->second;
            results[resultsIndex].push_back(data[0]);
            results[resultsIndex].push_back(data[1]);
            results[resultsIndex].push_back(data[2]);
        }else{ //new row
            nseqsToRow[nSeqs] = results.size();
            vector<double> theseResults;
            theseResults.push_back(numGroup); numGroup++;
            theseResults.push_back(nSeqs);
            theseResults.push_back(data[0]);
            theseResults.push_back(data[1]);
            theseResults.push_back(data[2]);
            results.push_back(theseResults);
        }
    }
	catch(exception& e) {
		m->errorOut(e, "SharedThreeColumnFile", "output");
		exit(1);
	}
}
/***********************************************************************/
void OneColumnFile::setLabelName(string label){
	try {
        if(!firstLabel)     { fileHeader += "\t" + label;            }
        else                { fileHeader = "numsampled\t" + label;   }
    }
	catch(exception& e) {
		m->errorOut(e, "OneColumnFile", "setLabelName");
		exit(1);
	}
}
/***********************************************************************/
void OneColumnFile::updateOutput(int nSeqs, vector<double> data){
	try {
        map<int, int>::iterator it = nseqsToRow.find(nSeqs);
        if (it != nseqsToRow.end()) { //new column in row for next label
            int resultsIndex = it->second;
            results[resultsIndex].push_back(data[0]);
        }else{
            //new row
            nseqsToRow[nSeqs] = results.size();
            vector<double> theseResults;
            theseResults.push_back(nSeqs);
            theseResults.push_back(data[0]);
            results.push_back(theseResults);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "OneColumnFile", "updateOutput");
		exit(1);
	}
}
/***********************************************************************/

void SharedOneColumnFile::setLabelName(string label){
	try {
        if(!firstLabel)     { fileHeader += "\t" + label;            }
        else                { fileHeader = "sampled\t" + label;   }
    }
	catch(exception& e) {
		m->errorOut(e, "SharedOneColumnFile", "setLabelName");
		exit(1);
	}
}
/***********************************************************************/
void SharedOneColumnFile::updateOutput(int nSeqs, vector<double> data){
	try {
        
        map<int, int>::iterator it = nseqsToRow.find(nSeqs);
        if (it != nseqsToRow.end()) { //new column in row for next label
            int resultsIndex = it->second;
            for (int i = 0; i < data.size(); i++) { results[resultsIndex].push_back(data[i]); }
        }else{
            //new row
            nseqsToRow[nSeqs] = results.size();
            vector<double> theseResults;
            theseResults.push_back(nSeqs);
            for (int i = 0; i < data.size(); i++) { theseResults.push_back(data[i]); }
            results.push_back(theseResults);
        }
    }
	catch(exception& e) {
		m->errorOut(e, "SharedOneColumnFile", "updateOutput");
		exit(1);
	}
}
/***********************************************************************/


