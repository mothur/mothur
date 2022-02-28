/*
 *  mothurOut.cpp
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"
#include "ordervector.hpp"
#include "sharedordervector.h"
#include "counttable.h"

/******************************************************/
MothurOut* MothurOut::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new MothurOut();
	}
	return _uniqueInstance;
}
/*********************************************************************************************/
void MothurOut::appendLogBuffer(string partialLog)  {
    try {
        buffer += partialLog;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "appendLogBuffer");
        exit(1);
    }
}
/*********************************************************************************************/
void MothurOut::setLogFileName(string filename, bool append)  {
	try {
        silenceWarnings = false;
        Utils util;
        if ((filename == "silent")) { silenceLog = true; }
        else {
            logFileName = filename;
            if (outLog != NULL) {
                closeLog();
                delete outLog; outLog = NULL;
            }
            outLog = new ofstream();
            silenceLog = false;
            if (append)     {
                util.openOutputFileAppend(filename, *outLog);
                *outLog << "\n\n************************************************************\n\n\n";
            }else            {
                bool opendLog = util.openOutputFile(filename, *outLog);
                if (!opendLog) { control_pressed = true; }
            }
        }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFileName");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::closeLog()  {
	try {
        if (buffer != "") { string output = buffer; buffer = ""; mothurOut(output);   }
        
        string outputLogName = "Logfile : " + logFileName + "\n\n";
        if (!silenceLog) { mothurOut(outputLogName); }
            
        if (numErrors != 0) {
            if (!silenceLog) {
                *outLog << "\n\n************************************************************\n";
                *outLog << "************************************************************\n";
                *outLog << "************************************************************\n";
                *outLog << "Detected " + toString(numErrors) + " [ERROR] messages, please review.\n";
                *outLog << "************************************************************\n";
                *outLog << "************************************************************\n";
                *outLog << "************************************************************\n";
            }
            logger() << "\n\n************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "Detected " + toString(numErrors) + " [ERROR] messages, please review.\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
        }
        
        if (numWarnings != 0) {
            if (!silenceLog) {
                *outLog << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                *outLog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                *outLog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                *outLog << "Detected " + toString(numWarnings) + " [WARNING] messages, please review.\n";
                *outLog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                *outLog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                *outLog << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            }
            logger() << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "Detected " + toString(numWarnings) + " [WARNING] messages, please review.\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
        }
        
        outLog->close();
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "closeLog");
		exit(1);
	}
}
/*********************************************************************************************/
MothurOut::~MothurOut() {
	try {
		_uniqueInstance = 0;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOut(string output) {
	try {
        if (outLog == NULL) { appendLogBuffer(output);  return; }
        
        if (buffer != "") { output = buffer + output; buffer = ""; }
        if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        bool savedSilenceLog = silenceLog;
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {
                    logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
        if (!quietMode) {
            if (!silenceLog) {
                if (silenceWarnings && containsWarning) {} //do not print warning to logfile if warnings are silenced
                else { *outLog << output;  }
            }
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                if (!silenceLog) { *outLog << output; }
                logger() << output;
            }
        }
        silenceLog = savedSilenceLog;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToScreen(string output) {
	try {
        if (buffer != "") { output = buffer + output; buffer = ""; }
		if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {  logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n"; }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
        if (!quietMode) {
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                logger() << output;
            }
        }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine() {
	try {
        if (outLog == NULL) { appendLogBuffer("\n"); return; }
        
		if (!quietMode) {
            if (!silenceLog) { *outLog << buffer << endl; }
            logger() << buffer << endl;
        }
        buffer = "";
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToLog(string output) {
	try {
        if (outLog == NULL) { appendLogBuffer(output); return; }
        
        if (buffer != "") { output = buffer + output; buffer = ""; }
        if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        
        bool savedSilenceLog = silenceLog;
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {
                    logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
        if (!quietMode) {
            if (!silenceLog) {
                if (silenceWarnings && containsWarning) {} //do not print warning to logfile if warnings are silenced
                else { *outLog << output; }
            }
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                if (!silenceLog) { *outLog << output; }
            }
        }
        silenceLog = savedSilenceLog;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutJustToLog");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::errorOut(exception& e, string object, string function) {
    numErrors++; 
	
    string errorType = toString(e.what());
    
    int pos = errorType.find("bad_alloc");
    mothurOut("[ERROR]: " + errorType);
    
    double ramUsed, total;
    Utils util;
    ramUsed = util.getRAMUsed(); total = util.getTotalRAM();
    mothurOut("RAM used: " + toString(ramUsed/(double)GIG) + "Gigabytes . Total Ram: " + toString(total/(double)GIG) + "Gigabytes.\n\n");
    
    if (pos == string::npos) { //not bad_alloc
        mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry\n");
    }else { //bad alloc
        if (object == "cluster"){
            mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory.  There are two common causes for this, file size and format.\n\nFile Size:\nThe cluster command loads your distance matrix into RAM, and your distance file is most likely too large to fit in RAM. There are two options to help with this. The first is to use a cutoff. By using a cutoff mothur will only load distances that are below the cutoff. If that is still not enough, there is a command called cluster.split, http://www.mothur.org/wiki/cluster.split which divides the distance matrix, and clusters the smaller pieces separately. You may also be able to reduce the size of the original distance matrix by using the commands outlined in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP. \n\nWrong Format:\nThis error can be caused by trying to read a column formatted distance matrix using the phylip parameter. By default, the dist.seqs command generates a column formatted distance matrix. To make a phylip formatted matrix set the dist.seqs command parameter output to lt.  \n\nIf you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        }else if (object == "shhh.flows"){
                mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory. The shhh.flows command is very memory intensive. This error is most commonly caused by trying to process a dataset too large, using multiple processors, or failing to run trim.flows before shhh.flows. If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Running trim.flows with an oligos file, and then shhh.flows with the file option may also resolve the issue. If for some reason you are unable to run shhh.flows with your data, a good alternative is to use the trim.seqs command using a 50-bp sliding window and to trim the sequence when the average quality score over that window drops below 35. Our results suggest that the sequencing error rates by this method are very good, but not quite as good as by shhh.flows and that the resulting sequences tend to be a bit shorter. If you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry. ");
        }else {
            mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory.  This is most commonly caused by trying to process a dataset too large, using multiple processors, or a file format issue. If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Also, you may be able to reduce the size of your dataset by using the commands outlined in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP. If you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        }
    }
}
/*********************************************************************************************/
void MothurOut::setHomePath(string pathname)  {
    try {
        if (pathname != "") {
            //add / to name if needed
            string lastChar = pathname.substr(pathname.length()-1);
            if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
        }
        homePath = pathname;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "setHomePath");
        exit(1);
    }
}
/*********************************************************************************************/
void MothurOut::setPaths(vector<string> pathVariables)  {
    try {
        for (int i = 0; i < pathVariables.size(); i++) {
            string pathname = pathVariables[i];
            if (pathname != "") {
                //add / to name if needed
                string lastChar = pathname.substr(pathname.length()-1);
                if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
            }
        }
        
        paths = pathVariables;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "setPaths");
        exit(1);
    }
}
/*********************************************************************************************/
void MothurOut::initialize()  {
    try {
        
        validAminoAcids.insert('A');
        validAminoAcids.insert('R');
        validAminoAcids.insert('N');
        validAminoAcids.insert('D');
        
        validAminoAcids.insert('B');
        validAminoAcids.insert('C');
        validAminoAcids.insert('Q');
        validAminoAcids.insert('E');
        
        validAminoAcids.insert('Z');
        validAminoAcids.insert('G');
        validAminoAcids.insert('H');
        validAminoAcids.insert('I');
        
        validAminoAcids.insert('L');
        validAminoAcids.insert('K');
        validAminoAcids.insert('M');
        validAminoAcids.insert('F');
        
        validAminoAcids.insert('P');
        validAminoAcids.insert('S');
        validAminoAcids.insert('T');
        validAminoAcids.insert('W');
        
        validAminoAcids.insert('Y');
        validAminoAcids.insert('V');
        validAminoAcids.insert('X');
        validAminoAcids.insert('-');
        validAminoAcids.insert('.');
        
        validAminoAcids.insert('*');
        validAminoAcids.insert('?');
        
        codons.clear(); codons.resize(4);
        for (int i = 0; i < codons.size(); i++) {
            codons[i].resize(4);
            for (int j = 0; j < codons[i].size(); j++) {
                codons[i][j].resize(4);
            }
        }
                
        //AAX
        codons[0][0][0] = 'K';   //AAA |  Lysine (K) -> 11. where 11 is the index into the aas enum.
        codons[0][0][1] = 'N';   //AAT |  Asparagine (N) -> 2.
        codons[0][0][2] = 'K';   //AAG |  Lysine (K) -> 11.
        codons[0][0][3] = 'N';   //AAC |  Asparagine (N) -> 2.
        
        //ATX
        codons[0][1][0] = 'I';   //ATA |  Isoleucine (I) -> 9.
        codons[0][1][1] = 'I';   //ATT |  Isoleucine (I) -> 9.
        codons[0][1][2] = 'M';   //ATG |  Methionine (M) -> 12.
        codons[0][1][3] = 'I';   //ATC |  Isoleucine (I) -> 9.
        
        //AGX
        codons[0][2][0] = 'R';   //AGA |  Arginine (R) -> 1.
        codons[0][2][1] = 'S';   //AGT |  Serine (S) -> 15.
        codons[0][2][2] = 'R';   //AGG |  Arginine (R) -> 1.
        codons[0][2][3] = 'S';   //AGC |  Serine (S) -> 15.
        
        //ACX
        codons[0][3][0] = 'T';    //ACA |  Threonine (T) -> 17.
        codons[0][3][1] = 'T';    //ACT |  Threonine (T) -> 17.
        codons[0][3][2] = 'T';    //ACG |  Threonine (T) -> 17.
        codons[0][3][3] = 'T';    //ACC |  Threonine (T) -> 17.
        
        
        //TAX
        codons[1][0][0] = '*';   //TAA | Termination (X) -> 22
        codons[1][0][1] = 'Y';   //TAT | Tyrosine (Y) -> 19
        codons[1][0][2] = '*';   //TAG | Termination (X) -> 22
        codons[1][0][3] = 'Y';   //TAC | Tyrosine (Y) -> 19
        
        //TTX
        codons[1][1][0] = 'L';    //TTA | Leucine (L) -> 10
        codons[1][1][1] = 'F';    //TTT | Phenylalanine (F) -> 13
        codons[1][1][2] = 'L';    //TTG | Leucine (L) -> 10
        codons[1][1][3] = 'F';    //TTC | Phenylalanine (F) -> 13
        
        //TGX
        codons[1][2][0] = '*';    //TGA | Termination (X) -> 22
        codons[1][2][1] = 'C';    //TGT | Cysteine (C) -> 4
        codons[1][2][2] = 'W';    //TGG | Tryptophan (W) -> 18
        codons[1][2][3] = 'C';    //TGC | Cysteine (C) -> 4
        
        //TCX
        codons[1][3][0] = 'S';    //TCA | Serine (S) -> 15
        codons[1][3][1] = 'S';    //TCT | Serine (S) -> 15
        codons[1][3][2] = 'S';    //TCG | Serine (S) -> 15
        codons[1][3][3] = 'S';    //TCC | Serine (S) -> 15
        
        //GAX
        codons[2][0][0] = 'E';   //GAA | Glutamate (E) -> 6
        codons[2][0][1] = 'D';   //GAT | Aspartate (D) -> 3
        codons[2][0][2] = 'E';   //GAG | Glutamate (E) -> 6
        codons[2][0][3] = 'D';   //GAC | Aspartate (D) -> 3
        
        //GTX
        codons[2][1][0] = 'V';    //GTA | Valine (V)
        codons[2][1][1] = 'V';    //GTT | Valine (V)
        codons[2][1][2] = 'V';    //GTG | Valine (V)
        codons[2][1][3] = 'V';    //GTC | Valine (V)
        
        //GGX
        codons[2][2][0] = 'G';    //GGA | Glycine (G)
        codons[2][2][1] = 'G';    //GGT | Glycine (G)
        codons[2][2][2] = 'G';    //GGG | Glycine (G)
        codons[2][2][3] = 'G';    //GGC | Glycine (G)
        
        //GCX
        codons[2][3][0] = 'A';    //GCA | Alanine (A)
        codons[2][3][1] = 'A';    //GCT | Alanine (A)
        codons[2][3][2] = 'A';    //GCG | Alanine (A)
        codons[2][3][3] = 'A';    //GCC | Alanine (A)
        
        //CAX
        codons[3][0][0] = 'Q';   //CAA | Glutamine (Q)
        codons[3][0][1] = 'H';   //CAT | Histidine (H)
        codons[3][0][2] = 'Q';   //CAG | Glutamine (Q)
        codons[3][0][3] = 'H';   //CAC | Histidine (H)
        
        //CTX
        codons[3][1][0] = 'L';    //CTA | Leucine (L)
        codons[3][1][1] = 'L';    //CTT | Leucine (L)
        codons[3][1][2] = 'L';    //CTG | Leucine (L)
        codons[3][1][3] = 'L';    //CTC | Leucine (L)
        
        //CGX
        codons[3][2][0] = 'R';    //CGA | Arginine (R)
        codons[3][2][1] = 'R';    //CGT | Arginine (R)
        codons[3][2][2] = 'R';    //CGG | Arginine (R)
        codons[3][2][3] = 'R';    //CGC | Arginine (R)
        
        //CCX
        codons[3][3][0] = 'P';    //CCA | Proline (P)
        codons[3][3][1] = 'P';    //CCT | Proline (P)
        codons[3][3][2] = 'P';    //CCG | Proline (P)
        codons[3][3][3] = 'P';    //CCC | Proline (P)
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "initialize");
        exit(1);
    }
}
/********************************************************************/


