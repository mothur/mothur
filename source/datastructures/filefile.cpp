//
//  filefile.cpp
//  Mothur
//
//  Created by Sarah Westcott on 12/5/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "filefile.hpp"
#include "utils.hpp"

/**************************************************************************************************/

FileFile::FileFile(string f, string md) : filename(f), mode(md) {
    try {
        m = MothurOut::getInstance();
        
        current = CurrentFile::getInstance();
        inputDir = current->getInputDir();
        mpath = current->getProgramPath();
        
        columnWithGroups = false;
        fileOption = 0;
        gz = false;
        hasIndex = false;
        read(f, mode);
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "FileFile");
        exit(1);
    }
}
/**************************************************************************************************

FileFile::FileFile(){
    try {
        m = MothurOut::getInstance();
        
        current = CurrentFile::getInstance();
        inputDir = current->getInputDir();
        mpath = current->getProgramPath();
        
        fileOption = 0;
        gz = false;
        hasIndex = false;
        filename = ""; mode = "";
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "FileFile");
        exit(1);
    }
}
/*
 file option 1
 
 sfffile1   oligosfile1
 sfffile2   oligosfile2
 ...
 
 file option 2
 
 fastqfile1 oligosfile1
 fastqfile2 oligosfile2
 ...
 
 file option 3
 
 ffastqfile1 rfastqfile1
 ffastqfile2 rfastqfile2
 ...
 
 file option 4
 
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 ...
 
 file option 5
 
 My.forward.fastq My.reverse.fastq none My.rindex.fastq //none is an option is no forward or reverse index file
 ...
 */

/**************************************************************************************************/
vector< vector<string> > FileFile::read(string f, string mode){
    try {
        filename = f;
        
        bool allGZ = true; bool allPlainTxt = true;
        
        ifstream in;
        util.openInputFile(filename, in);
        
        while(!in.eof()) {
            
            if (m->getControl_pressed()) { return files; }
            
            bool skip = false;
            string line = util.getline(in);  util.gobble(in);
            
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: " + line +"\n");  }
            
            if(line[0] == '#'){  } //ignore
            else {
                
                vector<string> pieces = util.splitWhiteSpace(line);
                
                string group = ""; string forward, reverse, findex, rindex;
                skip = validateFiles(pieces, forward, reverse, findex, rindex, group);
                
                if (!skip) { //good pair
                    file2Group[files.size()] = group;
                    if (((findex != "") || (rindex != ""))) { hasIndex = true; }
                    
                    if ((mode == "contigs") || (mode == "sra")) { setGZ(forward, reverse, findex, rindex, allGZ, allPlainTxt); }
                    
                    vector<string> pair;
                    pair.push_back(forward); pair.push_back(reverse); pair.push_back(findex); pair.push_back(rindex);
                    files.push_back(pair);
                }
            }
        }
        in.close();
        
        if ((mode == "contigs") || (mode == "sra")){ if (allGZ) { gz = true; } else { gz = false; } }
        
        if (files.size() == 0) { m->setControl_pressed(true); }
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "read");
        exit(1);
    }
}
/**************************************************************************************************/
bool FileFile::validateFiles(vector<string> pieces, string& forward, string& reverse, string& findex, string& rindex, string& group){
    try {
        bool skip = false; //innocent until proven guilty
        group = "";
        if (pieces.size() == 2) {
            if (mode == "parsefastqpacbio") {
                group = pieces[0];
                util.checkGroupName(group);
                forward = pieces[1];
                reverse = "";
            }else {
                forward = pieces[0];
                reverse = pieces[1];
                group = "";
            }
            findex = "";
            rindex = "";
            fileOption = 1;
        }else if (pieces.size() == 3) {
            group = pieces[0];
            util.checkGroupName(group);
            forward = pieces[1];
            reverse = pieces[2];
            findex = "";
            rindex = "";
            fileOption = 2;
            columnWithGroups = true;
        }else if (pieces.size() == 4) {
            forward = pieces[0];
            reverse = pieces[1];
            findex = pieces[2];
            rindex = pieces[3];
            if ((findex == "none") || (findex == "NONE")){ findex = "NONE"; }
            if ((rindex == "none") || (rindex == "NONE")){ rindex = "NONE"; }
            fileOption = 3;
        }else {
            m->mothurOut("[ERROR]: file lines can be 2, 3, or 4 columns. The forward fastq files in the first column and their matching reverse fastq files in the second column, or a groupName then forward fastq file and reverse fastq file, or forward fastq file then reverse fastq then forward index and reverse index file.  If you only have one index file add 'none' for the other one. \n"); m->setControl_pressed(true);
        }
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: group = " + group + ", forward = " + forward + ", reverse = " + reverse + ", forwardIndex = " + findex + ", reverseIndex = " + rindex + ".\n"); }
        
        if (inputDir != "") {
            string path = util.hasPath(forward);
            if (path == "") {  forward = inputDir + forward;  }
            
            if (reverse != "") { //could be blank if pacbio mode
                path = util.hasPath(reverse);
                if (path == "") {  reverse = inputDir + reverse;  }
            }
            
            if (findex != "") {
                path = util.hasPath(findex);
                if (path == "") {  findex = inputDir + findex;  }
            }
            
            if (rindex != "") {
                path = util.hasPath(rindex);
                if (path == "") {  rindex = inputDir + rindex;  }
            }
        }
        
        //check to make sure both are able to be opened
        bool openForward = util.checkLocations(forward, current->getLocations());
        if (openForward) {
            if (util.isBlank(forward)) { m->mothurOut("[WARNING]: " + forward + " is blank, skipping.\n"); skip=true; }
        }else { m->mothurOut("[WARNING]: can't find " + forward + ", ignoring pair.\n"); }
        
        bool openReverse = true;
        if (reverse != ""){
            openReverse = util.checkLocations(reverse, current->getLocations());
            if (openReverse) {
                if (util.isBlank(reverse)) { m->mothurOut("[WARNING]: " + reverse + " is blank, skipping.\n"); skip=true; }
            }else { m->mothurOut("[WARNING]: can't find " + reverse + ", ignoring pair.\n"); }
        }
        
        bool openFindex = true;
        if ((findex != "") && (findex != "NONE")){
            openFindex = util.checkLocations(findex, current->getLocations());
            if (openFindex) {
                if (util.isBlank(findex)) { m->mothurOut("[WARNING]: " + findex + " is blank, skipping.\n"); skip=true; }
            }else { m->mothurOut("[WARNING]: can't find " + findex + ", ignoring pair.\n"); }
        }
        
        bool openRindex = true;
        if ((rindex != "") && (rindex != "NONE")) {
            openRindex = util.checkLocations(rindex, current->getLocations());
            if (openRindex) {
                if (util.isBlank(rindex)) { m->mothurOut("[WARNING]: " + rindex + " is blank, skipping.\n"); skip=true; }
            }else { m->mothurOut("[WARNING]: can't find " + rindex + ", ignoring pair.\n"); }
        }

        if ((openForward) && (openReverse) && (openFindex) && (openRindex) && (!skip)) { //good pair
            return false;
        }else { return true; }
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "validateFiles");
        exit(1);
    }
}
/**************************************************************************************************/
void FileFile::setGZ(string forward, string reverse, string findex, string rindex, bool& allGZ, bool& allPlainTxt){
    try {
        
#ifdef USE_BOOST
        if (util.isGZ(forward)[1]) { allPlainTxt = false;  }
        else {   allGZ = false;  }
        if (util.isGZ(reverse)[1]) { allPlainTxt = false;  }
        else {   allGZ = false;  }
        if ((findex != "") && (findex != "NONE")) {
            if (util.isGZ(findex)[1]) { allPlainTxt = false;  }
            else {   allGZ = false;  }
        }
        if ((rindex != "") && (rindex != "NONE")) {
            if (util.isGZ(rindex)[1]) { allPlainTxt = false;  }
            else {   allGZ = false;  }
        }
        if (!allGZ && !allPlainTxt) { //mixed bag of files, uh oh...
            m->mothurOut("[ERROR]: Your files must all be in compressed .gz form or all in plain text form.  Please correct. \n"); m->setControl_pressed(true);
        }
#else
        allGZ=false;
#if defined NON_WINDOWS
#else
        string extension = util.getExtension(forward);
        if (extension == "gz") {  m->mothurOut("[ERROR]: You cannot use compressed .gz files as input with our windows version of mothur. \n"); m->setControl_pressed(true); }
        extension = util.getExtension(reverse);
        if (extension == "gz") {  m->mothurOut("[ERROR]: You cannot use compressed .gz files as input with our windows version of mothur. \n"); m->setControl_pressed(true); }
        if ((findex != "") && (findex != "NONE")) {
            extension = util.getExtension(findex);
            if (extension == "gz") {  m->mothurOut("[ERROR]: You cannot use compressed .gz files as input with our windows version of mothur. \n"); m->setControl_pressed(true); }
        }
        if ((rindex != "") && (rindex != "NONE")) {
            extension = util.getExtension(rindex);
            if (extension == "gz") {  m->mothurOut("[ERROR]: You cannot use compressed .gz files as input with our windows version of mothur. \n"); m->setControl_pressed(true); }
        }
#endif
#endif
        
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "setGZ");
        exit(1);
    }
}
/**************************************************************************************************/


