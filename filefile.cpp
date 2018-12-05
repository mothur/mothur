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
        
        gz = false;
        hasIndex = false;
        createFileGroup = false;
        read(f, mode);
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "FileFile");
        exit(1);
    }
}
/**************************************************************************************************/

FileFile::FileFile(){
    try {
        m = MothurOut::getInstance();
        
        current = CurrentFile::getInstance();
        inputDir = current->getInputDir();
        mpath = current->getProgramPath();
        
        gz = false;
        hasIndex = false;
        createFileGroup = false;
        filename = ""; mode = "";
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "FileFile");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<string> > FileFile::read(string f, string mode){
    try {
        filename = f;
        
        if (mode == "contigs") {
            string forward, reverse, findex, rindex;
            bool allGZ = true;
            bool allPlainTxt = true;
            
            ifstream in;
            util.openInputFile(filename, in);
            
            while(!in.eof()) {
                
                if (m->getControl_pressed()) { return files; }
                
                bool skip = false;
                string line = util.getline(in);  util.gobble(in);
                
                if (m->getDebug()) {  m->mothurOut("[DEBUG]: " + line +"\n");  }
                
                vector<string> pieces = util.splitWhiteSpace(line);
                
                string group = "";
                if (pieces.size() == 2) {
                    forward = pieces[0];
                    reverse = pieces[1];
                    group = "";
                    findex = "";
                    rindex = "";
                }else if (pieces.size() == 3) {
                    group = pieces[0];
                    util.checkGroupName(group);
                    forward = pieces[1];
                    reverse = pieces[2];
                    findex = "";
                    rindex = "";
                    createFileGroup = true;
                }else if (pieces.size() == 4) {
                    forward = pieces[0];
                    reverse = pieces[1];
                    findex = pieces[2];
                    rindex = pieces[3];
                    if ((findex == "none") || (findex == "NONE")){ findex = "NONE"; }
                    if ((rindex == "none") || (rindex == "NONE")){ rindex = "NONE"; }
                }else {
                    m->mothurOut("[ERROR]: file lines can be 2, 3, or 4 columns. The forward fastq files in the first column and their matching reverse fastq files in the second column, or a groupName then forward fastq file and reverse fastq file, or forward fastq file then reverse fastq then forward index and reverse index file.  If you only have one index file add 'none' for the other one. \n"); m->setControl_pressed(true);
                }
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: group = " + group + ", forward = " + forward + ", reverse = " + reverse + ", forwardIndex = " + findex + ", reverseIndex = " + rindex + ".\n"); }
                
                if (inputDir != "") {
                    string path = util.hasPath(forward);
                    if (path == "") {  forward = inputDir + forward;  }
                    
                    path = util.hasPath(reverse);
                    if (path == "") {  reverse = inputDir + reverse;  }
                    
                    if (findex != "") {
                        path = util.hasPath(findex);
                        if (path == "") {  findex = inputDir + findex;  }
                    }
                    
                    if (rindex != "") {
                        path = util.hasPath(rindex);
                        if (path == "") {  rindex = inputDir + rindex;  }
                    }
                }
                
                //look for mothur exe
                
                
                //check to make sure both are able to be opened
                ifstream in2;
                bool openForward = util.openInputFile(forward, in2, "noerror");
                
                string tryPath = forward;
                //if you can't open it, try default location
                if (!openForward) {
                    if (current->getDefaultPath() != "") { //default path is set
                        tryPath = current->getDefaultPath() + util.getSimpleName(forward);
                        m->mothurOut("Unable to open " + forward + ". Trying default " + tryPath); m->mothurOutEndLine();
                        ifstream in3;
                        openForward = util.openInputFile(tryPath, in3, "noerror");
                        in3.close();
                        forward = tryPath;
                    }
                }
                
                //if you can't open it, try output location
                if (!openForward) {
                    if (current->getOutputDir() != "") { //default path is set
                        tryPath = current->getOutputDir() + util.getSimpleName(forward);
                        m->mothurOut("Unable to open " + forward + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                        ifstream in4;
                        openForward = util.openInputFile(tryPath, in4, "noerror");
                        forward = tryPath;
                        in4.close();
                    }
                }
                
                //if you can't open it, try mothur's location
                if (!openForward) {
                    tryPath = mpath + util.getSimpleName(forward);
                    m->mothurOut("Unable to open " + forward + ". Trying mothur's executable directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openForward = util.openInputFile(tryPath, in4, "noerror");
                    forward = tryPath;
                    in4.close();
                }
                
                if (!openForward) { //can't find it
                    m->mothurOut("[WARNING]: can't find " + forward + ", ignoring pair.\n");
                }else{
                    if (util.isBlank(tryPath)) { m->mothurOut("[WARNING]: " + forward + " is blank, skipping.\n"); skip=true; }
                    in2.close();
                }
                
                ifstream in3;
                bool openReverse = util.openInputFile(reverse, in3, "noerror");
                tryPath = reverse;
                
                //if you can't open it, try default location
                if (!openReverse) {
                    if (current->getDefaultPath() != "") { //default path is set
                        tryPath = current->getDefaultPath() + util.getSimpleName(reverse);
                        m->mothurOut("Unable to open " + reverse + ". Trying default " + tryPath); m->mothurOutEndLine();
                        ifstream in3;
                        openReverse = util.openInputFile(tryPath, in3, "noerror");
                        in3.close();
                        reverse = tryPath;
                    }
                }
                
                //if you can't open it, try mothur's location
                if (!openReverse) {
                    tryPath = mpath + util.getSimpleName(reverse);
                    m->mothurOut("Unable to open " + reverse + ". Trying mothur's executable directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openReverse = util.openInputFile(tryPath, in4, "noerror");
                    reverse = tryPath;
                    in4.close();
                }
                
                //if you can't open it, try output location
                if (!openReverse) {
                    if (current->getOutputDir() != "") { //default path is set
                        tryPath = current->getOutputDir() + util.getSimpleName(reverse);
                        m->mothurOut("Unable to open " + reverse + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                        ifstream in4;
                        openReverse = util.openInputFile(tryPath, in4, "noerror");
                        reverse = tryPath;
                        in4.close();
                    }
                }
                
                if (!openReverse) { //can't find it
                    m->mothurOut("[WARNING]: can't find " + reverse + ", ignoring pair.\n");
                }else{  if (util.isBlank(tryPath)) { m->mothurOut("[WARNING]: " + reverse + " is blank, skipping.\n"); skip=true; } in3.close();  }
                
                bool openFindex = true;
                if ((findex != "") && (findex != "NONE")){
                    ifstream in4;
                    openFindex = util.openInputFile(findex, in4, "noerror"); in4.close();
                    tryPath = findex;
                    
                    //if you can't open it, try default location
                    if (!openFindex) {
                        if (current->getDefaultPath() != "") { //default path is set
                            tryPath = current->getDefaultPath() + util.getSimpleName(findex);
                            m->mothurOut("Unable to open " + findex + ". Trying default " + tryPath); m->mothurOutEndLine();
                            ifstream in5;
                            openFindex = util.openInputFile(tryPath, in5, "noerror");
                            in5.close();
                            findex = tryPath;
                        }
                    }
                    
                    //if you can't open it, try mothur's location
                    if (!openFindex) {
                        tryPath = mpath + util.getSimpleName(findex);
                        m->mothurOut("Unable to open " + findex + ". Trying mothur's executable directory " + tryPath); m->mothurOutEndLine();
                        ifstream in14;
                        openFindex = util.openInputFile(tryPath, in14, "noerror");
                        findex = tryPath;
                        in14.close();
                    }
                    
                    //if you can't open it, try output location
                    if (!openFindex) {
                        if (current->getOutputDir() != "") { //default path is set
                            tryPath = current->getOutputDir() + util.getSimpleName(findex);
                            m->mothurOut("Unable to open " + findex + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                            ifstream in6;
                            openFindex = util.openInputFile(tryPath, in6, "noerror");
                            findex = tryPath;
                            in6.close();
                        }
                    }
                    
                    if (!openFindex) { //can't find it
                        m->mothurOut("[WARNING]: can't find " + findex + ", ignoring pair.\n");
                    }else{
                        if (util.isBlank(tryPath)) { m->mothurOut("[WARNING]: " + findex + " is blank, skipping.\n"); skip=true; }
                    }
                }
                
                bool openRindex = true;
                if ((rindex != "") && (rindex != "NONE")) {
                    ifstream in7;
                    openRindex = util.openInputFile(rindex, in7, "noerror"); in7.close();
                    tryPath = rindex;
                    
                    //if you can't open it, try default location
                    if (!openRindex) {
                        if (current->getDefaultPath() != "") { //default path is set
                            tryPath = current->getDefaultPath() + util.getSimpleName(rindex);
                            m->mothurOut("Unable to open " + rindex + ". Trying default " + tryPath); m->mothurOutEndLine();
                            ifstream in8;
                            openRindex = util.openInputFile(tryPath, in8, "noerror");
                            in8.close();
                            rindex = tryPath;
                        }
                    }
                    
                    //if you can't open it, try mothur's location
                    if (!openRindex) {
                        tryPath = mpath + util.getSimpleName(rindex);
                        m->mothurOut("Unable to open " + rindex + ". Trying mothur's executable directory " + tryPath); m->mothurOutEndLine();
                        ifstream in14;
                        openRindex = util.openInputFile(tryPath, in14, "noerror");
                        rindex = tryPath;
                        in14.close();
                    }
                    
                    //if you can't open it, try output location
                    if (!openRindex) {
                        if (current->getOutputDir() != "") { //default path is set
                            tryPath = current->getOutputDir() + util.getSimpleName(rindex);
                            m->mothurOut("Unable to open " + rindex + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                            ifstream in9;
                            openRindex = util.openInputFile(tryPath, in9, "noerror");
                            rindex = tryPath;
                            in9.close();
                        }
                    }
                    
                    if (!openRindex) { //can't find it
                        m->mothurOut("[WARNING]: can't find " + rindex + ", ignoring pair.\n");
                    }else{
                        if (util.isBlank(tryPath)) { m->mothurOut("[WARNING]: " + rindex + " is blank, skipping.\n"); skip=true; }
                    }
                }
                
                
                
                if ((openForward) && (openReverse) && (openFindex) && (openRindex) && (!skip)) { //good pair
                    file2Group[files.size()] = group;
                    vector<string> pair;
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
                    pair.push_back(forward);
                    pair.push_back(reverse);
                    pair.push_back(findex);
                    pair.push_back(rindex);
                    if (((findex != "") || (rindex != ""))) { hasIndex = true; }
                    files.push_back(pair);
                }
            }
            in.close();
            
            if (allGZ) {
                gz = true;
            }else { gz = false; }
        }
        if (files.size() == 0) { m->setControl_pressed(true); }
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "FileFile", "read");
        exit(1);
    }
}
/**************************************************************************************************/


