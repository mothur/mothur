//
//  utils.cpp
//  Mothur
//
//  Created by Sarah Westcott on 11/13/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "utils.hpp"
#include "ordervector.hpp"
#include "sharedordervector.h"

/***********************************************************************/
string getLabelTag(string label){
    
    string tag = "";
    
    //remove OTU or phylo tag
    string newLabel1 = "";
    for (int i = 0; i < label.length(); i++) {
        if(label[i]>47 && label[i]<58) { //is a digit
        }else {  tag += label[i];  }
    }
    
    return tag;
}
/***********************************************************************/
Utils::Utils(){
    try {

        m = MothurOut::getInstance();  modifyNames = m->getChangedSeqNames();
        long long s = m->getRandomSeed();
        mersenne_twister_engine.seed(s); srand(s);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }
}
/***********************************************************************/
float Utils::randomUniform() {
    try {
        uniform_real_distribution<float> unif;
        return (unif(mersenne_twister_engine));
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "randomUniform");
        exit(1);
    }
}
/***********************************************************************/
float Utils::randomExp() {
    try {
        exponential_distribution<float> unif;
        return (unif(mersenne_twister_engine));
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "randomExp");
        exit(1);
    }
}
/***********************************************************************/
float Utils::randomNorm() {
    try {
        normal_distribution<float> unif;
        return (unif(mersenne_twister_engine));
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "randomNorm");
        exit(1);
    }
}
/***********************************************************************/
float Utils::randomGamma(float range) {
    try {
        gamma_distribution<float> unif(range, range);
        return (unif(mersenne_twister_engine));
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "randomGamma");
        exit(1);
    }
}
/***********************************************************************/
vector<float> Utils::randomDirichlet(vector<float> alphas) {
    try {
        int nAlphas = (int)alphas.size();
        vector<float> dirs(nAlphas, 0);

        float sum = 0.0000;
        for(int i=0;i<nAlphas;i++){
            dirs[i] = randomGamma(alphas[i]);
						while(isinf(dirs[i])) { dirs[i] = randomGamma(alphas[i]);	}
            sum += dirs[i];
        }

        for(int i=0;i<nAlphas;i++){ dirs[i] /= sum; }

        return dirs;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "randomDirichlet");
        exit(1);
    }
}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<int>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<weightedSeq>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }
    
}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<long long>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }
    
}
/***********************************************************************/
void Utils::mothurRandomShuffle(OrderVector& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<SharedRAbundVector*>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(SharedOrderVector& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<string>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<item>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<PCell*>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<colDist>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }
    
}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector<PDistCellMin>& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
void Utils::mothurRandomShuffle(vector< vector<double> >& randomize){
    try {
        shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRandomShuffle");
        exit(1);
    }

}
/***********************************************************************/
long long Utils::getRandomIndex(long long highest){
    try {
        
        if (highest == 0) {  return 0; }
        
        uniform_int_distribution<long long> dis(0, highest);
        
        long long random = dis(mersenne_twister_engine);
        
        return random;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getRandomIndex");
        exit(1);
    }
}
/***********************************************************************/

int Utils::getRandomIndex(int highest){
    try {
        if (highest == 0) { return 0; }

        uniform_int_distribution<int> dis(0, highest);
        int random = dis(mersenne_twister_engine);
        return random;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getRandomIndex");
        exit(1);
    }

}
/***********************************************************************/
int Utils::getRandomNumber(){
    try {
        uniform_int_distribution<int> dis;

        int random = dis(mersenne_twister_engine);

        return random;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getRandomNumber");
        exit(1);
    }

}
/***********************************************************************/
double Utils::getRandomDouble0to1(){
    try {
        uniform_real_distribution<double> dis(0, 1);

        double random = dis(mersenne_twister_engine);

        return random;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getRandomNumber");
        exit(1);
    }

}

/*********************************************************************************************/
unsigned long long Utils::getRAMUsed() {
    try {

#if defined (__APPLE__) || (__MACH__)
        /* Mac: ru_maxrss gives the size in bytes */
        struct rusage r_usage;
        getrusage(RUSAGE_SELF, & r_usage);

        return r_usage.ru_maxrss;
#elif (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        /* Linux: ru_maxrss gives the size in kilobytes  */
        struct rusage r_usage;
        getrusage(RUSAGE_SELF, & r_usage);
        return r_usage.ru_maxrss * 1024;
#else
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        return (size_t)(status.ullTotalPhys - status.ullAvailPhys);
#endif
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getMemoryUsed");
        exit(1);
    }
}
/*********************************************************************************************/
unsigned long long Utils::getTotalRAM() {
    try {

#if defined NON_WINDOWS
#if defined _SC_PHYS_PAGES && defined _SC_PAGESIZE
        /* This works on linux-gnu, solaris2 and cygwin.  */
        double pages = sysconf (_SC_PHYS_PAGES);
        double pagesize = sysconf (_SC_PAGESIZE);
        if (0 <= pages && 0 <= pagesize)
            return pages * pagesize;
#else
        m->mothurOut("[WARNING]: Cannot determine amount of RAM");
#endif

#elif defined (_WIN32)
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        return (size_t)status.ullTotalPhys;
#else
        struct sysinfo si;
        if (sysinfo(&si))
            mothurOut("[WARNING]: Cannot determine amount of RAM");
        return si.totalram * si.mem_unit;

#endif
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getTotalRAM");
        exit(1);
    }
}
//********************************************************************/
string Utils::reverseOligo(string oligo){
    try {
        string reverse = "";

        for(int i=oligo.length()-1;i>=0;i--){

            if(oligo[i] == 'A')		{	reverse += 'T';	}
            else if(oligo[i] == 'T'){	reverse += 'A';	}
            else if(oligo[i] == 'U'){	reverse += 'A';	}

            else if(oligo[i] == 'G'){	reverse += 'C';	}
            else if(oligo[i] == 'C'){	reverse += 'G';	}

            else if(oligo[i] == 'R'){	reverse += 'Y';	}
            else if(oligo[i] == 'Y'){	reverse += 'R';	}

            else if(oligo[i] == 'M'){	reverse += 'K';	}
            else if(oligo[i] == 'K'){	reverse += 'M';	}

            else if(oligo[i] == 'W'){	reverse += 'W';	}
            else if(oligo[i] == 'S'){	reverse += 'S';	}

            else if(oligo[i] == 'B'){	reverse += 'V';	}
            else if(oligo[i] == 'V'){	reverse += 'B';	}

            else if(oligo[i] == 'D'){	reverse += 'H';	}
            else if(oligo[i] == 'H'){	reverse += 'D';	}

            else						{	reverse += 'N';	}
        }


        return reverse;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "reverseOligo");
        exit(1);
    }
}

/*********************************************************************************************/
bool Utils::fileExists(string name)  {
    try {
        bool exists = false;
        name = getFullPathName(name);

    #if defined NON_WINDOWS
        ifstream in;
        openInputFile(name, in, "");

        //if this file exists
        if (in) { in.close(); exists = true;  }
    #else

        DWORD attributes = GetFileAttributes(name.c_str());
        exists = (attributes != INVALID_FILE_ATTRIBUTES && !(attributes & FILE_ATTRIBUTE_DIRECTORY));
    #endif

        return exists;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "fileExists");
        exit(1);
    }
}
/***********************************************************************/
string Utils::getFullPathName(string fileName){
    try{
        string path = hasPath(fileName);
        string newFileName;
        int pos;
        vector<string> dirs;
        int index = 0;

        if (path == "") { return fileName; } //its a simple name
        else { //we need to complete the pathname
            // ex. ../../../filename
            // cwd = /user/work/desktop
            //get current working directory
            string cwd;
            if (path.find("~") != string::npos) { //go to home directory
#if defined NON_WINDOWS
                string homeDir;
                char *homepath = NULL; homepath = getenv ("HOME");
                if ( homepath != NULL) { homeDir = homepath; }
                else { homeDir = "";  }
#else
                string homeDir = getenv ("HOMEPATH");
#endif
                newFileName = homeDir + fileName.substr(fileName.find("~")+1);
                return newFileName;
            }else { //find path
                string pattern = ".";  pattern += PATH_SEPARATOR;
                if (path.rfind(pattern) == string::npos) { return fileName; } //already complete name
                else { newFileName = fileName.substr(fileName.rfind(pattern)+2); } //save the complete part of the name

                char *cwdpath = NULL; cwdpath = getcwd(NULL, 0); // or _getcwd
                if (cwdpath != NULL)    { cwd = cwdpath;    }
                else                    { cwd = "";         }

                //rip off first '/'
                string simpleCWD = cwd;
#if defined NON_WINDOWS
                if (cwd.length() > 0) { simpleCWD = cwd.substr(1); }
#endif
                //break apart the current working directory
                while (simpleCWD.find_first_of(PATH_SEPARATOR) != string::npos) {
                    string dir = simpleCWD.substr(0,simpleCWD.find_first_of(PATH_SEPARATOR));
                    simpleCWD = simpleCWD.substr(simpleCWD.find_first_of(PATH_SEPARATOR)+1, simpleCWD.length());
                    dirs.push_back(dir);
                }
                //get last one              // ex. ../../../filename = /user/work/desktop/filename
                dirs.push_back(simpleCWD);  //ex. dirs[0] = user, dirs[1] = work, dirs[2] = desktop

                index = dirs.size()-1;
                string searchString = "."; searchString += PATH_SEPARATOR;
                while((pos = path.rfind(searchString)) != string::npos) { //while you don't have a complete path
                    if (pos == 0) { break;  //you are at the end
                    }else if (path[(pos-1)] == '.') { //you want your parent directory ../
                        path = path.substr(0, pos-1);
                        index--;
                        if (index == 0) {  break; }
                    }else if (path[(pos-1)] == '/') { //you want the current working dir ./
                        path = path.substr(0, pos);
                    }else if (pos == 1) { break;  //you are at the end
                    }else { m->mothurOut("[ERROR}: Can not resolve path for " +  fileName + "\n"); m->setControl_pressed(true); return fileName;  }
                }
            }

            for (int i = index; i >= 0; i--) { newFileName = dirs[i] +  PATH_SEPARATOR + newFileName; }

#if defined NON_WINDOWS
            newFileName =  PATH_SEPARATOR +  newFileName;
#endif
            return newFileName;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getFullPathName");
        exit(1);
    }
}
/********************************************************************/
string Utils::findProgramPath(string programName){
    try {
        string pPath = "";

        //look in ./
        //is this the programs path?
        string tempIn = ".";
        tempIn += PATH_SEPARATOR;

        //if this file exists
        if (fileExists(tempIn+programName)) { pPath = getFullPathName(tempIn); if (m->getDebug()) { m->mothurOut("[DEBUG]: found it, programPath = " + pPath + "\n"); } return pPath;   }

        string envPath = getenv("PATH");

        //delimiting path char
        char delim;
#if defined NON_WINDOWS
        delim = ':';
#else
        delim = ';';
#endif

        //break apart path variable by ':'
        vector<string> dirs;
        splitAtChar(envPath, dirs, delim);

        if (m->getDebug()) { m->mothurOut("[DEBUG]: dir's in path: \n"); }

        //get path related to mothur
        pPath = "";
        for (int i = 0; i < dirs.size(); i++) {

            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + dirs[i] + "\n"); }

            //to lower so we can find it
            string tempLower = "";
            for (int j = 0; j < dirs[i].length(); j++) {  tempLower += tolower(dirs[i][j]);  }

            //is this mothurs path?
            if (tempLower.find(programName) != -1) {  pPath = dirs[i]; break;  }
        }

        if (m->getDebug()) { m->mothurOut("[DEBUG]: programPath = " + pPath + "\n"); }

        //add programName so it looks like what argv would look like
        if (pPath != "") { pPath += PATH_SEPARATOR;  }
        else {
            //okay programName is not in the path, so the folder programName is in must be in the path
            //lets find out which one

            //get path related to the program
            for (int i = 0; i < dirs.size(); i++) {

                if (m->getDebug()) { m->mothurOut("[DEBUG]: looking in " + dirs[i] + " for " + programName + " \n"); }

                //is this the programs path?
                string tempIn = dirs[i] + PATH_SEPARATOR;

                //if this file exists
                if (fileExists(tempIn + programName)) { pPath = getFullPathName(tempIn); if (m->getDebug()) { m->mothurOut("[DEBUG]: found it, programPath = " + pPath + "\n"); } break;   }
            }
        }

#if defined NON_WINDOWS
#else
        if (pPath == "") {
            char buffer[MAX_PATH];
            GetModuleFileName(NULL, buffer, MAX_PATH) ;

            pPath = buffer;
            pPath = getPathName(pPath);

            //if this file exists
            if (fileExists(pPath + programName)) { pPath = getFullPathName(pPath); if (m->getDebug()) { m->mothurOut("[DEBUG]: found it, programPath = " + pPath + "\n"); } }
        }
#endif
        return pPath;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "findProgramPath");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::checkLocations(string& filename, vector<string> locations){
    try {
        filename = getFullPathName(filename);
        string inputDir = locations[0];
        string outputDir = locations[1];
        string defaultPath = locations[2];
        string mothurPath = locations[3];

        bool ableToOpen;
        ifstream in;
        ableToOpen = openInputFile(filename, in, "noerror");
        in.close();

        //if you can't open it, try input location
        if (!ableToOpen) {
            if (inputDir != "") { //default path is set
                string tryPath = inputDir + getSimpleName(filename);
                m->mothurOut("Unable to open " + filename + ". Trying input directory " + tryPath + ".\n");
                ifstream in2; ableToOpen = openInputFile(tryPath, in2, "noerror"); in2.close();
                filename = tryPath;
            }
        }

        //if you can't open it, try output location
        if (!ableToOpen) {
            if (outputDir != "") { //default path is set
                string tryPath = outputDir + getSimpleName(filename);
                m->mothurOut("Unable to open " + filename + ". Trying output directory " + tryPath+ ".\n");
                ifstream in2; ableToOpen = openInputFile(tryPath, in2, "noerror"); in2.close();
                filename = tryPath;
            }
        }


        //if you can't open it, try default location
        if (!ableToOpen) {
            if (defaultPath != "") { //default path is set
                string tryPath = defaultPath + getSimpleName(filename);
                m->mothurOut("Unable to open " + filename + ". Trying default " + tryPath+ ".\n");
                ifstream in2; ableToOpen = openInputFile(tryPath, in2, "noerror"); in2.close();
                filename = tryPath;
            }
        }

        //if you can't open it its not in current working directory or inputDir, try mothur excutable location
        if (!ableToOpen) {
            string tryPath = mothurPath + getSimpleName(filename);
            m->mothurOut("Unable to open " + filename + ". Trying mothur's executable location " + tryPath+ ".\n");
            ifstream in2; ableToOpen = openInputFile(tryPath, in2, "noerror"); in2.close();
            filename = tryPath;
        }


        if (!ableToOpen) { m->mothurOut("Unable to open " + filename + ".\n");  return false;  }

        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "checkLocations");
        exit(1);
    }
}
/***********************************************************************/

bool Utils::openInputFile(string fileName, ifstream& fileHandle, string mode){
    try {
        //get full path name
        string completeFileName = getFullPathName(fileName);

        fileHandle.open(completeFileName.c_str());
        if(!fileHandle) { return false;  }
        else {
            //check for blank file
            zapGremlins(fileHandle);
            gobble(fileHandle);
            return true;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openInputFile - no Error");
        exit(1);
    }
}
/***********************************************************************/

bool Utils::openInputFile(string fileName, ifstream& fileHandle){
    try {

        //get full path name
        string completeFileName = getFullPathName(fileName);

        fileHandle.open(completeFileName.c_str());
        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + completeFileName + "\n"); return false; }
        else {
            //check for blank file
            zapGremlins(fileHandle);
            gobble(fileHandle);
            if (fileHandle.eof()) { m->mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct.\n");   }
            return true;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openInputFile");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::openInputFileBinary(string fileName, ifstream& fileHandle){
    try {

        //get full path name
        string completeFileName = getFullPathName(fileName);

        fileHandle.open(completeFileName.c_str(), ios::binary);
        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + completeFileName+ "\n");  return false; }
        else {
            //check for blank file
            zapGremlins(fileHandle);
            gobble(fileHandle);
            if (fileHandle.eof()) { m->mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct.\n");  }
            return true;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openInputFileBinary");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::openInputFileBinary(string fileName, ifstream& fileHandle, string noerror){
    try {

        //get full path name
        string completeFileName = getFullPathName(fileName);

        fileHandle.open(completeFileName.c_str(), ios::binary);
        if(!fileHandle) { return false; }
        else {
            //check for blank file
            zapGremlins(fileHandle);
            gobble(fileHandle);
            return true;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openInputFileBinary - no error");
        exit(1);
    }
}
/***********************************************************************/
#ifdef USE_BOOST
bool Utils::openInputFileBinary(string fileName, ifstream& file, boost::iostreams::filtering_istream& in){
    try {

        //get full path name
        string completeFileName = getFullPathName(fileName);

        file.open(completeFileName.c_str(), ios_base::in | ios_base::binary);

        if(!file) { m->mothurOut("[ERROR]: Could not open " + completeFileName + "\n"); return false; }
        else {
            //check for blank file
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
            if (file.eof()) { m->mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct.\n");  }
            return true;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openInputFileGZBinary");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::openInputFileBinary(string fileName, ifstream& file, boost::iostreams::filtering_istream& in, string noerror){
    try {

        //get full path name
        string completeFileName = getFullPathName(fileName);

        file.open(completeFileName.c_str(), ios_base::in | ios_base::binary);

        if(!file) { return false; }
        else { //check for blank file
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
            return true;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openInputFileGZBinary - no error");
        exit(1);
    }
}
#endif
/***********************************************************************/
//results[0] = allGZ, results[1] = allNotGZ
vector<bool> Utils::allGZFiles(vector<string> & files){
    try {
        vector<bool> results;
        bool allGZ = true;
        bool allNOTGZ = true;

        for (int i = 0; i < files.size(); i++) {
            if (m->getControl_pressed()) { break; }

            //ignore none and blank filenames
            if ((files[i] != "") || (files[i] != "NONE")) {
                if (isGZ(files[i])[1]) { allNOTGZ = false;  }
                else {  allGZ = false;  }
            }
        }

        if (!allGZ && !allNOTGZ) { //mixed bag
            m->mothurOut("[ERROR]: Cannot mix .gz and non compressed files. Please decompress your files and rerun.\n"); m->setControl_pressed(true);
        }

        results.push_back(allGZ);
        results.push_back(allNOTGZ);

        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "areGZFiles");
        exit(1);
    }
}

/***********************************************************************/
vector<bool> Utils::isGZ(string filename){
    try {
        vector<bool> results; results.resize(2, false);
#ifdef USE_BOOST
        ifstream fileHandle;
        boost::iostreams::filtering_istream gzin;

        if (getExtension(filename) != ".gz") { return results; } // results[0] = false; results[1] = false;

        bool ableToOpen = openInputFileBinary(filename, fileHandle, gzin, ""); //no error
        if (!ableToOpen) { return results; } // results[0] = false; results[1] = false;
        else {  results[0] = true;  }

        char c;
        try
        {
            gzin >> c;
            results[1] = true;
        }
        catch ( boost::iostreams::gzip_error & e )
        {
            gzin.pop();
            fileHandle.close();
            return results;  // results[0] = true; results[1] = false;
        }
        fileHandle.close();
#else
        m->mothurOut("[ERROR]: cannot test for gz format without enabling boost libraries.\n"); m->setControl_pressed(true);
#endif
        return results; //results[0] = true; results[1] = true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isGZ");
        exit(1);
    }
}

/***********************************************************************/

int Utils::renameFile(string oldName, string newName){
    try {
        if(m->getDebug()) { m->mothurOut("[DEBUG]: renaming " + oldName + " to " + newName + "\n"); }
        
        if (oldName == newName) { return 0; }

        ifstream inTest;
        bool exist = openInputFile(newName, inTest, "");
        inTest.close();

#if defined NON_WINDOWS
        if (exist) { //you could open it so you want to delete it
            if(m->getDebug()) { m->mothurOut("[DEBUG]: removing old copy of " + newName + "\n"); }
            string command = "rm " + newName;
            system(command.c_str());
        }
        
        if(m->getDebug()) { m->mothurOut("[DEBUG]: mv " + oldName + " to " + newName + "\n"); }
        
        string command = "mv " + oldName + " " + newName;
        
        if(m->getDebug()) { m->mothurOut("[DEBUG]: running system command mv " + oldName + " " + newName + "\n"); }
        
        int returnCode = system(command.c_str());
        
        if(m->getDebug()) { m->mothurOut("[DEBUG]: system command mv " + oldName + " " + newName + " returned " + toString(returnCode) + "\n"); }
        
        if (returnCode != 0) {
            int renameOk = rename(oldName.c_str(), newName.c_str());
        
            if(m->getDebug()) { m->mothurOut("[DEBUG]: rename " + oldName + " " + newName + " returned " + toString(renameOk) + "\n"); }
        }
#else
        mothurRemove(newName);
        int renameOk = rename(oldName.c_str(), newName.c_str());
        
        if(m->getDebug()) { m->mothurOut("[DEBUG]: rename " + oldName + " " + newName + " returned " + toString(renameOk) + "\n"); }
#endif
        return 0;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "renameFile");
        exit(1);
    }
}

/***********************************************************************/

bool Utils::openOutputFile(string fileName, ofstream& fileHandle){
    try {
        string completeFileName = getFullPathName(fileName);
        fileHandle.open(completeFileName.c_str(), ios::trunc);

        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + completeFileName + "\n"); return false; }
        else { return true; }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openOutputFile");
        exit(1);
    }

}
/***********************************************************************/

bool Utils::openOutputFileBinary(string fileName, ofstream& fileHandle){
    try {
        string completeFileName = getFullPathName(fileName);
        fileHandle.open(completeFileName.c_str(), ios::trunc | ios::binary);

        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + completeFileName + "\n");  return false; }
        else { return true; }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openOutputFileBinary");
        exit(1);
    }
}
/**************************************************************************************************/
int Utils::appendFiles(string temp, string filename) {
    try{
        ofstream output;
        ifstream input;

        //open output file in append mode
        openOutputFileBinaryAppend(filename, output);
        bool ableToOpen = openInputFileBinary(temp, input, "no error");
        //bool ableToOpen = openInputFile(temp, input);

        int numLines = 0;
        if (ableToOpen) { //you opened it
            char buffer[4096];
            while (!input.eof()) {
                input.read(buffer, 4096);
                output.write(buffer, input.gcount());
                //count number of lines
                for (int i = 0; i < input.gcount(); i++) {  if (buffer[i] == '\n') {numLines++;} }
            }
            input.close();
        }
        output.close();

        return numLines;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "appendFiles");
        exit(1);
    }
}
/**************************************************************************************************/
int Utils::appendFilesFront(string temp, string filename) {
    try{
        ofstream output;
        ifstream input;

        //open output file in append mode
        openOutputFileBinaryAppend(temp, output);
        bool ableToOpen = openInputFileBinary(filename, input, "no error");

        if (ableToOpen) { //you opened it
            char buffer[4096];
            while (!input.eof()) {
                input.read(buffer, 4096);
                output.write(buffer, input.gcount());
            }
            input.close();
        }
        output.close();

        mothurRemove(filename);
        renameFile(temp, filename);
        mothurRemove(temp);

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "appendFiles");
        exit(1);
    }
}

/**************************************************************************************************/
bool Utils::appendBinaryFiles(string temp, string filename) {
    try{
        ofstream output;
        ifstream input;

        //open output file in append mode
        openOutputFileBinaryAppend(filename, output);
        bool ableToOpen = openInputFileBinary(temp, input, "no error");

        if (ableToOpen) { //you opened it

            char buffer[4096];
            while (!input.eof()) {
                input.read(buffer, 4096);
                output.write(buffer, input.gcount());
            }
            input.close();
        }
        output.close();

        return ableToOpen;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "appendBinaryFiles");
        exit(1);
    }
}
/**************************************************************************************************/
bool Utils::appendSFFFiles(string temp, string filename) {
    try{
        ofstream output;
        ifstream input;
        bool ableToOpen = 0;

        //open output file in append mode
        string fullFileName = getFullPathName(filename);

        output.open(fullFileName.c_str(), ios::app | ios::binary);
        if(!output) { m->mothurOut("[ERROR]: Could not open " + fullFileName + "\n"); return false;  }
        else {
            //get full path name
            string completeFileName = getFullPathName(temp);
            input.open(completeFileName.c_str(), ios::binary);
            if(!input) { return false; }
            else {
                char buffer[4096];
                while (!input.eof()) {
                    input.read(buffer, 4096);
                    output.write(buffer, input.gcount());
                }
                input.close();
            }
            output.close();
        }

        return ableToOpen;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "appendSFFFiles");
        exit(1);
    }
}
/**************************************************************************************************/
int Utils::appendFilesWithoutHeaders(string temp, string filename) {
    try{
        ofstream output;
        ifstream input;

        //open output file in append mode
        openOutputFileAppend(filename, output);
        bool ableToOpen = openInputFile(temp, input, "no error");

        int numLines = 0;
        if (ableToOpen) { //you opened it

            string headers = getline(input); gobble(input);
            char buffer[4096];
            while (!input.eof()) {
                input.read(buffer, 4096);
                output.write(buffer, input.gcount());
                //count number of lines
                for (int i = 0; i < input.gcount(); i++) {  if (buffer[i] == '\n') {numLines++;} }
            }
            input.close();
        }

        output.close();

        return numLines;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "appendFiles");
        exit(1);
    }
}
/**************************************************************************************************/
string Utils::sortFile(string distFile, string outputDir){
    try {

        //if (outputDir == "") {  outputDir += hasPath(distFile);  }
        string outfile = getRootName(distFile) + "sorted.dist";


        //if you can, use the unix sort since its been optimized for years
#if defined NON_WINDOWS
        string command = "sort -n -k +3 " + distFile + " -o " + outfile;
        system(command.c_str());
#else //you are stuck with my best attempt...
        //windows sort does not have a way to specify a column, only a character in the line
        //since we cannot assume that the distance will always be at the the same character location on each line
        //due to variable sequence name lengths, I chose to force the distance into first position, then sort and then put it back.

        //read in file line by file and put distance first
        string tempDistFile = distFile + ".temp";
        ifstream input;
        ofstream output;
        openInputFile(distFile, input);
        openOutputFile(tempDistFile, output);

        string firstName, secondName;
        float dist;
        while (!input.eof()) {
            input >> firstName >> secondName >> dist;
            output << dist << '\t' << firstName << '\t' << secondName << endl;
            gobble(input);
        }
        input.close();
        output.close();


        //sort using windows sort
        string tempOutfile = outfile + ".temp";
        string command = "sort " + tempDistFile + " /O " + tempOutfile;
        system(command.c_str());

        //read in sorted file and put distance at end again
        ifstream input2;
        ofstream output2;
        openInputFile(tempOutfile, input2);
        openOutputFile(outfile, output2);

        while (!input2.eof()) {
            input2 >> dist >> firstName >> secondName;
            output2 << firstName << '\t' << secondName << '\t' << dist << endl;
            gobble(input2);
        }
        input2.close();
        output2.close();

        //remove temp files
        mothurRemove(tempDistFile);
        mothurRemove(tempOutfile);
#endif

        return outfile;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "sortFile");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::openOutputFileAppend(string fileName, ofstream& fileHandle){
    try {
        fileName = getFullPathName(fileName);

        fileHandle.open(fileName.c_str(), ios::app);
        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + fileName + "\n"); return 1; }else { return false; }
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openOutputFileAppend");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::openOutputFileBinaryAppend(string fileName, ofstream& fileHandle){
    try {
        fileName = getFullPathName(fileName);

        fileHandle.open(fileName.c_str(), ios::app | ios::binary);
        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + fileName + "\n"); return 1; }else { return false; }
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "openOutputFileAppend");
        exit(1);
    }
}

/***********************************************************************/
void Utils::gobble(istream& f){
    try {

        char d;
        while(isspace(d=f.get()))		{ ;}
        if(!f.eof()) { f.putback(d); }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "gobble");
        exit(1);
    }
}
/***********************************************************************/
void Utils::gobble(istringstream& f){
    try {
        char d;
        while(isspace(d=f.get()))		{;}
        if(!f.eof()) { f.putback(d); }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "gobble");
        exit(1);
    }
}
/***********************************************************************/
void Utils::zapGremlins(istream& f){
    try {

        char d;
        while('\0'==(d=f.get()))		{ ;}
        if(!f.eof()) { f.putback(d); }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "zapGremlins");
        exit(1);
    }
}
/***********************************************************************/
void Utils::zapGremlins(istringstream& f){
    try {
        char d;
        while('\0'==(d=f.get()))		{ ;}
        if(!f.eof()) { f.putback(d); }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "zapGremlins");
        exit(1);
    }
}

/***********************************************************************/
string Utils::getline(istringstream& fileHandle) {
    try {
        string line = "";
        while (!fileHandle.eof())	{
            //get next character
            char c = fileHandle.get();

            //are you at the end of the line
            if ((c == '\n') || (c == '\r') || (c == '\f')){  break;	}
            else {		line += c;		}
        }

        return line;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getline");
        exit(1);
    }
}
/***********************************************************************/
string Utils::getline(ifstream& fileHandle) {
    try {
        string line = "";
        while (fileHandle)	{
            //get next character
            char c = fileHandle.get();

            //are you at the end of the line
            if ((c == '\n') || (c == '\r') || (c == '\f') || (c == EOF)){  break;	}
            else {		line += c;		}
        }

        return line;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getline");
        exit(1);
    }
}
#ifdef USE_BOOST
/***********************************************************************/
string Utils::getline(boost::iostreams::filtering_istream& fileHandle) {
    try {
        string line = "";
        while (fileHandle)	{
            //get next character
            char c = fileHandle.get();
            
            //are you at the end of the line
            if ((c == '\n') || (c == '\r') || (c == '\f') || (c == EOF)){ break; }
            else {		line += c;		}
        }
        
        return line;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getline");
        exit(1);
    }
}
#endif
/**********************************************************************/
string Utils::getPathName(string longName){
    try {
        string rootPathName = longName;

        if(longName.find_last_of("/\\") != longName.npos){
            int pos = longName.find_last_of("/\\")+1;
            rootPathName = longName.substr(0, pos);
        }

        return rootPathName;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getPathName");
        exit(1);
    }
}
/***********************************************************************/
string Utils::getRootName(string longName){
    try {

        string rootName = longName;

        if(rootName.find_last_of(".") != rootName.npos){
            int pos = rootName.find_last_of('.')+1;
            rootName = rootName.substr(0, pos);
        }

        return rootName;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getRootName");
        exit(1);
    }
}
/***********************************************************************/

string Utils::getSimpleName(string longName){
    try {
        string simpleName = longName;

        size_t found; found=longName.find_last_of("/\\");

        if(found != longName.npos){ simpleName = longName.substr(found+1); }

        return simpleName;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getSimpleName");
        exit(1);
    }
}
//**********************************************************************************************************************
string Utils::getStringFromVector(vector<string>& list, string delim){
    try {
        string result = "";

        if (list.size() == 0) { return result; }

        result = list[0];

        for (int i = 1; i < list.size(); i++) {
            if (m->getControl_pressed()) { break;  }
            result += delim + list[i];
        }

        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getStringFromVector");
        exit(1);
    }
}
//**********************************************************************************************************************
string Utils::getStringFromVector(vector<int>& list, string delim){
    try {
        string result = "";

        if (list.size() == 0) { return result; }

        result = toString(list[0]);

        for (int i = 1; i < list.size(); i++) {
            if (m->getControl_pressed()) { break;  }
            string temp = toString(list[i]);
            result += delim + temp;
        }

        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getStringFromVector");
        exit(1);
    }
}
//**********************************************************************************************************************
string Utils::getStringFromVector(vector<double>& list, string delim){
    try {
        string result = "";

        if (list.size() == 0) { return result; }

        result = toString(list[0]);

        for (int i = 1; i < list.size(); i++) {
            if (m->getControl_pressed()) { break;  }
            string temp = toString(list[i]);
            result += delim + temp;
        }

        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getStringFromVector");
        exit(1);
    }
}
//**********************************************************************************************************************
string Utils::getStringFromSet(set<int>& list, string delim){
    try {
        string result = "";
        
        if (list.size() == 0) { return result; }
        
        vector<int> vlist;
        for (set<int>::iterator it = list.begin(); it != list.end(); it++) {
            if (m->getControl_pressed()) { break;  }
            int value = *it;
            vlist.push_back(value);
        }
        result = getStringFromVector(vlist, delim);
        
        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getStringFromVector");
        exit(1);
    }
}
//**********************************************************************************************************************
//NOTE: assumes questions.size() == qanswers.size(), issues.size() == ianswers.size(), howtos.size() == hanswers.size()
string Utils::getFormattedHelp(vector<string> questions, vector<string> qanswers, vector<string> issues, vector<string> ianswers, vector<string> howtos,vector<string> hanswers) {
    try {
        
        string commonQuestions = ""; vector<string> headers;
        string header = "\nCommon Questions: \n"; headers.push_back(header);
        header = "\nCommon Issues: \n"; headers.push_back(header);
        header = "\nHow To: \n"; headers.push_back(header);
        
        commonQuestions += headers[0]+"\n";
#if defined NON_WINDOWS
        cout << BOLDGREEN << headers[0]; cout << RESET << endl;
#endif
        
        for (int i = 0; i < questions.size(); i++) {
            commonQuestions += toString(i+1) + ". " + questions[i]+"\n"+qanswers[i]+"\n";
#if defined NON_WINDOWS
            cout << BOLDBLUE << toString(i+1)+". "+questions[i]; cout << RESET << endl << qanswers[i] << endl;
#endif
        }
        
        if (questions.size() == 0) {
            commonQuestions += "Can't find your question? Please feel free to ask questions on our forum, https://forum.mothur.org.\n\n";
#if defined NON_WINDOWS
            cout << RESET "Can't find your question? Please feel free to ask questions on our forum, https://forum.mothur.org.\n\n";
#endif

        }
        
        commonQuestions += headers[1]+"\n";
#if defined NON_WINDOWS
        cout << BOLDGREEN << headers[1]; cout << RESET << endl;
#endif
        
        for (int i = 0; i < issues.size(); i++) {
            commonQuestions += toString(i+1)+". "+issues[i]+"\n"+ianswers[i]+"\n";
#if defined NON_WINDOWS
            cout << BOLDBLUE << toString(i+1)+". "+issues[i]; cout << RESET << endl << ianswers[i] << endl;
#endif
        }
        
        if (issues.size() == 0) {
            commonQuestions += "Can't find your issue? Please feel free to ask questions on our forum, https://forum.mothur.org or send bug reports to mothur.bugs@gmail.com.\n\n";
#if defined NON_WINDOWS
            cout << RESET "Can't find your issue? Please feel free to ask questions on our forum, https://forum.mothur.org or send bug reports to mothur.bugs@gmail.com.\n\n";
#endif
            
        }


        commonQuestions += headers[2]+"\n";
#if defined NON_WINDOWS
        cout << BOLDGREEN << headers[2]; cout << RESET << endl;
#endif
        
        for (int i = 0; i < howtos.size(); i++) {
            commonQuestions += toString(i+1) + ". " + howtos[i]+"\n"+hanswers[i]+"\n";
#if defined NON_WINDOWS
            cout << BOLDBLUE << toString(i+1)+". "+howtos[i]; cout << RESET << endl << hanswers[i] << endl;
#endif
        }
        
        if (howtos.size() == 0) {
            commonQuestions += "Not sure how to do what you want? Please feel free to ask questions on our forum, https://forum.mothur.org.\n\n";
#if defined NON_WINDOWS
            cout << RESET "Not sure how to do what you want? Please feel free to ask questions on our forum, https://forum.mothur.org.\n\n";
#endif
            
        }
        
#if defined NON_WINDOWS
        m->mothurOutJustToLog(commonQuestions);
        
        cout << BOLDMAGENTA << "\nFor further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki.\n"; cout << RESET << endl;
        m->mothurOutJustToLog("\nFor further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki.\n");
#else
        m->mothurOut(commonQuestions + "\nFor further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki.\n");
#endif

        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getFormattedHelp");
        exit(1);
    }
}
//**********************************************************************************************************************
string Utils::removeNs(string seq){
    try {
        string newSeq = "";
        for (int i = 0; i < seq.length(); i++) { if (seq[i] != 'N') {  newSeq += seq[i]; } }
        return newSeq;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "removeNs");
        exit(1);
    }
}
/***********************************************************************/
int Utils::getOTUNames(vector<string>& currentLabels, int numBins, string tagHeader){
    try {

        if (currentLabels.size() == numBins) {  return 0; }

        int maxLabelNumber = 0;
        if (currentLabels.size() < numBins) {
            string snumBins = toString(numBins);

            for (int i = 0; i < numBins; i++) {
                string binLabel = tagHeader;
                if (i < currentLabels.size()) { //label exists
                    if (getLabelTag(currentLabels[i]) == tagHeader) { //adjust 0's??
                        string sbinNumber = getSimpleLabel(currentLabels[i]);
                        int tempBinNumber; mothurConvert(sbinNumber, tempBinNumber);
                        if (tempBinNumber > maxLabelNumber) { maxLabelNumber = tempBinNumber; }
                        if (sbinNumber.length() < snumBins.length()) {
                            int diff = snumBins.length() - sbinNumber.length();
                            for (int h = 0; h < diff; h++) { binLabel += "0"; }
                        }
                        binLabel += sbinNumber;
                        currentLabels[i] = binLabel;
                    }
                }else{ //create new label
                    string sbinNumber = toString(maxLabelNumber+1); maxLabelNumber++;
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    currentLabels.push_back(binLabel);
                }
            }
        }
        return currentLabels.size();

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getOTUNames");
        exit(1);
    }
}
/**************************************************************************************/
void Utils::getCombos(vector<string>& groupComb, vector<string> userGroups, int& numComp) { //groupcomb, Groups, numcomb
    try {
        sort(userGroups.begin(), userGroups.end());

        //calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
        numComp = 0;
        for (int i=0; i< userGroups.size(); i++) {
            numComp += i;
            for (int l = 0; l < i; l++) {  //set group comparison labels
                if (userGroups[i] > userGroups[l])  { groupComb.push_back(userGroups[l] + "-" + userGroups[i]);     }
                else                                { groupComb.push_back(userGroups[i] + "-" + userGroups[l]);     }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getCombos");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::dirCheck(string& dirName){
    try {

        if (dirName == "") { return false; }

        //add / to name if needed
        string lastChar = dirName.substr(dirName.length()-1);
        if (lastChar != PATH_SEPARATOR) { dirName += PATH_SEPARATOR; }

        //test to make sure directory exists
        dirName = getFullPathName(dirName);
        string outTemp = dirName + "temp"+ toString(time(NULL));
        ofstream out;
        out.open(outTemp.c_str(), ios::trunc);
        if(!out) { m->mothurOut(dirName + " directory does not exist or is not writable.\n");  }
        else{ out.close(); mothurRemove(outTemp); return true; }

        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "dirCheck");
        exit(1);
    }

}
/***********************************************************************/

bool Utils::dirCheck(string& dirName, string noError){
    try {

        if (dirName == "") { return false; }

        //add / to name if needed
        string lastChar = dirName.substr(dirName.length()-1);
        if (lastChar != PATH_SEPARATOR) { dirName += PATH_SEPARATOR; }

        //test to make sure directory exists
        dirName = getFullPathName(dirName);
        string outTemp = dirName  + "temp"+ toString(time(NULL));
        ofstream out;
        out.open(outTemp.c_str(), ios::trunc);
        if(out) { out.close();  mothurRemove(outTemp); return true; }

        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "dirCheck - noError");
        exit(1);
    }

}
/***********************************************************************/
//returns true it exits or if we can make it
bool Utils::mkDir(string& dirName){
    try {
        bool dirExist = dirCheck(dirName, "noError");
        if (dirExist) { return true; }

        string makeDirectoryCommand = "mkdir -p \"" + dirName + "\"";
        system(makeDirectoryCommand.c_str());
        if (dirCheck(dirName)) { return true; }

        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mkDir");
        exit(1);
    }

}
//***********************************************************************
map<string, vector<string> > Utils::parseClasses(string classes){
    try {
        map<string, vector<string> > parts;

        //treatment<Early|Late>-age<young|old>
        vector<string> pieces; splitAtDash(classes, pieces); // -> treatment<Early|Late>, age<young|old>

        for (int i = 0; i < pieces.size(); i++) {
            string category = ""; string value = "";
            bool foundOpen = false;
            for (int j = 0; j < pieces[i].length(); j++) {
                if (m->getControl_pressed()) { return parts; }

                if (pieces[i][j] == '<')        { foundOpen = true;         }
                else if (pieces[i][j] == '>')   { j += pieces[i].length();  }
                else {
                    if (!foundOpen) { category += pieces[i][j]; }
                    else { value += pieces[i][j]; }
                }
            }
            vector<string> values; splitAtChar(value, values, '|');
            parts[category] = values;
        }

        return parts;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "parseClasses");
        exit(1);
    }
}
/***********************************************************************/
string Utils::hasPath(string longName){
    try {
        string path = "";
        size_t found;
        found=longName.find_last_of("~/\\");

        if(found != longName.npos){ path = longName.substr(0, found+1); }

        return path;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "hasPath");
        exit(1);
    }
}
/***********************************************************************/
void Utils::getCurrentDate(string& thisYear, string& thisMonth, string& thisDay){
    try {
        time_t rawtime;
        struct tm * timeinfo;

        time (&rawtime);
        timeinfo = localtime(&rawtime);

        char buffer[80];
        strftime(buffer,sizeof(buffer),"%Y",timeinfo);
        string year(buffer); thisYear = year;

        strftime(buffer,sizeof(buffer),"%m",timeinfo);
        string Month(buffer); thisMonth = Month;

        strftime(buffer,sizeof(buffer),"%d",timeinfo);
        string Day(buffer); thisDay = Day;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getCurrentDate");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::isUTF_8(string& input){
    try {
        string::iterator end_it = utf8::find_invalid(input.begin(), input.end());

        if (end_it != input.end()) {

            //try to convert it
            // Get the line length (at least for the valid part)
            //int length = utf8::distance(input.begin(), end_it);
            m->mothurOut("[WARNING]: trying to convert " + input + " to UTF-8 format..."); cout.flush();

            // Convert it to utf-16
            vector<unsigned short> utf16line;
            utf8::utf8to16(input.begin(), end_it, back_inserter(utf16line));

            // And back to utf-8
            string utf8line;
            utf8::utf16to8(utf16line.begin(), utf16line.end(), back_inserter(utf8line));

            // Confirm that the conversion went OK:
            if (utf8line != string(input.begin(), end_it)) {  m->mothurOut(" unsuccessful.\n"); return false; }
            else { m->mothurOut(" successful. Replacing " + input + " with " + utf8line + ".\n"); input = utf8line; }

        }

        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isUTF_8");
        exit(1);
    }
}
/***********************************************************************/
string Utils::getExtension(string longName){
    try {
        string extension = "";

        if(longName.find_last_of('.') != longName.npos){
            int pos = longName.find_last_of('.');
            extension = longName.substr(pos, longName.length());
        }

        return extension;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getExtension");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::isBlank(string fileName){
    try {

        fileName = getFullPathName(fileName);

        ifstream fileHandle;
        fileHandle.open(fileName.c_str());
        if(!fileHandle) { m->mothurOut("[ERROR]: Could not open " + fileName + "\n");  }
        else {  //check for blank file
            zapGremlins(fileHandle);
            gobble(fileHandle);
            if (fileHandle.eof()) { fileHandle.close(); return true;  }
            fileHandle.close();
        }
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isBlank");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::stringBlank(string input){
    try {
        for (int i = 0; i < input.length(); i++) { if (!isspace(input[i])) { return false; } }
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isBlank");
        exit(1);
    }
}
/**************************************************************************************************/
vector<unsigned long long> Utils::setFilePosFasta(string filename, long long& num, char delim) {
    try {
        vector<unsigned long long> positions;
        ifstream inFASTA;
        string completeFileName = getFullPathName(filename);
        inFASTA.open(completeFileName.c_str(), ios::binary);
        int nameLine = 2;
        if (delim == '@') { nameLine = 4; }
        else if (delim == '>') { nameLine = 2; }
        else { m->mothurOut("[ERROR]: unknown file deliminator, quitting.\n"); m->setControl_pressed(true); }

        unsigned long long count = 0;
        long long numLines = 0;
        while(!inFASTA.eof()){
            char c = inFASTA.get(); count++;
            string input = ""; input += c;
            while ((c != '\n') && (c != '\r') && (c != '\f') && (c != EOF)) {
                c = inFASTA.get(); count++;
                input += c;
            }
            numLines++;
            //gobble
            while(isspace(c=inFASTA.get()))		{ input += c; count++;}
            if(!inFASTA.eof()) { inFASTA.putback(c); count--;  }


            if (input.length() != 0) {
                if((input[0] == delim) && (((numLines-1)%nameLine) == 0)){ //this is a name line
                    positions.push_back(count+numLines-input.length());
                }else if (int(c) == -1) { break; }
                else { input = ""; }
            }
        }
        inFASTA.close();

        num = positions.size();

        FILE * pFile;
        unsigned long long size;

        //get num bytes in file
        pFile = fopen (completeFileName.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

        positions.push_back(size);
        positions[0] = 0;

        return positions;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "setFilePosFasta");
        exit(1);
    }
}
/**************************************************************************************************/
vector<unsigned long long> Utils::setFilePosFasta(string filename, long long& num) {
    try {
        vector<unsigned long long> positions;
        ifstream inFASTA;
        //openInputFileBinary(filename, inFASTA);
        string completeFileName = getFullPathName(filename);
        inFASTA.open(completeFileName.c_str(), ios::binary);

        string input;
        unsigned long long count = 0;
        while(!inFASTA.eof()){
            char c = inFASTA.get(); count++;
            if (c == '>') { positions.push_back(count-1); }
        }
        inFASTA.close();

        num = positions.size();

        FILE * pFile;
        unsigned long long size;

        //get num bytes in file
        pFile = fopen (completeFileName.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

        positions.push_back(size);
        positions[0] = 0;

        return positions;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "setFilePosFasta");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<consTax> Utils::readConsTax(string inputfile){
    try {

        vector<consTax> taxes;

        ifstream in;
        openInputFile(inputfile, in);

        //read headers
        getline(in);

        while (!in.eof()) {

            if (m->getControl_pressed()) { break; }

            string otu = ""; string tax = "unknown";
            int size = 0;

            in >> otu; gobble(in);
            in >> size; gobble(in);
            tax = getline(in); gobble(in);

            consTax temp(otu, tax, size);
            taxes.push_back(temp);
        }
        in.close();

        return taxes;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readConsTax");
        exit(1);
    }
}
//**********************************************************************************************************************
int Utils::readConsTax(string inputfile, map<int, consTax2>& taxes){
    try {
        ifstream in;
        openInputFile(inputfile, in);

        //read headers
        getline(in);

        while (!in.eof()) {

            if (m->getControl_pressed()) { break; }

            string otu = ""; string tax = "unknown";
            int size = 0;

            in >> otu; gobble(in);
            in >> size; gobble(in);
            tax = getline(in); gobble(in);

            consTax2 temp(otu, tax, size);
            string simpleBin = getSimpleLabel(otu);
            int bin;
            convert(simpleBin, bin);
            taxes[bin] = temp;
        }
        in.close();

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readConsTax");
        exit(1);
    }
}
/**************************************************************************************************/
vector<unsigned long long> Utils::setFilePosEachLine(string filename, long long& num) {
    try {
        filename = getFullPathName(filename);

        vector<unsigned long long> positions;
        ifstream in;
        //openInputFile(filename, in);
        openInputFileBinary(filename, in);

        string input;
        unsigned long long count = 0;
        positions.push_back(0);

        while(!in.eof()){
            //getline counting reads
            char d = in.get(); count++;
            while ((d != '\n') && (d != '\r') && (d != '\f') && (d != in.eof()))	{
                //get next character
                d = in.get();
                count++;
            }

            if (!in.eof()) {
                d=in.get(); count++;
                while(isspace(d) && (d != in.eof()))		{ d=in.get(); count++;}
            }
            positions.push_back(count-1);
            
        }
        in.close();

        num = positions.size()-1;

        FILE * pFile;
        unsigned long long size;

        //get num bytes in file
        pFile = fopen (filename.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

        positions[(positions.size()-1)] = size;

        return positions;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "setFilePosEachLine");
        exit(1);
    }
}
/**************************************************************************************************/
vector<unsigned long long> Utils::setFilePosEachLine(string filename, unsigned long long& num) {
    try {
        filename = getFullPathName(filename);

        vector<unsigned long long> positions;
        ifstream in;
        //openInputFile(filename, in);
        openInputFileBinary(filename, in);

        string input;
        unsigned long long count = 0;
        positions.push_back(0);

        while(!in.eof()){
            //getline counting reads
            char d = in.get(); count++;
            while ((d != '\n') && (d != '\r') && (d != '\f') && (d != in.eof()))	{
                //get next character
                d = in.get();
                count++;
            }

            if (!in.eof()) {
                d=in.get(); count++;
                while(isspace(d) && (d != in.eof()))		{ d=in.get(); count++;}
            }
            positions.push_back(count-1);
        }
        in.close();

        num = positions.size()-1;

        FILE * pFile;
        unsigned long long size;

        //get num bytes in file
        pFile = fopen (filename.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

        positions[(positions.size()-1)] = size;

        return positions;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "setFilePosEachLine");
        exit(1);
    }
}

/**************************************************************************************************/

vector<unsigned long long> Utils::divideFile(string filename, int& proc) {
    try{
        vector<unsigned long long> filePos;
        filePos.push_back(0);

        FILE * pFile;
        unsigned long long size;

        filename = getFullPathName(filename);

        //get num bytes in file
        pFile = fopen (filename.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

#if defined NON_WINDOWS

        //estimate file breaks
        unsigned long long chunkSize = 0;
        chunkSize = size / proc;

        //file to small to divide by processors
        if (chunkSize == 0)  {  proc = 1;	filePos.push_back(size); return filePos;	}

        if (proc > 1) {
            //for each process seekg to closest file break and search for next '>' char. make that the filebreak
            for (int i = 0; i < proc; i++) {
                unsigned long long spot = (i+1) * chunkSize;

                ifstream in;
                openInputFile(filename, in);
                in.seekg(spot);

                //look for next '>'
                unsigned long long newSpot = spot;
                while (!in.eof()) {
                    char c = in.get();

                    if (c == '>') {   in.putback(c); newSpot = in.tellg(); break;  }
                    else if (int(c) == -1) { break; }

                }

                //there was not another sequence before the end of the file
                unsigned long long sanityPos = in.tellg();

                if (sanityPos == -1) {	break;  }
                else {  filePos.push_back(newSpot);  }

                in.close();
            }
        }
        //save end pos
        filePos.push_back(size);

        //sanity check filePos
        for (int i = 0; i < (filePos.size()-1); i++) {
            if (filePos[(i+1)] <= filePos[i]) {  filePos.erase(filePos.begin()+(i+1)); i--; }
        }

        proc = (filePos.size() - 1);
#else
        m->mothurOut("[ERROR]: Windows version should not be calling the divideFile function.\n");
        proc=1;
        filePos.push_back(size);
#endif
        return filePos;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "divideFile");
        exit(1);
    }
}
/**************************************************************************************************/

vector<unsigned long long> Utils::divideFile(string filename, int& proc, char delimChar) {
    try{
        vector<unsigned long long> filePos;
        filePos.push_back(0);

        FILE * pFile;
        unsigned long long size;

        filename = getFullPathName(filename);

        //get num bytes in file
        pFile = fopen (filename.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

        char secondaryDelim = '>';
        if (delimChar == '@') { secondaryDelim = '+'; }

#if defined NON_WINDOWS

        //estimate file breaks
        unsigned long long chunkSize = 0;
        chunkSize = size / proc;

        //file to small to divide by processors
        if (chunkSize == 0)  {  proc = 1;	filePos.push_back(size); return filePos;	}

        //for each process seekg to closest file break and search for next delimChar char. make that the filebreak
        for (int i = 0; i < proc; i++) {
            unsigned long long spot = (i+1) * chunkSize;

            ifstream in;
            openInputFile(filename, in);
            in.seekg(spot);

            getline(in); //get to end of line in case you jump into middle of line where the delim char happens to fall.

            //look for next delimChar
            unsigned long long newSpot = spot;
            while (!in.eof()) {
                char c = in.get();
                string input = ""; input += c;
                while ((c != '\n') && (c != '\r') && (c != '\f') && (c != EOF)) {
                    c = in.get();
                    input += c;
                }

                if (input.length() != 0) {
                    if(input[0] == delimChar){ //this is a potential name line
                        newSpot = in.tellg();
                        newSpot -=input.length();
                        //get two lines and look for secondary delim
                        //inf a fasta file this would be a new sequence, in fastq it will be the + line, if this was a nameline.
                        getline(in); gobble(in);
                        if (!in.eof()) {
                            string secondInput = getline(in); gobble(in);
                            if (secondInput[0] == secondaryDelim) { break; } //yes, it was a nameline so stop
                            else { input = ""; gobble(in); } //nope it was a delim at the beginning of a non nameline, keep looking.
                        }
                    }else if (int(c) == -1) { break; }
                    else {  input = ""; gobble(in); }
                }
            }

            //there was not another sequence before the end of the file
            unsigned long long sanityPos = in.tellg();

            if (sanityPos == -1) {	break;  }
            else {  filePos.push_back(newSpot);  }

            in.close();
        }

        //save end pos
        filePos.push_back(size);

        //sanity check filePos
        for (int i = 0; i < (filePos.size()-1); i++) {
            if (filePos[(i+1)] <= filePos[i]) {  filePos.erase(filePos.begin()+(i+1)); i--; }
        }

        proc = (filePos.size() - 1);
#else
        m->mothurOut("[ERROR]: Windows version should not be calling the divideFile function.\n");
        proc=1;
        filePos.push_back(size);
#endif
        return filePos;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "divideFile");
        exit(1);
    }
}

/**************************************************************************************************/

vector<unsigned long long> Utils::divideFilePerLine(string filename, int& proc) {
    try{
        vector<unsigned long long> filePos;
        filePos.push_back(0);

        FILE * pFile;
        unsigned long long size;

        filename = getFullPathName(filename);

        //get num bytes in file
        pFile = fopen (filename.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }

#if defined NON_WINDOWS
        //estimate file breaks
        unsigned long long chunkSize = 0;
        chunkSize = size / proc;

        //file to small to divide by processors
        if (chunkSize == 0)  {  proc = 1;	filePos.push_back(size); return filePos;	}

        //for each process seekg to closest file break and search for next '>' char. make that the filebreak
        for (int i = 0; i < proc; i++) {
            unsigned long long spot = (i+1) * chunkSize;

            ifstream in;
            openInputFile(filename, in);
            in.seekg(spot);

            //look for next line break
            unsigned long long newSpot = spot;
            while (!in.eof()) {
                char c = in.get();

                if ((c == '\n') || (c == '\r') || (c == '\f'))	{ gobble(in); newSpot = in.tellg(); break; }
                else if (int(c) == -1) { break; }
            }

            //there was not another line before the end of the file
            unsigned long long sanityPos = in.tellg();

            if (sanityPos == -1) {	break;  }
            else {  filePos.push_back(newSpot);  }

            in.close();
        }

        //save end pos
        filePos.push_back(size);

        //sanity check filePos
        for (int i = 0; i < (filePos.size()-1); i++) {
            if (filePos[(i+1)] <= filePos[i]) {  filePos.erase(filePos.begin()+(i+1)); i--; }
        }

        proc = (filePos.size() - 1);
#else
        m->mothurOut("[ERROR]: Windows version should not be calling the divideFile function.\n");
        proc=1;
        filePos.push_back(size);
#endif
        return filePos;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "divideFile");
        exit(1);
    }
}
/**************************************************************************************************/
int Utils::divideFile(string filename, int& proc, vector<string>& files) {
    try{

        vector<unsigned long long> filePos = divideFile(filename, proc);

        for (int i = 0; i < (filePos.size()-1); i++) {

            //read file chunk
            ifstream in;
            openInputFile(filename, in);
            in.seekg(filePos[i]);
            unsigned long long size = filePos[(i+1)] - filePos[i];
            char* chunk = new char[size];
            in.read(chunk, size);
            in.close();

            //open new file
            string fileChunkName = filename + "." + toString(i) + ".tmp";
            ofstream out;
            openOutputFile(fileChunkName, out);

            out << chunk << endl;
            out.close();
            delete[] chunk;

            //save name
            files.push_back(fileChunkName);
        }

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "divideFile");
        exit(1);
    }
}
/***********************************************************************/

bool Utils::isTrue(string f){
    try {

        for (int i = 0; i < f.length(); i++) { f[i] = toupper(f[i]); }

        if ((f == "TRUE") || (f == "T")) {	return true;	}
        else {	return false;  }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isTrue");
        exit(1);
    }
}

/***********************************************************************/

float Utils::roundDist(float dist, int precision){
    try {
        return int(dist * precision + 0.5)/float(precision);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "roundDist");
        exit(1);
    }
}
/***********************************************************************/

float Utils::ceilDist(float dist, int precision){
    try {
        return int(ceil(dist * precision))/float(precision);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "ceilDist");
        exit(1);
    }
}
/***********************************************************************/

vector<string> Utils::splitWhiteSpace(string& rest, char buffer[], int size){
    try {
        vector<string> pieces;

        for (int i = 0; i < size; i++) {
            if (!isspace(buffer[i]))  { rest += buffer[i];  }
            else {
                if (rest != "") { pieces.push_back(rest);  rest = ""; }
                while (i < size) {  //gobble white space
                    if (isspace(buffer[i])) { i++; }
                    else { rest = buffer[i];  break; }
                }
            }
        }

        return pieces;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitWhiteSpace");
        exit(1);
    }
}
/***********************************************************************/
vector<string> Utils::splitWhiteSpace(string input){
    try {
        vector<string> pieces;
        string rest = "";

        for (int i = 0; i < input.length(); i++) {
            if (!isspace(input[i]))  { rest += input[i];  }
            else {
                if (rest != "") { pieces.push_back(rest);  rest = ""; }
                while (i < input.length()) {  //gobble white space
                    if (isspace(input[i])) { i++; }
                    else { rest = input[i];  break; }
                }
            }
        }

        if (rest != "") { pieces.push_back(rest); }

        return pieces;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitWhiteSpace");
        exit(1);
    }
}
/***********************************************************************/
int Utils::splitWhiteSpace(string input, vector<float>& pieces, int index){
    try {
        pieces.clear();
        string rest = "";
        int count = 0;

        for (int i = 0; i < input.length(); i++) {
            if (!isspace(input[i]))  { rest += input[i];  }
            else {
                if (rest != "") { float tdist; mothurConvert(rest, tdist); pieces.push_back(tdist); count++; rest = ""; }
                while (i < input.length()) {  //gobble white space
                    if (isspace(input[i])) { i++; }
                    else { rest = input[i];  break; }
                }
                if (count > index) { return 0; }
            }
        }

        if (rest != "") { float tdist; mothurConvert(rest, tdist); count++; pieces.push_back(tdist); }

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitWhiteSpace");
        exit(1);
    }
}
/***********************************************************************/
vector<string> Utils::splitWhiteSpaceWithQuotes(string input){
    try {
        vector<string> pieces;
        string rest = "";

        int pos = input.find('\'');
        int pos2 = input.find('\"');

        if ((pos == string::npos) && (pos2 == string::npos)) { return splitWhiteSpace(input); } //no quotes to worry about
        else {
            for (int i = 0; i < input.length(); i++) {
                if ((input[i] == '\'') || (input[i] == '\"') || (rest == "\'") || (rest == "\"")) { //grab everything til end or next ' or "
                    rest += input[i];
                    for (int j = i+1; j < input.length(); j++) {
                        if ((input[j] == '\'') || (input[j] == '\"')) {  //then quit
                            rest += input[j];
                            i = j+1;
                            j+=input.length();
                        }else { rest += input[j]; }
                    }
                }else if (!isspace(input[i]))  { rest += input[i];  }
                else {
                    if (rest != "") { pieces.push_back(rest);  rest = ""; }
                    while (i < input.length()) {  //gobble white space
                        if (isspace(input[i])) { i++; }
                        else { rest = input[i];  break; } 
                    }
                }
            }

            if (rest != "") { pieces.push_back(rest); }
        }
        return pieces;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitWhiteSpace");
        exit(1);
    }
}
//**********************************************************************************************************************
int Utils::readTax(string taxfile, map<string, string>& taxMap, bool removeConfidence) {
    try {
        //open input file
        ifstream in;
        openInputFile(taxfile, in);

        bool error = false;
        string name, taxonomy;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in >> name; gobble(in);
            taxonomy = getline(in); gobble(in);

            checkName(name);
            
            //are there confidence scores, if so remove them
            if (removeConfidence) {  if (taxonomy.find_first_of('(') != -1) {  removeConfidences(taxonomy);	} }
            map<string, string>::iterator itTax = taxMap.find(name);

            if(itTax == taxMap.end()) {
                bool ignore = false;
                if (taxonomy != "") { if (taxonomy[taxonomy.length()-1] != ';') { m->mothurOut("[ERROR]: " + name + " is missing the final ';', ignoring.\n"); ignore=true; }
                }
                if (!ignore) { taxMap[name] = taxonomy; }
            }else { m->mothurOut("[ERROR]: " + name + " is already in your taxonomy file, names must be unique./n"); error = true; }
        }
        in.close();

        if (error) { m->setControl_pressed(true); }

        return taxMap.size();

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readTax");
        exit(1);
    }
}
/**********************************************************************************************************************/
int Utils::readNames(string namefile, map<string, string>& nameMap, bool redund) {
    try {
        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);

                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    pairDone = false;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);

                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {   nameMap[theseNames[i]] = firstCol;  }
                    pairDone = false;
                }
            }
        }

        return nameMap.size();

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
/**********************************************************************************************************************/
int Utils::readNames(string namefile, map<string, string>& nameMap, int flip) {
    try {
        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[secondCol] = firstCol;
                    pairDone = false;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[secondCol] = firstCol;
                    pairDone = false;
                }
            }
        }

        return nameMap.size();

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
/**********************************************************************************************************************/
int Utils::readNames(string namefile, map<string, string>& nameMap, map<string, int>& nameCount) {
    try {
        nameMap.clear(); nameCount.clear();
        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    nameCount[firstCol] = theseNames.size();
                    pairDone = false;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    nameCount[firstCol] = theseNames.size();
                    pairDone = false;
                }
            }

        }
        return nameMap.size();

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
/**********************************************************************************************************************/
int Utils::readNames(string namefile, map<string, string>& nameMap) {
    try {
        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[firstCol] = secondCol; pairDone = false; }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[firstCol] = secondCol; pairDone = false; }
            }
        }

        return nameMap.size();

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
/**********************************************************************************************************************/
int Utils::readNames(string namefile, map<string, string>& nameMap, set<string>& namesToInclude) {
    try {
        //open input file
        ifstream in;
        openInputFile(namefile, in);
        
        string firstCol, secondCol;
        
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            in >> firstCol; gobble(in);
            in >> secondCol; gobble(in);
            
            checkName(firstCol);
            checkName(secondCol);
            
            vector<string> secondNames; splitAtComma(secondCol, secondNames);
            
            secondCol = ""; firstCol = "";
            
            for (int i = 0; i < secondNames.size(); i++) {
                if (namesToInclude.count(secondNames[i]) != 0) { //we want to include you
                    secondCol += secondNames[i] + ",";
                    if (firstCol == "") {   firstCol = secondNames[i]; }
                }
            }
            
            if (secondCol != "") {
                //remove last comma
                secondCol = secondCol.substr(0,secondCol.length()-1);
            
                nameMap[firstCol] = secondCol;
            }
            
        }
        in.close();
        
        
        return nameMap.size();
        
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}

/**********************************************************************************************************************/
int Utils::readNames(string namefile, map<string, vector<string> >& nameMap) {
    try {
        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    vector<string> temp;
                    splitAtComma(secondCol, temp);
                    nameMap[firstCol] = temp;
                    pairDone = false;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    vector<string> temp;
                    splitAtComma(secondCol, temp);
                    nameMap[firstCol] = temp;
                    pairDone = false;
                }
            }
        }

        return nameMap.size();
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
/**********************************************************************************************************************/
map<string, int> Utils::readNames(string namefile) {
    try {
        map<string, int> nameMap;

        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;
                }
            }
        }

        return nameMap;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
/**********************************************************************************************************************/
void Utils::readNames(string namefile, map<string, long long>& nameMap) {
    try {
        //open input file
        ifstream in; openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    long long num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;
                }
            }
        }
        in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    long long num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;
                }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}

/**********************************************************************************************************************/
map<string, int> Utils::readNames(string namefile, unsigned long int& numSeqs) {
    try {
        map<string, int> nameMap;
        numSeqs = 0;

        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;
                    numSeqs += num;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;
                    numSeqs += num;
                }
            }
        }

        return nameMap;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
//**********************************************************************************************************************
int Utils::printVsearchFile(vector<seqPriorityNode>& nameMapCount, string filename, string tag, string tag2){
    try {

        sort(nameMapCount.begin(), nameMapCount.end(), compareSeqPriorityNodes);

        ofstream out;
        openOutputFile(filename, out);

        //print new file in order of
        for (int i = 0; i < nameMapCount.size(); i++) {
            if (m->getControl_pressed()) {break;}
            out << ">" << nameMapCount[i].name  << tag << nameMapCount[i].numIdentical << tag2 << endl << nameMapCount[i].seq << endl;
        }
        out.close();

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "printVsearchFile");
        exit(1);
    }
}
/************************************************************/
int Utils::checkName(string& name) {
    try {
        if (modifyNames) {
            for (int i = 0; i < name.length(); i++) {
                if (name[i] == ':') { name[i] = '_'; m->setChangedSeqNames(true); }
            }
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "checkName");
        exit(1);
    }
}
/************************************************************/
bool Utils::checkGroupName(string name) {
    try {

        bool goodName = true;
        for (int i = 0; i < name.length(); i++) {
            if (name[i] == ':') {  goodName = false; break;  }
            else if (name[i] == '-') {  goodName = false; break;  }
            else if (name[i] == '/') {  goodName = false; break;  }
        }

        if (!goodName) {
            m->mothurOut("\n[WARNING]: group " + name + " contains illegal characters in the name. Group names should not include :, -, or / characters.  The ':' character is a special character used in trees. Using ':' will result in your tree being unreadable by tree reading software.  The '-' character is a special character used by mothur to parse group names.  Using the '-' character will prevent you from selecting groups. The '/' character will created unreadable filenames when mothur includes the group in an output filename.  \n\n");
        }

        return goodName;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "checkGroupName");
        exit(1);
    }
}
/**********************************************************************************************************************/
int Utils::readNames(string namefile, vector<seqPriorityNode>& nameVector, map<string, string>& fastamap) {
    try {
        int error = 0;

        //open input file
        ifstream in;
        openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);

                    map<string, string>::iterator it = fastamap.find(firstCol);
                    if (it == fastamap.end()) {
                        error = 1;
                        m->mothurOut("[ERROR]: " + firstCol + " is not in your fastafile, but is in your namesfile, please correct.\n");
                    }else {
                        seqPriorityNode temp(num, it->second, firstCol);
                        nameVector.push_back(temp);
                    }

                    pairDone = false;
                }
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);

            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }

                if (pairDone) {
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);

                    map<string, string>::iterator it = fastamap.find(firstCol);
                    if (it == fastamap.end()) {
                        error = 1;
                        m->mothurOut("[ERROR]: " + firstCol + " is not in your fastafile, but is in your namesfile, please correct.\n");
                    }else {
                        seqPriorityNode temp(num, it->second, firstCol);
                        nameVector.push_back(temp);
                    }

                    pairDone = false;
                }
            }
        }
        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readNames");
        exit(1);
    }
}
//**********************************************************************************************************************
set<string> Utils::readAccnos(string accnosfile){
    try {
        set<string> names;
        ifstream in;
        bool ableToOpen = openInputFile(accnosfile, in, "");
        if (!ableToOpen) {  m->mothurOut("[ERROR]: Could not open " + accnosfile + "\n"); return names; }
        string name;

        string rest = "";
        char buffer[4096];

        unsigned long long count = 0;
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]);
                names.insert(pieces[i]);
                count++;
            }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.insert(pieces[i]);  count++; }
        }

        return names;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readAccnos");
        exit(1);
    }
}
//**********************************************************************************************************************
int Utils::readAccnos(string accnosfile, vector<string>& names){
    try {
        names.clear();
        ifstream in;
        openInputFile(accnosfile, in);
        string name;

        string rest = "";
        char buffer[4096];

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readAccnos");
        exit(1);
    }
}
//**********************************************************************************************************************
int Utils::readAccnos(string accnosfile, vector<string>& names, string noerror){
    try {
        names.clear();
        ifstream in;
        openInputFile(accnosfile, in, noerror);
        string name;

        string rest = "";
        char buffer[4096];

        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }

            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());

            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }
        in.close();

        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readAccnos");
        exit(1);
    }
}
/***********************************************************************/

int Utils::getNumNames(string names){
    try {
        int count = 0;

        if(names != ""){
            count = 1;
            for(int i=0;i<names.size();i++){
                if(names[i] == ','){
                    count++;
                }
            }
        }

        return count;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getNumNames");
        exit(1);
    }
}
/***********************************************************************/

int Utils::getNumChar(string line, char c){
    try {
        int count = 0;

        if(line != ""){
            for(int i=0;i<line.size();i++){
                if(line[i] == c){
                    count++;
                }
            }
        }

        return count;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getNumChar");
        exit(1);
    }
}
/***********************************************************************/
string Utils::getSimpleLabel(string label){
    try {
        string simple = "";

        //remove OTU or phylo tag
        string newLabel1 = "";
        for (int i = 0; i < label.length(); i++) {
            if(label[i]>47 && label[i]<58) { //is a digit
                newLabel1 += label[i];
            }
        }

        int num1;

        mothurConvert(newLabel1, num1);

        simple = toString(num1);

        return simple;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getSimpleLabel");
        exit(1);
    }
}
/***********************************************************************/

bool Utils::isLabelEquivalent(string label1,  string label2){
    try {
        bool same = false;

        //remove OTU or phylo tag
        string newLabel1 = "";
        for (int i = 0; i < label1.length(); i++) {
            if(label1[i]>47 && label1[i]<58) { //is a digit
                newLabel1 += label1[i];
            }
        }

        string newLabel2 = "";
        for (int i = 0; i < label2.length(); i++) {
            if(label2[i]>47 && label2[i]<58) { //is a digit
                newLabel2 += label2[i];
            }
        }

        int num1, num2;
        mothurConvert(newLabel1, num1);
        mothurConvert(newLabel2, num2);

        if (num1 == num2) { same = true; }

        return same;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isLabelEquivalent");
        exit(1);
    }
}
//**********************************************************************************************************************
bool Utils::isSubset(vector<string> bigset, vector<string> subset) {
    try {


        if (subset.size() > bigset.size()) { return false;  }

        //check if each guy in subset is also in bigset
        for (int i = 0; i < subset.size(); i++) {
            bool match = false;
            for (int j = 0; j < bigset.size(); j++) {
                if (subset[i] == bigset[j]) { match = true; break; }
            }

            //you have a guy in subset that had no match in bigset
            if (!match) { return false; }
        }

        return true;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isSubset");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::mothurRemove(string filename){
    try {
        filename = getFullPathName(filename);
        int error = remove(filename.c_str());
        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurRemove");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::mothurConvert(string item, int& num){
    try {
        bool error = false;

        if (isNumeric1(item)) { convert(item, num); }
        else {
            num = 0;
            error = true;
            m->mothurOut("[ERROR]: cannot convert " + item + " to an integer.\n");
            m->setControl_pressed(true);
        }

        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-int");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::mothurConvert(char item, int& num){
    try {
        bool error = false;

        if (isdigit(item)) {
            string mystring; mothurConvert(item, mystring);
            mothurConvert(mystring, num);
        }else {
            num = 0;
            error = true;
            m->mothurOut("[ERROR]: cannot convert " + toString(item) + " to an integer.\n");
            m->setControl_pressed(true);
        }

        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-int");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::mothurConvert(char item, string& output){
    try {

        stringstream ss;
        ss << item;
        ss >> output;
        return true;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-char");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::mothurConvert(string item, intDist& num){
    try {
        bool error = false;

        if (isNumeric1(item)) {
            convert(item, num);
        }else {
            num = 0;
            error = true;
            m->mothurOut("[ERROR]: cannot convert " + item + " to an integer.\n");
            m->setControl_pressed(true);
        }

        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-intDist");
        exit(1);
    }
}
/***********************************************************************/
set<string> Utils::mothurConvert(vector<string>& input){
    try {
        set<string> output(input.begin(), input.end());
        
        
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-vectorToSet");
        exit(1);
    }
}
/***********************************************************************/
vector<string> Utils::mothurConvert(set<string>& input){
    try {
        vector<string> output(input.begin(), input.end());
        
        
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-SetToVector");
        exit(1);
    }
}
/**************************************************************************************************/
string Utils::addUnclassifieds(string tax, int maxlevel, bool probs) {
    try{
        string newTax, taxon;

        string savedTax = tax;
        vector<string> taxons; splitAtChar(tax, taxons, ';'); taxons.pop_back();
        vector<int> confidences;

        if (taxons.size() == maxlevel) { return savedTax; }

        int index = 0;
        int confidence = 0;
        int level = 1;
        for (int i = 0; i < taxons.size(); i++) {
            index = i;
            string thisTax = taxons[i]+";";
            confidence = removeConfidences(thisTax);
            confidences.push_back(confidence);

            if (thisTax == "unclassified;"){ index--; break; }
            else{ newTax += taxons[i] + ";";  }
        }
        level = index+1;

        string thisTax = taxons[index]+";";

        removeConfidences(thisTax);
        taxon = thisTax.substr(0, thisTax.length()-1);

        string cTax = "";
        if (probs)  { cTax = taxon + "_unclassified(" + toString(confidences[index]) + ");";     }
        else        { cTax = taxon + "_unclassified;";          }

        //add "unclassified" until you reach maxLevel
        while (level < maxlevel) {
            newTax += cTax;
            level++;
        }

        return newTax;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "addUnclassifieds");
        exit(1);
    }
}
/**************************************************************************************************/
string Utils::trimTax(string tax, int trimLevel) {
    try{
        string newTax = "";
        string savedTax = tax;
        vector<string> taxons; splitAtChar(tax, taxons, ';'); taxons.pop_back();
    
        if (taxons.size() == trimLevel) { return savedTax; }
        else {
            int level = 0;
            for (int i = 0; i < taxons.size(); i++) {
                newTax += taxons[i] +";";
                level++;
                if (level == trimLevel) { break; }
            }
        }
        
        return newTax;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "trimTax");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::isNumeric1(string stringToCheck){
    try {
        bool numeric = false;

        if (stringToCheck == "") { numeric = false;  }
        else if(stringToCheck.find_first_not_of("0123456789.-") == string::npos) { numeric = true; }

        return numeric;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isNumeric1");
        exit(1);
    }

}
/***********************************************************************/
bool Utils::isPositiveNumeric(string stringToCheck){
    try {
        bool numeric = false;
        
        if (stringToCheck == "") { numeric = false;  }
        else if(stringToCheck.find_first_not_of("0123456789.") == string::npos) { numeric = true; }
        
        return numeric;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isPositiveNumeric");
        exit(1);
    }
    
}
/***********************************************************************/
bool Utils::isEqual(float num1, float num2){
    try {
        bool equal = false;
        
        if (fabs(num1-num2) <= fabs(num1 * 0.001)) { equal = true; }
        
        return equal;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isEqual");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::isEqual(double num1, double num2){
    try {
        bool equal = false;
        
        if (fabs(num1-num2) <= fabs(num1 * 0.001)) { equal = true; }
        
        return equal;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isEqual");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::allSpaces(string stringToCheck){
    try {

        for (int i = 0; i < stringToCheck.length(); i++) {
            char c = stringToCheck[i];
            if (!isspace(c)) { return false; }
        }

        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isNumeric1");
        exit(1);
    }

}
/***********************************************************************/
bool Utils::isInteger(string stringToCheck){
    try {
        bool isInt = false;

        if(stringToCheck.find_first_not_of("0123456789-") == string::npos) { isInt = true; }

        return isInt;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isInteger");
        exit(1);
    }

}
/***********************************************************************/
bool Utils::containsAlphas(string stringToCheck){
    try {
        bool containsAlpha = false;

        if(stringToCheck.find_first_of("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz") != string::npos) { containsAlpha = true; }

        return containsAlpha;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "containsAlphas");
        exit(1);
    }

}
/***********************************************************************/
bool Utils::isAllAlphas(string stringToCheck){
    try {
        bool allAlphas = true;

        if(stringToCheck.find_first_not_of("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz") != string::npos) { allAlphas = false; }

        return allAlphas;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isAllAlphas");
        exit(1);
    }

}
/***********************************************************************/
bool Utils::isAllAlphaNumerics(string stringToCheck){
    try {
        bool allAlphaNumerics = true;

        if(stringToCheck.find_first_not_of("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz0123456789") != string::npos) { allAlphaNumerics = false; }

        return allAlphaNumerics;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isAllAlphas");
        exit(1);
    }

}
/***********************************************************************/
bool Utils::mothurConvert(string item, float& num){
    try {
        bool error = false;

        if (isNumeric1(item)) {
            convert(item, num);
        }else {
            num = 0;
            error = true;
            m->mothurOut("[ERROR]: cannot convert " + item + " to a float.\n");
            m->setControl_pressed(true);
        }

        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-float");
        exit(1);
    }
}
/***********************************************************************/
bool Utils::mothurConvert(string item, double& num){
    try {
        bool error = false;

        if (isNumeric1(item)) {
            convert(item, num);
        }else {
            num = 0;
            error = true;
            m->mothurOut("[ERROR]: cannot convert " + item + " to a double.\n");
            m->setControl_pressed(true);
        }

        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "mothurConvert-double");
        exit(1);
    }
}
/**************************************************************************************************/

vector<vector<double> > Utils::binomial(int maxOrder){
    try {
        vector<vector<double> > binomial(maxOrder+1);

        for(int i=0;i<=maxOrder;i++){
            binomial[i].resize(maxOrder+1);
            binomial[i][0]=1;
            binomial[0][i]=0;
        }
        binomial[0][0]=1;

        binomial[1][0]=1;
        binomial[1][1]=1;

        for(int i=2;i<=maxOrder;i++){
            binomial[1][i]=0;
        }

        for(int i=2;i<=maxOrder;i++){
            for(int j=1;j<=maxOrder;j++){
                if(i==j){	binomial[i][j]=1;									}
                if(j>i)	{	binomial[i][j]=0;									}
                else	{	binomial[i][j]=binomial[i-1][j-1]+binomial[i-1][j];	}
            }
        }

        return binomial;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "binomial");
        exit(1);
    }
}
/**************************************************************************************************/
unsigned int Utils::fromBase36(string base36){
    try {
        unsigned int num = 0;

        map<char, int> converts;
        converts['A'] = 0;
        converts['a'] = 0;
        converts['B'] = 1;
        converts['b'] = 1;
        converts['C'] = 2;
        converts['c'] = 2;
        converts['D'] = 3;
        converts['d'] = 3;
        converts['E'] = 4;
        converts['e'] = 4;
        converts['F'] = 5;
        converts['f'] = 5;
        converts['G'] = 6;
        converts['g'] = 6;
        converts['H'] = 7;
        converts['h'] = 7;
        converts['I'] = 8;
        converts['i'] = 8;
        converts['J'] = 9;
        converts['j'] = 9;
        converts['K'] = 10;
        converts['k'] = 10;
        converts['L'] = 11;
        converts['l'] = 11;
        converts['M'] = 12;
        converts['m'] = 12;
        converts['N'] = 13;
        converts['n'] = 13;
        converts['O'] = 14;
        converts['o'] = 14;
        converts['P'] = 15;
        converts['p'] = 15;
        converts['Q'] = 16;
        converts['q'] = 16;
        converts['R'] = 17;
        converts['r'] = 17;
        converts['S'] = 18;
        converts['s'] = 18;
        converts['T'] = 19;
        converts['t'] = 19;
        converts['U'] = 20;
        converts['u'] = 20;
        converts['V'] = 21;
        converts['v'] = 21;
        converts['W'] = 22;
        converts['w'] = 22;
        converts['X'] = 23;
        converts['x'] = 23;
        converts['Y'] = 24;
        converts['y'] = 24;
        converts['Z'] = 25;
        converts['z'] = 25;
        converts['0'] = 26;
        converts['1'] = 27;
        converts['2'] = 28;
        converts['3'] = 29;
        converts['4'] = 30;
        converts['5'] = 31;
        converts['6'] = 32;
        converts['7'] = 33;
        converts['8'] = 34;
        converts['9'] = 35;

        int i = 0;
        while (i < base36.length()) {
            char c = base36[i];
            num = 36 * num + converts[c];
            i++;
        }

        return num;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "fromBase36");
        exit(1);
    }
}
/***********************************************************************/
string  Utils::findEdianness() {
    try {
        // find real endian type
        string endianType = "unknown";
        int num = 1;
        if(*(char *)&num == 1)
        {
            endianType = "LITTLE_ENDIAN";
        }
        else
        {
            endianType = "BIG_ENDIAN";
        }
        return endianType;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "findEdianness");
        exit(1);
    }
}
/***********************************************************************/
double  Utils::median(vector<double> x) {
    try {
        double value = 0.0;

        if (x.size() == 0) { } //error
        else {
            //For example, if a < b < c, then the median of the list {a, b, c} is b, and, if a < b < c < d, then the median of the list {a, b, c, d} is the mean of b and c; i.e., it is (b + c)/2.
            sort(x.begin(), x.end());
            //is x.size even?
            if ((x.size()%2) == 0) { //size() is even. median = average of 2 midpoints
                int midIndex1 = (x.size()/2)-1;
                int midIndex2 = (x.size()/2);
                value = (x[midIndex1]+ x[midIndex2]) / 2.0;
            }else {
                int midIndex = (x.size()/2);
                value = x[midIndex];
            }
        }
        return value;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "median");
        exit(1);
    }
}
/***********************************************************************/
int  Utils::median(vector<int> x) {
    try {
        double value = 0;

        if (x.size() == 0) { } //error
        else {
            //For example, if a < b < c, then the median of the list {a, b, c} is b, and, if a < b < c < d, then the median of the list {a, b, c, d} is the mean of b and c; i.e., it is (b + c)/2.
            sort(x.begin(), x.end());
            //is x.size even?
            if ((x.size()%2) == 0) { //size() is even. median = average of 2 midpoints
                int midIndex1 = (x.size()/2)-1;
                int midIndex2 = (x.size()/2);
                value = (x[midIndex1]+ x[midIndex2]) / 2.0;
            }else {
                int midIndex = (x.size()/2);
                value = x[midIndex];
            }
        }
        return (int) value;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "median - int");
        exit(1);
    }
}
/***********************************************************************/
int  Utils::average(vector<int> x) {
    try {
        int value = 0;

        for (int i = 0; i < x.size(); i++) {
            if (m->getControl_pressed()) { break; }
            value += x[i];
        }

        return ((int) value / x.size());
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "average - int");
        exit(1);
    }
}
int Utils::factorial(int num){
    try {
        int total = 1;

        for (int i = 1; i <= num; i++) {
            total *= i;
        }

        return total;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "factorial");
        exit(1);
    }
}
/***********************************************************************/
int Utils::getAlignmentLength(string file){
    try {
        ifstream in; openInputFile(file, in);
        
        Sequence seq(in);
        
        in.close();
        
        return seq.getAlignLength();
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAlignmentLength");
        exit(1);
    }
}

/***********************************************************************/

int Utils::getNumSeqs(ifstream& file){
    try {
        int numSeqs = count(istreambuf_iterator<char>(file),istreambuf_iterator<char>(), '>');
        file.seekg(0);
        return numSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getNumSeqs");
        exit(1);
    }
}
/***********************************************************************/
void Utils::getNumSeqs(ifstream& file, int& numSeqs){
    try {
        string input;
        numSeqs = 0;
        while(!file.eof()){
            input = getline(file);
            if (input.length() != 0) {
                if(input[0] == '>'){ numSeqs++;	}
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getNumSeqs");
        exit(1);
    }
}
/***********************************************************************/

//This function parses the estimator options and puts them in a vector
void Utils::splitAtChar(string& estim, vector<string>& container, char symbol) {
    try {

        if (symbol == '-') { splitAtDash(estim, container); return; }

        string individual = "";
        int estimLength = estim.size();
        for(int i=0;i<estimLength;i++){
            if(estim[i] == symbol){
                container.push_back(individual);
                individual = "";
            }
            else{
                individual += estim[i];
            }
        }
        container.push_back(individual);

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtChar");
        exit(1);
    }
}

/***********************************************************************/

//This function parses the estimator options and puts them in a vector
void Utils::splitAtDash(string& estim, vector<string>& container) {
    try {
        string individual = "";
        int estimLength = estim.size();
        bool prevEscape = false;

        for(int i=0;i<estimLength;i++){
            if(estim[i] == '-'){
                if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
                else {
                    container.push_back(individual);
                    individual = "";
                }
            }else if(estim[i] == '\\'){
                if (i < estimLength-1) {
                    if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                    else { individual += estim[i]; prevEscape = false;  } //if no, add in
                }else { individual += estim[i]; }
            }else {
                individual += estim[i];
            }
        }



        container.push_back(individual);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtDash");
        exit(1);
    }
}

/***********************************************************************/
//This function parses the label options and puts them in a set
void Utils::splitAtDash(string& estim, set<string>& container) {
    try {
        string individual = "";
        int estimLength = estim.size();
        bool prevEscape = false;

        for(int i=0;i<estimLength;i++){
            if(estim[i] == '-'){
                if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
                else {
                    container.insert(individual);
                    individual = "";
                }
            }else if(estim[i] == '\\'){
                if (i < estimLength-1) {
                    if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                    else { individual += estim[i]; prevEscape = false;  } //if no, add in
                }else { individual += estim[i]; }
            }else {
                individual += estim[i];
            }
        }
        container.insert(individual);

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtDash");
        exit(1);
    }
}
/***********************************************************************/
//This function parses the line options and puts them in a set
void Utils::splitAtDash(string& estim, set<int>& container) {
    try {
        string individual = "";
        int lineNum;
        int estimLength = estim.size();
        bool prevEscape = false;

        for(int i=0;i<estimLength;i++){
            if(estim[i] == '-'){
                if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
                else {
                    convert(individual, lineNum); //convert the string to int
                    container.insert(lineNum);
                    individual = "";
                }
            }else if(estim[i] == '\\'){
                if (i < estimLength-1) {
                    if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                    else { individual += estim[i]; prevEscape = false;  } //if no, add in
                }else { individual += estim[i]; }
            }else {
                individual += estim[i];
            }
        }

        convert(individual, lineNum); //convert the string to int
        container.insert(lineNum);
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtDash");
        exit(1);
    }
}

/***********************************************************************/
string Utils::makeList(vector<string>& names) {
    try {
        string list = "";

        if (names.size() == 0) { return list; }

        for (int i = 0; i < names.size()-1; i++) { list += names[i] + ",";  }

        //get last name
        list += names[names.size()-1];

        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "makeList");
        exit(1);
    }
}

/***********************************************************************/
//This function parses the a string and puts peices in a vector
void Utils::splitAtComma(string& estim, vector<string>& container) {
    try {
        string individual = "";
        int estimLength = estim.size();
        for(int i=0;i<estimLength;i++){
            if(estim[i] == ','){
                container.push_back(individual);
                individual = "";
            }
            else{
                individual += estim[i];
            }
        }
        container.push_back(individual);

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtComma");
        exit(1);
    }
}
/***********************************************************************/
//This function parses the a string and puts peices in a vector
void Utils::splitAtComma(string& estim, vector<int>& convertedContainer) {
    try {
        string individual = "";
        vector<string> container;
        int estimLength = estim.size();
        for(int i=0;i<estimLength;i++){
            if(estim[i] == ','){
                container.push_back(individual);
                individual = "";
            }
            else{
                individual += estim[i];
            }
        }
        container.push_back(individual);

        for (int i = 0; i < container.size(); i++) {
            int temp;
            if (mothurConvert(container[i], temp)) { convertedContainer.push_back(temp); }
        }

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtComma");
        exit(1);
    }
}
/***********************************************************************/
//This function splits up the various option parameters
void Utils::splitAtChar(string& prefix, string& suffix, char c){
    try {

        string individual = "";
        int estimLength = prefix.size();
        for(int i=0;i<estimLength;i++){
            if(prefix[i] == c){
                suffix = prefix.substr(i+1);
                prefix = individual;
                break;
            }
            else{
                individual += prefix[i];
            }
        }

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtChar");
        exit(1);
    }
}

/***********************************************************************/

//This function splits up the various option parameters
void Utils::splitAtComma(string& prefix, string& suffix){
    try {
        prefix = suffix.substr(0,suffix.find_first_of(','));
        if ((suffix.find_first_of(',')+2) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
            suffix = suffix.substr(suffix.find_first_of(',')+1, suffix.length());
            string space = " ";
            while(suffix.at(0) == ' ')
                suffix = suffix.substr(1, suffix.length());
        }else {  suffix = "";  }

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtComma");
        exit(1);
    }
}
/***********************************************************************/

//This function separates the key value from the option value i.e. dist=96_...
void Utils::splitAtEquals(string& key, string& value){
    try {
        if(value.find_first_of('=') != -1){
            key = value.substr(0,value.find_first_of('='));
            if ((value.find_first_of('=')+1) <= value.length()) {
                value = value.substr(value.find_first_of('=')+1, value.length());
            }
        }else{
            key = value;
            value = 1;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "splitAtEquals");
        exit(1);
    }
}

/**************************************************************************************************/

bool Utils::inUsersGroups(string groupname, vector<string> Groups) {
    try {
        for (int i = 0; i < Groups.size(); i++) {
            if (groupname == Groups[i]) { return true; }
        }
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "inUsersGroups");
        exit(1);
    }
}
/**************************************************************************************************/

bool Utils::inUsersGroups(vector<int> set, vector< vector<int> > sets) {
    try {
        for (int i = 0; i < sets.size(); i++) {
            if (set == sets[i]) { return true; }
        }
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "inUsersGroups");
        exit(1);
    }
}
/**************************************************************************************************/

bool Utils::inUsersGroups(int groupname, vector<int> Groups) {
    try {
        for (int i = 0; i < Groups.size(); i++) {
            if (groupname == Groups[i]) { return true; }
        }
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "inUsersGroups");
        exit(1);
    }
}

/**************************************************************************************************/
//returns true if any of the strings in first vector are in second vector
bool Utils::inUsersGroups(vector<string> groupnames, vector<string> Groups) {
    try {

        for (int i = 0; i < groupnames.size(); i++) {
            if (inUsersGroups(groupnames[i], Groups)) { return true; }
        }
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "inUsersGroups");
        exit(1);
    }
}

/**************************************************************************************************/
string Utils::getTag(string filename) {
    try {
        string tag = "Otu";
        int pos = filename.find_first_of(".tx.");
        if (pos != string::npos) { tag = "Phylo"; }
        return tag;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getTag");
        exit(1);
    }
}
/**************************************************************************************************/
//removes entries that are only white space
int Utils::removeBlanks(vector<string>& tempVector) {
    try {
        vector<string> newVector;
        for (int i = 0; i < tempVector.size(); i++) {
            bool isBlank = true;
            for (int j = 0; j < tempVector[i].length(); j++) {
                if (!isspace(tempVector[i][j])) { isBlank = false; j+= tempVector[i].length(); } //contains non space chars, break out and save
            }
            if (!isBlank) { newVector.push_back(tempVector[i]); }
        }
        tempVector = newVector;
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "removeBlanks");
        exit(1);
    }
}
/***********************************************************************/
//this function determines if the user has given us labels that are smaller than the given label.
//if so then it returns true so that the calling function can run the previous valid distance.
//it's a "smart" distance function.  It also checks for invalid labels.
bool Utils::anyLabelsToProcess(string label, set<string>& userLabels, string errorOff) {
    try {

        set<string>::iterator it;
        vector<float> orderFloat;
        map<string, float> userMap;  //the conversion process removes trailing 0's which we need to put back
        map<string, float>::iterator it2;
        float labelFloat;
        bool smaller = false;

        //unique is the smallest line
        if (label == "unique") {  return false;  }
        else {
            if (convertTestFloat(label, labelFloat)) {
                convert(label, labelFloat);
            }else { //cant convert
                return false;
            }
        }

        //go through users set and make them floats
        for(it = userLabels.begin(); it != userLabels.end();) {

            float temp;
            if ((*it != "unique") && (convertTestFloat(*it, temp) )){
                convert(*it, temp);
                orderFloat.push_back(temp);
                userMap[*it] = temp;
                it++;
            }else if (*it == "unique") {
                orderFloat.push_back(-1.0);
                userMap["unique"] = -1.0;
                it++;
            }else {
                if (errorOff == "") {  cout << (*it + " is not a valid label.\n");  }
                userLabels.erase(it++);
            }
        }

        //sort order
        sort(orderFloat.begin(), orderFloat.end());

        /*************************************************/
        //is this label bigger than any of the users labels
        /*************************************************/
        
        //loop through order until you find a label greater than label
        for (int i = 0; i < orderFloat.size(); i++) {
            if (orderFloat[i] < labelFloat) {
                smaller = true;
                if (isEqual(orderFloat[i], -1)) {
                    if (errorOff == "") { cout << ("Your file does not include the label unique.\n"); }
                    userLabels.erase("unique");
                }
                else {
                    if (errorOff == "") { cout << ("Your file does not include the label. \n");  }
                    string s = "";
                    for (it2 = userMap.begin(); it2!= userMap.end(); it2++) {
                        if (isEqual(it2->second, orderFloat[i])) {
                            s = it2->first;
                            //remove small labels
                            userLabels.erase(s);
                            break;
                        }
                    }
                    if (errorOff == "") {cout << ( s +  ". I will use the next smallest distance. \n"); }
                }
                //since they are sorted once you find a bigger one stop looking
            }else { break; }
        }

        return smaller;

    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "anyLabelsToProcess");
        exit(1);
    }
}

/**************************************************************************************************/
bool Utils::checkReleaseVersion(string line, string version) {
    try {

        bool good = true;

        //before we added this check
        if (line[0] != '#') {  good = false;  }
        else {
            //rip off #
            line = line.substr(1);

            vector<string> versionVector;
            splitAtChar(version, versionVector, '.');

            //check file version
            vector<string> linesVector;
            splitAtChar(line, linesVector, '.');

            if (versionVector.size() != linesVector.size()) { good = false; }
            else {
                for (int j = 0; j < versionVector.size(); j++) {
                    int num1, num2;
                    convert(versionVector[j], num1);
                    convert(linesVector[j], num2);

                    //if mothurs version is newer than this files version, then we want to remake it
                    if (num1 > num2) {  good = false; break;  }
                }
            }

        }
        return good;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "checkReleaseVersion");
        exit(1);
    }
}
/**************************************************************************************************/
int Utils::getTimeStamp(string filename) {
    try {
        int timeStamp = 0;

#if defined NON_WINDOWS
        struct stat st;
        int errorCode = stat (filename.c_str(), &st);
        if (errorCode != 0) {
            m->mothurOut("[ERROR]: Can't find timestamp for " + filename + "\n"); m->setControl_pressed(true);
        }else {
            timeStamp = st.st_mtime;
        }
#else
        HANDLE hFile;

        hFile = CreateFile(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL,
                           OPEN_EXISTING, 0, NULL);

        if(hFile == INVALID_HANDLE_VALUE) {
            m->mothurOut("[ERROR]: Can't find timestamp for " + filename + "\n"); m->setControl_pressed(true);
            CloseHandle(hFile); return timeStamp;
        }

        FILETIME ftCreate, ftAccess, ftWrite;
        SYSTEMTIME stUTC;
        DWORD dwRet;

        // Retrieve the file times for the file.
        bool success = GetFileTime(hFile, &ftCreate, &ftAccess, &ftWrite);

        if (success) {
            FileTimeToSystemTime(&ftWrite, &stUTC);

            tm time;
            time.tm_sec = stUTC.wSecond;
            time.tm_min = stUTC.wMinute;
            time.tm_hour = stUTC.wHour;
            time.tm_mday = stUTC.wDay;
            time.tm_mon = stUTC.wMonth - 1;
            time.tm_year = stUTC.wYear - 1900;
            time.tm_isdst = -1;
            time_t t = mktime(&time);

            timeStamp = t;
        }
        else { m->mothurOut("[ERROR]: Can't find timestamp for " + filename + "\n"); m->setControl_pressed(true); }
        CloseHandle(hFile);
#endif

        return timeStamp;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getTimeStamp");
        exit(1);
    }
}
/**************************************************************************************************/
vector<double> Utils::getAverages(vector< vector<double> >& dists) {
    try{
        vector<double> averages; //averages.resize(numComp, 0.0);
        for (int i = 0; i < dists[0].size(); i++) { averages.push_back(0.0); }

        for (int thisIter = 0; thisIter < dists.size(); thisIter++) {
            for (int i = 0; i < dists[thisIter].size(); i++) {
                averages[i] += dists[thisIter][i];
            }
        }

        //finds average.
        for (int i = 0; i < averages.size(); i++) {  averages[i] /= (double) dists.size(); }

        return averages;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverages");
        exit(1);
    }
}
/**************************************************************************************************/
double Utils::getAverage(vector<double> dists) {
    try{
        double average = 0;

        for (int i = 0; i < dists.size(); i++) {
            average += dists[i];
        }

        //finds average.
        average /= (double) dists.size();

        return average;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverage");
        exit(1);
    }
}

/**************************************************************************************************/
vector<double> Utils::getStandardDeviation(vector< vector<double> >& dists) {
    try{

        vector<double> averages = getAverages(dists);

        //find standard deviation
        vector<double> stdDev; //stdDev.resize(numComp, 0.0);
        for (int i = 0; i < dists[0].size(); i++) { stdDev.push_back(0.0); }

        for (int thisIter = 0; thisIter < dists.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int j = 0; j < dists[thisIter].size(); j++) {
                stdDev[j] += ((dists[thisIter][j] - averages[j]) * (dists[thisIter][j] - averages[j]));
            }
        }
        for (int i = 0; i < stdDev.size(); i++) {
            stdDev[i] /= (double) dists.size();
            stdDev[i] = sqrt(stdDev[i]);
        }

        return stdDev;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverages");
        exit(1);
    }
}
/**************************************************************************************************/
vector<double> Utils::getStandardDeviation(vector< vector<double> >& dists, vector<double>& averages) {
    try{
        //find standard deviation
        vector<double> stdDev; //stdDev.resize(numComp, 0.0);
        for (int i = 0; i < dists[0].size(); i++) { stdDev.push_back(0.0); }

        for (int thisIter = 0; thisIter < dists.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int j = 0; j < dists[thisIter].size(); j++) {
                stdDev[j] += ((dists[thisIter][j] - averages[j]) * (dists[thisIter][j] - averages[j]));
            }
        }
        for (int i = 0; i < stdDev.size(); i++) {
            stdDev[i] /= (double) dists.size();
            stdDev[i] = sqrt(stdDev[i]);
        }

        return stdDev;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getStandardDeviation");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<seqDist> > Utils::getAverages(vector< vector< vector<seqDist> > >& calcDistsTotals, string mode) {
    try{

        vector< vector<seqDist>  > calcAverages; //calcAverages.resize(calcDistsTotals[0].size());
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            //calcAverages[i].resize(calcDistsTotals[0][i].size());
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            calcAverages.push_back(temp);
        }

        if (mode == "average") {
            for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //sum all groups dists for each calculator
                for (int i = 0; i < calcAverages.size(); i++) {  //initialize sums to zero.
                    for (int j = 0; j < calcAverages[i].size(); j++) {
                        calcAverages[i][j].dist += calcDistsTotals[thisIter][i][j].dist;
                    }
                }
            }

            for (int i = 0; i < calcAverages.size(); i++) {  //finds average.
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    calcAverages[i][j].dist /= (float) calcDistsTotals.size();
                }
            }
        }else { //find median
            for (int i = 0; i < calcAverages.size(); i++) { //for each calc
                for (int j = 0; j < calcAverages[i].size(); j++) {  //for each comparison
                    vector<double> dists;
                    for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //for each subsample
                        dists.push_back(calcDistsTotals[thisIter][i][j].dist);
                    }
                    sort(dists.begin(), dists.end());
                    calcAverages[i][j].dist = dists[(calcDistsTotals.size()/2)];
                }
            }
        }

        return calcAverages;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverages");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<seqDist> > Utils::getAverages(vector< vector< vector<seqDist> > >& calcDistsTotals) {
    try{

        vector< vector<seqDist>  > calcAverages; //calcAverages.resize(calcDistsTotals[0].size());
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            //calcAverages[i].resize(calcDistsTotals[0][i].size());
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            calcAverages.push_back(temp);
        }


        for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //sum all groups dists for each calculator
            for (int i = 0; i < calcAverages.size(); i++) {  //initialize sums to zero.
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    calcAverages[i][j].dist += calcDistsTotals[thisIter][i][j].dist;
                }
            }
        }

        for (int i = 0; i < calcAverages.size(); i++) {  //finds average.
            for (int j = 0; j < calcAverages[i].size(); j++) {
                calcAverages[i][j].dist /= (float) calcDistsTotals.size();
            }
        }

        return calcAverages;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverages");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<seqDist> > Utils::getStandardDeviation(vector< vector< vector<seqDist> > >& calcDistsTotals) {
    try{

        vector< vector<seqDist> > calcAverages = getAverages(calcDistsTotals);

        //find standard deviation
        vector< vector<seqDist>  > stdDev;
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            stdDev.push_back(temp);
        }

        for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int i = 0; i < stdDev.size(); i++) {
                for (int j = 0; j < stdDev[i].size(); j++) {
                    stdDev[i][j].dist += ((calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist) * (calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist));
                }
            }
        }

        for (int i = 0; i < stdDev.size(); i++) {  //finds average.
            for (int j = 0; j < stdDev[i].size(); j++) {
                stdDev[i][j].dist /= (float) calcDistsTotals.size();
                stdDev[i][j].dist = sqrt(stdDev[i][j].dist);
            }
        }

        return stdDev;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverages");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<seqDist> > Utils::getStandardDeviation(vector< vector< vector<seqDist> > >& calcDistsTotals, vector< vector<seqDist> >& calcAverages) {
    try{
        //find standard deviation
        vector< vector<seqDist>  > stdDev;
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            stdDev.push_back(temp);
        }

        for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int i = 0; i < stdDev.size(); i++) {
                for (int j = 0; j < stdDev[i].size(); j++) {
                    stdDev[i][j].dist += ((calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist) * (calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist));
                }
            }
        }

        for (int i = 0; i < stdDev.size(); i++) {  //finds average.
            for (int j = 0; j < stdDev[i].size(); j++) {
                stdDev[i][j].dist /= (float) calcDistsTotals.size();
                stdDev[i][j].dist = sqrt(stdDev[i][j].dist);
            }
        }

        return stdDev;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getAverages");
        exit(1);
    }
}

/**************************************************************************************************/
bool Utils::isContainingOnlyDigits(string input) {
    try{

        //are you a digit in ascii code
        for (int i = 0;i < input.length(); i++){
            if( input[i]>47 && input[i]<58){}
            else { return false; }
        }

        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "isContainingOnlyDigits");
        exit(1);
    }
}
/**************************************************************************************************/
/*M02352_41_000000000-AT06G_1_2104_18738_21630 Eukaryota(100);Archaeplastida(100);Chloroplastida(100);Chlorophyta(100);Mamiellophyceae(100);Mamiellales(100);Ostreococcus(100);Ostreococcus tauri(100);

 When I run remove.lineage with:
 taxon=Chloroplast-Mitochondria-unknown-Bacteria-Archaea-Metazoa-Charophyta

 The word "Chloroplast" in the taxon string gets matched to the lineage Chloroplastida in the taxonomy (above) and wipes out all of the green algae.*/

bool Utils::findTaxon(string tax, string stax) {
    try {
        string taxon = "";
        //int taxLength = tax.length();
        string searchTax = stax;
        if (searchTax[searchTax.length()-1] == ';') { searchTax = stax.substr(0, searchTax.length()-1); }

        int pos = tax.find(searchTax);
        if (pos != string::npos) {  //we found a match, but is it complete
            int endOfTaxon = pos + searchTax.length();
            if (tax[endOfTaxon] == ';') {  //we found a complete match
                return true;
            }
        }
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "findTaxon");
        exit(1);
    }
}
/**************************************************************************************************/
bool Utils::searchTax(string noQuotesTax, vector<string> listOfTaxons, vector<bool> taxonsHasConfidence, vector<string> noConfidenceTaxons, vector< vector< map<string, float> > > searchTaxons) {
    try {
        for (int j = 0; j < listOfTaxons.size(); j++) {

            string newtax = noQuotesTax;

            //if the users file contains confidence scores we want to ignore them when searching for the taxons, unless the taxon has them
            if (!taxonsHasConfidence[j]) {
                int hasConfidences = noQuotesTax.find_first_of('(');
                if (hasConfidences != string::npos) { newtax = noQuotesTax; removeConfidences(newtax); }

                if (findTaxon(newtax, noConfidenceTaxons[j])) { //this sequence contains the taxon the user wants
                    return true;  //since you belong to at least one of the taxons we want you are included so no need to search for other
                }
            }else{//if listOfTaxons[i] has them and you don't them remove taxons
                int hasConfidences = noQuotesTax.find_first_of('(');
                if (hasConfidences == string::npos) {

                    if (findTaxon(newtax, noConfidenceTaxons[j])) { //this sequence contains the taxon the user wants
                        return true;  //since you belong to at least one of the taxons we want you are included so no need to search for other
                    }
                }else { //both have confidences so we want to make sure the users confidences are greater then or equal to the taxons
                    //first remove confidences from both and see if the taxonomy exists

                    string noNewTax = noQuotesTax;
                    int hasConfidences = noQuotesTax.find_first_of('(');
                    if (hasConfidences != string::npos) { noNewTax = noQuotesTax; removeConfidences(noNewTax); }

                    if (findTaxon(noNewTax, noConfidenceTaxons[j])) { //if yes, then are the confidences okay

                        bool good = true; bool hasCon;
                        vector< map<string, float> > usersTaxon = getTaxons(newtax, hasCon);

                        //the usersTaxon is most likely longer than the searchTaxons, and searchTaxon[0] may relate to userTaxon[4]
                        //we want to "line them up", so we will find the the index where the searchstring starts
                        int index = 0;
                        for (int i = 0; i < usersTaxon.size(); i++) {

                            if (usersTaxon[i].begin()->first == searchTaxons[j][0].begin()->first) {
                                index = i;
                                int spot = 0;
                                bool goodspot = true;
                                //is this really the start, or are we dealing with a taxon of the same name?
                                while ((spot < searchTaxons[j].size()) && ((i+spot) < usersTaxon.size())) {
                                    if (usersTaxon[i+spot].begin()->first != searchTaxons[j][spot].begin()->first) { goodspot = false; break; }
                                    else { spot++; }
                                }

                                if (goodspot) { break; }
                            }
                        }

                        for (int i = 0; i < searchTaxons[j].size(); i++) {

                            if ((i+index) < usersTaxon.size()) { //just in case, should never be false
                                if (usersTaxon[i+index].begin()->second < searchTaxons[j][i].begin()->second) { //is the users cutoff less than the search taxons
                                    good = false;
                                    break;
                                }
                            }else { good = false; break; }
                        }

                        //passed the test so add you
                        if (good) { return true; }
                    }
                }
            }
        }

        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "searchTax");
        exit(1);
    }
}

/**************************************************************************************************/
vector< map<string, float> > Utils::getTaxons(string tax, bool& hasConfidence) {
    try {

        vector< map<string, float> > t;
        string taxon = "";
        int taxLength = tax.length();

        for(int i=0;i<taxLength;i++){
            if(tax[i] == ';'){
                string newtaxon = taxon; float confidence = 0;
                hasConfidence = hasConfidenceScore(newtaxon, confidence);

                map<string, float> temp; temp[newtaxon] = confidence; t.push_back(temp);
                taxon = "";
            }
            else{ taxon += tax[i]; }
        }

        return t;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getTaxons");
        exit(1);
    }
}
/**************************************************************************************************/
bool Utils::hasConfidenceScore(string& taxon, float& confidence) {
    try {
        int openParen = taxon.find_last_of('(');
        int closeParen = taxon.find_last_of(')');

        if ((openParen != string::npos) && (closeParen != string::npos)) {
            string confidenceScore = taxon.substr(openParen+1, (closeParen-(openParen+1)));
            if (isPositiveNumeric(confidenceScore)) {  //its a confidence
                taxon = taxon.substr(0, openParen); //rip off confidence
                mothurConvert(confidenceScore, confidence);
                return true;
            }
            else { confidence = 0; }//its part of the taxon
        }
        else{ confidence = 0; }

        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "hasConfidenceScore");
        exit(1);
    }
}
/**************************************************************************************************/
float Utils::removeConfidences(string& tax) {
    try {
        string temp = tax; float dummy; if (!hasConfidenceScore(temp, dummy)) { return 0; }

        string taxon;
        string newTax = "";
        string confidenceScore = "0";

        //remove last ";"
        if (tax.length() > 1) { tax = tax.substr(0, tax.length()-1); }
        vector<string> taxons; splitAtChar(tax, taxons, ';');

        for (int i = 0; i < taxons.size(); i++) {

            if (m->getControl_pressed()) { return 0; }

            taxon = taxons[i];

            int pos = taxon.find_last_of('(');
            if (pos != -1) {
                //is it a number?
                int pos2 = taxon.find_last_of(')');
                if (pos2 != -1) {
                    string temp = taxon.substr(pos+1, (pos2-(pos+1)));
                    if (isPositiveNumeric(temp)) {
                        taxon = taxon.substr(0, pos); //rip off confidence
                        confidenceScore = temp;
                    }
                }
            }
            taxon += ";";

            newTax += taxon;
        }

        tax = newTax;

        float confidence = 0; mothurConvert(confidenceScore, confidence);

        return confidence;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "removeConfidences");
        exit(1);
    }
}
/**************************************************************************************************/
string Utils::removeQuotes(string tax) {
    try {

        string taxon;
        string newTax = "";

        for (int i = 0; i < tax.length(); i++) {

            if (m->getControl_pressed()) { return newTax; }

            if ((tax[i] != '\'') && (tax[i] != '\"')) { newTax += tax[i]; }

        }

        return newTax;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "removeQuotes");
        exit(1);
    }
}
/**************************************************************************************************/
// function for calculating standard deviation
double Utils::getStandardDeviation(vector<int>& featureVector){
    try {
        //finds sum
        double average = 0;
        for (int i = 0; i < featureVector.size(); i++) { average += featureVector[i]; }
        average /= (double) featureVector.size();

        //find standard deviation
        double stdDev = 0;
        for (int i = 0; i < featureVector.size(); i++) { //compute the difference of each dist from the mean, and square the result of each
            stdDev += ((featureVector[i] - average) * (featureVector[i] - average));
        }

        stdDev /= (double) featureVector.size();
        stdDev = sqrt(stdDev);

        return stdDev;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "getStandardDeviation");
        exit(1);
    }
}
/*****************************************************************/
//this code is a mess and should be rethought...-slw
vector<string> Utils::parseTreeFile(string filename) {

    //only takes names from the first tree and assumes that all trees use the same names.
    try {
        //string filename = current->getTreeFile();
        ifstream filehandle;
        Utils util; util.openInputFile(filename, filehandle);
        int comment;
        char c;
        comment = 0;

        vector<string> Treenames;
        if((c = filehandle.peek()) != '#') {  //ifyou are not a nexus file

            while ((c = filehandle.peek()) != ';') {
                if (m->getControl_pressed()) {  filehandle.close(); return Treenames; }
                // get past comments
                if(c == '[')    { comment = 1; }
                if(c == ']')    { comment = 0; }
                if((c == '(') && (comment != 1)){ break; }
                filehandle.get();
            }

            Treenames = readTreeString(filehandle);

        }else if((c = filehandle.peek()) == '#') { //ifyou are a nexus file
            string holder = "";

            // get past comments
            while(holder != "translate" && holder != "Translate"){
                if (m->getControl_pressed()) {  filehandle.close(); return Treenames; }
                if(holder == "[" || holder == "[!") { comment = 1; }
                if(holder == "]")                   { comment = 0; }
                filehandle >> holder;

                //if there is no translate then you must read tree string otherwise use translate to get names
                if((holder == "tree") && (comment != 1)){
                    //pass over the "tree rep.6878900 = "
                    while (((c = filehandle.get()) != '(') && ((c = filehandle.peek()) != EOF)) {;}

                    if(c == EOF) { break; }
                    filehandle.putback(c);  //put back first ( of tree.
                    Treenames = readTreeString(filehandle);

                    break;
                }

                if (filehandle.eof()) { break; }
            }

            //use nexus translation rather than parsing tree to save time
            if((holder == "translate") || (holder == "Translate")) {

                string number, name, h;
                h = ""; // so it enters the loop the first time
                while((h != ";") && (number != ";")) {
                    if (m->getControl_pressed()) {  filehandle.close(); return Treenames; }
                    filehandle >> number;
                    filehandle >> name;

                    //c = , until done with translation then c = ;
                    h = name.substr(name.length()-1, name.length());
                    name.erase(name.end()-1);  //erase the comma
                    Treenames.push_back(number);
                }
                if(number == ";") { Treenames.pop_back(); }  //in case ';' from translation is on next line instead of next to last name
            }
        }
        filehandle.close();

        return Treenames;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "parseTreeFile");
        exit(1);
    }
}
/*******************************************************/
vector<string> Utils::readTreeString(ifstream& filehandle)	{
    try {
        char c;
        string name;  //, k
        vector<string> Treenames;

        while((c = filehandle.peek()) != ';') {
            if (m->getControl_pressed()) {  return Treenames; }

            if(c == ')')  {
                //to pass over labels in trees
                c=filehandle.get();
                while((c!=',') && (c != -1) && (c!= ':') && (c!=';')){ c=filehandle.get(); }
                filehandle.putback(c);
            }
            if(c == ';') { return Treenames; }
            if(c == -1) { return Treenames; }
            //if you are a name
            if((c != '(') && (c != ')') && (c != ',') && (c != ':') && (c != '\n') && (c != '\t') && (c != 32)) { //32 is space
                name = "";
                c = filehandle.get();
                while ((c != '(') && (c != ')') && (c != ',') && (c != ':')  && (c != '\n') && (c != 32) && (c != '\t')) {
                    name += c;
                    c = filehandle.get();
                }

                if (name != "\r" ) { Treenames.push_back(name);   }
                filehandle.putback(c);
            }

            if(c  == ':') { //read until you reach the end of the branch length
                while ((c != '(') && (c != ')') && (c != ',') && (c != ';') && (c != '\n') && (c != '\t') && (c != 32)) { c = filehandle.get(); }
                filehandle.putback(c);
            }
            c = filehandle.get();
            if(c == ';') { return Treenames; }
            if(c == ')') { filehandle.putback(c); }
            if (filehandle.eof()) { break; }
        }
        return Treenames;
    }
    catch(exception& e) {
        m->errorOut(e, "Utils", "readTreeString");
        exit(1);
    }
}
/*********************************************************************************************/
