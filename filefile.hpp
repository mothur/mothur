//
//  filefile.hpp
//  Mothur
//
//  Created by Sarah Westcott on 12/5/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef filefile_hpp
#define filefile_hpp

#include "mothurout.h"
#include "utils.hpp"
#include "currentfile.h"

/**************************************************************************************************/

class FileFile {
    
public:
    FileFile(string, string); //provide file file and read file, mode (ie. contigs, ...)
    FileFile();
    ~FileFile() {}
    
    vector< vector<string> > read(string, string); //read file, used with () constructor
    vector< vector<string> > getFiles() { return files; }
    bool getCreateGroup() { return createFileGroup; }
    bool isGZ() { return gz; }
    bool containsIndexFiles() { return hasIndex; }
    map<int, string> getFile2Group() { return file2Group; }
    
protected:
    
    MothurOut* m;
    CurrentFile* current;
    Utils util;
    string filename, mode, inputDir, mpath;
    bool createFileGroup, gz, hasIndex;
    
    vector< vector<string> > files;
    map<int, string> file2Group;
};

/**************************************************************************************************/



#endif /* filefile_hpp */
