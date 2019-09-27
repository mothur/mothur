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
 
 file option 4 - only vaild if mode is set to parsefastqpacbio
 
 group1 pacBiofastqfile1
 group2 pacBiofastqfile2
 ...
 
 file option 5
 
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 ...
 
 file option 6
 
 My.forward.fastq My.reverse.fastq none My.rindex.fastq //none is an option is no forward or reverse index file
 ...
 
 
 
 ********* fileOption; //1 -> 2 column(4 forms of 2 column), 2 -> 3 column, 3 -> 4 column ******************
 */
/**************************************************************************************************/

class FileFile {
    
public:
    FileFile(string, string); //provide file file and read file, mode (ie. contigs) mode options include: mimarks,contigs,parseFastq,renameSeqs,sra,parsefastqpacbio
    ~FileFile() {}
    
    vector< vector<string> > getFiles() { return files; }
    
    bool isColumnWithGroupNames() { return columnWithGroups; }
    int getFileFormat() { return fileOption; }
    
    bool isGZ() { return gz; } //are files listed in file compressed
    bool containsIndexFiles() { return hasIndex; } //indicates oligos file is required
    map<int, string> getFile2Group() { return file2Group; } //fileIndex2GroupName, files[0]'s group is -> file2Group[0]
    
protected:
    
    MothurOut* m;
    CurrentFile* current;
    Utils util;
    string filename, mode, inputDir, mpath;
    bool gz, hasIndex, columnWithGroups;
    int fileOption; //1 -> 2 column(3 forms of 2 column), 2 -> 3 column, 3 -> 4 column
    vector< vector<string> > files;
    map<int, string> file2Group;
    
    vector< vector<string> > read(string, string); //read file, used with () constructor
    bool validateFiles(vector<string> pieces, string& forward, string& reverse, string& findex, string& rindex, string& group); //checks locations, abletoOPen, fileOPtion
    void setGZ(string forward, string reverse, string findex, string rindex, bool&, bool&);
};

/**************************************************************************************************/



#endif /* filefile_hpp */
