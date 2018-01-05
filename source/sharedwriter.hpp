//
//  SharedWriter.hpp
//  Mothur
//
//  Created by Sarah Westcott on 12/7/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef SharedWriter_hpp
#define SharedWriter_hpp

#include "mothurout.h"
#include "utils.hpp"
#include "sequence.hpp"

/***********************************************************************/
class SynchronizedOutputFile {
public:
    SynchronizedOutputFile (const string& p) : path(p) { util.openOutputFile(p, out); }
    
    void write (const string& dataToWrite) {
        std::lock_guard<std::mutex> lock((writerMutex)); // Ensure that only one thread can execute at a time
        out << dataToWrite;
    }
    
    void setFixedShowPoint() {  out.setf(ios::fixed, ios::showpoint);  }
    void setPrecision(int p)  { out << setprecision(p);  }
    
private:
    string path;
    std::mutex writerMutex;
    Utils util;
    ofstream out;
};

/***********************************************************************
class SynchronizedInputFile {
public:
    SynchronizedInputFile (const string& p) : path(p) { util.openInputFile(p, in); }
    
    void read (string& dataToWrite) {
        std::lock_guard<std::mutex> lock((readerMutex)); // Ensure that only one thread can execute at a time
        in >> dataToWrite;
    }
    
    void read (Sequence& seq) {
        std::lock_guard<std::mutex> lock((readerMutex)); // Ensure that only one thread can execute at a time
        Sequence temp(in); util.gobble(in);
        seq = temp;
    }
    
    bool endOfFile() { return (in.eof()); }
    
private:
    string path;
    std::mutex readerMutex;
    Utils util;
    ifstream in;
};

/***********************************************************************/

#endif /* SharedWriter_hpp */
