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

/***********************************************************************/
class SynchronizedOutputFile {
public:
    SynchronizedOutputFile (const string& p)                : path(p) { util.openOutputFile(p, out);        }
    SynchronizedOutputFile (const string& p, bool append)   : path(p) { util.openOutputFileAppend(p, out);  }
    ~SynchronizedOutputFile() { if (out.is_open()) { out.close(); } } //if we forgot to close()
    
    void write (const string& dataToWrite) {
        std::lock_guard<std::mutex> lock((writerMutex)); // Ensure that only one thread can execute at a time
        out << dataToWrite;
    }
    void close() { if (out.is_open()) { out.close(); } }
    
    void setFixedShowPoint()    {  out.setf(ios::fixed, ios::showpoint);    }
    void setPrecision(int p)    {  out << setprecision(p);                  }
    
private:
    string path;
    std::mutex writerMutex;
    Utils util;
    ofstream out;
};

/***********************************************************************/

#endif
