//
//  writer.h
//  Mothur
//
//  Created by Sarah Westcott on 12/7/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef writer_h
#define writer_h

#include "sharedwriter.hpp"
#include "sequence.hpp"

/***********************************************************************/

class OutputWriter {
public:
    OutputWriter (std::shared_ptr<SynchronizedOutputFile> s) : sf(s) {}
    
    void write (const string& dataToWrite) { sf->write(dataToWrite); }
    void write (Sequence& dataToWrite) { sf->write(dataToWrite); }
private:
    std::shared_ptr<SynchronizedOutputFile> sf;
};
/***********************************************************************/

class InputReader {
public:
    InputReader (std::shared_ptr<SynchronizedInputFile> s) : sf(s) {}
    
    bool endOfFile() { return sf->endOfFile(); }
    void read (string& data) { sf->read(data); }
    void read (Sequence& data) { sf->read(data); }
private:
    std::shared_ptr<SynchronizedInputFile> sf;
};

/***********************************************************************/



#endif /* writer_h */
