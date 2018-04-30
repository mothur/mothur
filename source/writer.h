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

/***********************************************************************/

class OutputWriter {
public:

    OutputWriter (std::shared_ptr<SynchronizedOutputFile> s) : sf(s) {}
    
    void write (const string& dataToWrite) { sf->write(dataToWrite); }
private:
    std::shared_ptr<SynchronizedOutputFile> sf;
};
/***********************************************************************/

#endif 
