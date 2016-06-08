//
//  distcdataset.h
//  Mothur
//
//  Created by Sarah Westcott on 6/8/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__distcdataset__
#define __Mothur__distcdataset__

#include "mothurout.h"

class DistCDataSet  {
    
public:
    
    DistCDataSet();
    ~DistCDataSet() {}
    string getColumnFile() { return writeColumnFile(); }
    string getCountFile()  { return writeCountFile();  }
    
private:
    MothurOut* m;
    string writeColumnFile();
    string writeCountFile();
    
};


#endif /* defined(__Mothur__distcdataset__) */
