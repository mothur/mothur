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
#include "utils.hpp"
#include "currentfile.h"

class DistCDataSet  {
    
public:
    
    DistCDataSet();
    ~DistCDataSet() {}
    string getColumnFile() { return columnFile; }
    vector<string> getFiles(int);
    string getCountFile()  { return countFile;  }
    
private:
    MothurOut* m;
    Utils util;
    string columnFile, countFile;
    CurrentFile* current;
    
};


#endif /* defined(__Mothur__distcdataset__) */
