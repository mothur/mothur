//
//  distdataset.h
//  Mothur
//
//  Created by Sarah Westcott on 6/6/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__distdataset__
#define __Mothur__distdataset__

#include "fastqread.h"

class DistDataSet  {
    
public:
    
    DistDataSet();
    ~DistDataSet() {}
    string getColumnFile() { return writeColumnFile(); }
    string getPhylipFile() { return writePhylipFile(); }
    string getNameFile()   { return writeNameFile();   }
    string getCountFile()  { return writeCountFile();  }
    
private:
    MothurOut* m;
    string writeColumnFile();
    string writePhylipFile();
    string writeNameFile();
    string writeCountFile();
    
};


#endif /* defined(__Mothur__distdataset__) */
