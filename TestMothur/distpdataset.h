//
//  distdataset.h
//  Mothur
//
//  Created by Sarah Westcott on 6/6/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__distdataset__
#define __Mothur__distdataset__

#include "mothurout.h"

class DistPDataSet  {
    
public:
    
    DistPDataSet();
    ~DistPDataSet() {}
    
    string getPhylipFile() { return writePhylipFile(); }
    
private:
    MothurOut* m;
    string writePhylipFile();
    
};


#endif /* defined(__Mothur__distdataset__) */
