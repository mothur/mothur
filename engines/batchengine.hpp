//
//  batchengine.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef batchengine_hpp
#define batchengine_hpp

#include "engine.hpp"

class BatchEngine : public Engine {
public:
    
    BatchEngine(string, string);
    ~BatchEngine();
    
    virtual bool getInput();
    
private:
    ifstream inputBatchFile;
    string getNextCommand(ifstream&);
    string batchFileName;
    bool openedBatch;

};

#endif /* batchengine_hpp */
