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
    
    BatchEngine(string, string, map<string, string>);
    ~BatchEngine();
    
    virtual bool getInput();
    bool getOpenedBatch() { return openedBatch; }
    
private:
    ifstream inputBatchFile;
    string getNextCommand(ifstream&);
    string findType(string);
    string batchFileName;
    bool openedBatch;
    time_t bstart;
    int numBatches;
    

};

#endif /* batchengine_hpp */
