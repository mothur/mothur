//
//  interactengine.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef interactengine_hpp
#define interactengine_hpp

#include "engine.hpp"

class InteractEngine : public Engine {
public:
    InteractEngine(string);
    ~InteractEngine();
    virtual bool getInput();
    
private:
    string getCommand();
    
};


#endif /* interactengine_hpp */
