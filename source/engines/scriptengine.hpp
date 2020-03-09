//
//  scriptengine.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/21/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef scriptengine_hpp
#define scriptengine_hpp

#include "engine.hpp"

class ScriptEngine : public Engine {
public:
    ScriptEngine(string, string, map<string, string>);
    ~ScriptEngine();
    virtual bool getInput();
    bool openedBatch;
private:
    string listOfCommands;
    string getNextCommand(string&);

};


#endif /* scriptengine_hpp */
