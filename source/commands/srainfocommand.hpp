//
//  srainfocommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/29/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef srainfocommand_hpp
#define srainfocommand_hpp

#include "command.hpp"

/**************************************************************************************************/

class SRAInfoCommand : public Command {
public:
    SRAInfoCommand(string);
    SRAInfoCommand();
    ~SRAInfoCommand(){}
    
    vector<string> setParameters();
    string getCommandName()            { return "sra.info";             }
    string getCommandCategory()        { return "Sequence Processing";  }
    
    string getOutputPattern(string);
    
    string getHelpString();
    string getCitation()    { return ".... Add reference for NCBI .... http://www.mothur.org/wiki/sra.info"; }
    string getDescription() { return "extracts fastq or sff files from sra file using fasterq_dump or sff_dump programs written by NCBI"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    
    bool abort;
    vector<string> outputNames;
    string srafile, outputDir, outputType, fasterQLocation;
    int processors;
    
    void runFastqDump();
};

/**************************************************************************************************/




#endif /* srainfocommand_hpp */
