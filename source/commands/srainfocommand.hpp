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
    string getCitation()    { return "Wrapper for prefetch and fasterq_dump programs written by NCBI https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software http://www.mothur.org/wiki/sra.info"; }
    string getDescription() { return "extracts fastq files from samples using prefetch and fasterq_dump program written by NCBI"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    
    bool abort, compressGZ;
    vector<string> outputNames;
    string accnosfile, outputDir, outputType, fasterQLocation, prefetchLocation;
    int processors, maxSize;
    
    string runPreFetch(string);
    bool runFastqDump(string, vector<string>&);
    void runSystemCommand(string);
    bool checkVersion(string versionNeeded, string versionProvided);

};

/**************************************************************************************************/




#endif /* srainfocommand_hpp */
