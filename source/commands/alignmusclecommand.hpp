//
//  alignmusclecommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 2/16/22.
//  Copyright © 2022 Schloss Lab. All rights reserved.
//

#ifndef alignmusclecommand_hpp
#define alignmusclecommand_hpp

#include "command.hpp"

/***********************************************************/

class AlignMuscleCommand : public Command {
    
public:
    AlignMuscleCommand(string);
    ~AlignMuscleCommand() {}
    
    vector<string> setParameters();
    string getCommandName()            { return "align.muscle";        }
    string getCommandCategory()        { return "Sequence Processing"; }
    
    string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);
    string getCitation() { return "muscle Edgar, R.C. MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics 5, 113 (2004). https://doi.org/10.1186/1471-2105-5-113\nhttp://www.mothur.org/wiki/align.muscle\n"; }
    string getDescription()        { return "align protein sequences"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, usePerturb, stratified, diversified, useReplicates, useConsiters, useRefineiters, usePerm;
    string fastafile, perturb, method, perm, replicates, consiters, refineiters, muscleLocation;
    int processors;
    vector<string> outputNames;
    
    void driver();
    
};
/**************************************************************************************************/


#endif /* alignmusclecommand_hpp */
