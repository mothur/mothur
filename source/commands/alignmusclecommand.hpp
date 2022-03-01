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
    string getCitation() { return "R.C. Edgar (2021) MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. (https://www.biorxiv.org/content/10.1101/2021.06.20.449169v1.full.pdf)\nhttp://www.mothur.org/wiki/align.muscle\n"; }
    string getDescription()        { return "align protein sequences"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, usePerturb, stratified, diversified, useReplicates, useConsiters, useRefineiters, usePerm;
    string fastafile, perturb, method, perm, replicates, consiters, refineiters, muscleLocation;
    int processors;
    vector<string> outputNames;
    
    void wrapperFunction();
    
};
/**************************************************************************************************/


#endif /* alignmusclecommand_hpp */
