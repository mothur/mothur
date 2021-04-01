//
//  splitkmerdist.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/30/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "splitkmerdist.hpp"
#include "getseqscommand.h"
#include "kmerdist.hpp"

/***********************************************************************/

SplitKmerDistance::SplitKmerDistance(string f, string o, double c, int ksize){
    try{
        m = MothurOut::getInstance();
    
        cutoff = c;
        fastafile = f;
        kmerSize = ksize;
        string outputDir = o;
        if (outputDir == "") {  outputDir += util.hasPath(fastafile);  }
        
        //read fasta file into seqDB
        ifstream inFASTA; util.openInputFile(fastafile, inFASTA);
        SequenceDB* seqDB = new SequenceDB(inFASTA); inFASTA.close();
        
        vector<vector<long long> > groups = split(seqDB);
        
        for (int i = 0; i < groups.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            string outputFileName = outputDir + util.getRootName(util.getSimpleName(fastafile)) + toString(i) + ".temp";
            
            printGroup(seqDB, groups[i], outputFileName);
            parsedFiles.push_back(outputFileName);
        }
    
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "SplitKmerDistance");
        exit(1);
    }
}
/***********************************************************************/
vector<vector<long long> > SplitKmerDistance::split(SequenceDB*& db){
    try {
        
        int numGroups = 0;

        long long numSeqs = db->getNumSeqs();
        
        vector<int> seqsGroupAssignment(numSeqs, -1); //assign all seqs no group
        
        //when we merge groups, rather than reassigning all the seqs in that group, let's reassign the group
        //mergedGroups[0] = group assignment for group 0.
        vector<int> mergedGroups;
        
        KmerDist distCalculator(cutoff, kmerSize);
        
        for(long long i = 0; i < numSeqs; i++){
            
            Sequence seq = db->get(i);
            
            for(long long j = 0; j < i; j++){
                
                if (m->getControl_pressed()) {  break;  }
                
                Sequence seqJ = db->get(j);
                   
                double dist = distCalculator.calcDist(seq, seqJ);
                
                if(dist <= cutoff){
                    
                    int groupIDA = seqsGroupAssignment[i]; int groupIDB = seqsGroupAssignment[i]; int groupID = -1;
                    
                    if(groupIDA != -1 && groupIDB != -1){ //both are already assigned to a group, so merge the groups. set to groupIDA
                        int thisGroup = mergedGroups[groupIDA];
                        groupID = thisGroup;
                        mergedGroups[groupIDB] = thisGroup;
                    
                    }else if(groupIDA != -1 && groupIDB == -1){ //seqA has a group and seqB is unassigned
                        groupID = mergedGroups[groupIDA];
                        
                    }else if(groupIDA == -1 && groupIDB != -1) {  //seqB has a group and seqA is unassigned
                        
                        groupID = mergedGroups[groupIDB];
                        
                    }else {  //both seqs have no group
                        groupID = numGroups; //we need a new group
                        mergedGroups.push_back(numGroups); //assign group to merge with self
                        numGroups++;
                    }
                    
                    seqsGroupAssignment[i] = groupID;
                    seqsGroupAssignment[j] = groupID;
                }//end if cutoff
            }//end j
        }//end i
        
        vector<vector<long long> > groups;
        for (int i = 0; i < seqsGroupAssignment.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            if (seqsGroupAssignment[i] != -1) { //you have a distance below the cutoff
                
                
                
            }
        }
        
        return groups;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "split");
        exit(1);
    }
}
/***********************************************************************/
void SplitKmerDistance::printGroup(SequenceDB*& db, vector<long long> group, string outputFileName){
    try {
        
        ofstream out; util.openOutputFile(outputFileName, out);
        
        for (int i = 0; i < group.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            Sequence seq = db->get(group[i]);
            
            seq.printSequence(out);
        }
        out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "printGroup");
        exit(1);
    }
}
/***********************************************************************/
