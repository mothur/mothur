//
//  splitkmerdist.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/30/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "splitkmerdist.hpp"
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
        
        //seqsseqsGroupAssignment[0] = group assignment for seq 0 [db->get(0)].
        vector<int> seqsGroupAssignment(numSeqs, -1); //assign all seqs no group
        
        //when we merge groups, rather than reassigning all the seqs in that group, let's reassign the group
        //mergedGroups[0] = group assignment for group 0.
        vector<int> mergedGroups;
        
        KmerDist distCalculator(cutoff, kmerSize);
        
        map<int, long long> distanceSummary; //kmerdist rounded to 2 places, count of that dist in matrix
        map<int, long long>::iterator itDist;
        
        for(long long i = 0; i < numSeqs; i++){
            
            Sequence seq = db->getSeq(i);
            
            for(long long j = 0; j < i; j++){
                
                if (m->getControl_pressed()) {  break;  }
                
                Sequence seqJ = db->getSeq(j);
                   
                double dist = distCalculator.calcDist(seq, seqJ);
                
                int value = (int)(dist * 100 + .5);
                dist = value;
                
                itDist = distanceSummary.find(dist);
                if (itDist == distanceSummary.end()) {
                    distanceSummary[dist] = 1;
                }else {
                    distanceSummary[dist]++;
                }
                
               /* if(dist <= cutoff){
                    
                    int groupIDA = findRootGroup(mergedGroups, seqsGroupAssignment[i]);
                    int groupIDB = findRootGroup(mergedGroups, seqsGroupAssignment[i]);
                    int groupID = -1;
                    
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
               */ //}//end if cutoff
            }//end j
        }//end i
        
        
        for (itDist = distanceSummary.begin(); itDist != distanceSummary.end(); itDist++) {
            cout << (itDist->first / 100) << '\t' << itDist->second << endl;
        }
        
        vector<vector<long long> > groups;
        map<int, int> groupIndex;
        
        for (long long i = 0; i < seqsGroupAssignment.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            if (seqsGroupAssignment[i] != -1) { //you have a distance below the cutoff
                int group = findRootGroup(mergedGroups, seqsGroupAssignment[i]);
                
                map<int, int>::iterator it = groupIndex.find(group);
                if (it != groupIndex.end()) {
                    groups[it->second].push_back(i);
                }else{
                    vector<long long> temp; temp.push_back(i);
                    groupIndex[group] = groups.size();
                    groups.push_back(temp);
                }
                
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
int SplitKmerDistance::findRootGroup(vector<int>& mergedGroups, int pos){
    try {
        
        //if unassigned
        if (pos == -1) { return pos; }
        
        //if the mergedGroups[10] = 10 then you are at the root
        //if mergedGroups[10] != 10, find parent group
        //mergedGroups[10] = 5, then check mergeGroups[5].
        
        int rootGroup = mergedGroups[pos];
        
        if (rootGroup != pos) { //you need to look at your parent
            rootGroup = findRootGroup(mergedGroups, rootGroup);
            mergedGroups[pos] = rootGroup;
        }//else you are the root
        
        return rootGroup;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "findRootGroup");
        exit(1);
    }
}
/***********************************************************************/
void SplitKmerDistance::printGroup(SequenceDB*& db, vector<long long> group, string outputFileName){
    try {
        
        ofstream out; util.openOutputFile(outputFileName, out);
        
        for (int i = 0; i < group.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            Sequence seq = db->getSeq(group[i]);
            
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
