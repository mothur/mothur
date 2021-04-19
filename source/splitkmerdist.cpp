//
//  splitkmerdist.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/30/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "splitkmerdist.hpp"
#include "kmerdist.hpp"
#include "sequence.hpp"
#include "kmer.hpp"

/***********************************************************************/

SplitKmerDistance::SplitKmerDistance(string f, string o, double c, int ksize){
    try{
        m = MothurOut::getInstance();
    
        cutoff = c;
        fastafile = f;
        kmerSize = ksize;
        string outputDir = o;
        if (outputDir == "") {  outputDir += util.hasPath(fastafile);  }
        
        //reads fasta file and fills kmerDB
        fillKmerDB();
        
        vector<vector<long long> > groups = split();
        
        for (int i = 0; i < groups.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            string outputFileName = outputDir + util.getRootName(util.getSimpleName(fastafile)) + toString(i) + ".temp";
            
            //printGroup(seqDB, groups[i], outputFileName);
            parsedFiles.push_back(outputFileName);
        }
    
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "SplitKmerDistance");
        exit(1);
    }
}
/***********************************************************************/
vector<vector<long long> > SplitKmerDistance::split(){
    try {
        cutoff = 85;
        int numGroups = 0;
        
        //seqsseqsGroupAssignment[0] = group assignment for seq 0 [db->get(0)].
        vector<int> seqsGroupAssignment(numSeqs, -1); //assign all seqs no group
        
        //when we merge groups, rather than reassigning all the seqs in that group, let's reassign the group
        //mergedGroups[0] = group assignment for group 0.
        vector<int> mergedGroups;
        
        KmerDist distCalculator(kmerSize);
        
        map<int, long long> distanceSummary; //kmerdist rounded to 2 places, count of that dist in matrix
        map<int, long long>::iterator itDist;
        
        for(long long i = 0; i < numSeqs; i++){
        
            vector<int> seqA = getUniqueKmers(i);
            
            for(long long j = 0; j < i; j++){
                
                if (m->getControl_pressed()) {  break;  }
                
                vector<bool> seqB = kmerDB[j];
                
                int length = min(lengths[i], lengths[j]);
                   
                double dist = distCalculator.calcDist(seqA, seqB, length);
                
                int value = (int)(dist * 100 + .5);
                dist = value;
                
                itDist = distanceSummary.find(dist);
                if (itDist == distanceSummary.end()) {
                    distanceSummary[dist] = 1;
                }else {
                    distanceSummary[dist]++;
                }
                
                if(dist <= cutoff){
                    
                    int groupIDA = findRootGroup(mergedGroups, seqsGroupAssignment[i]);
                    int groupIDB = findRootGroup(mergedGroups, seqsGroupAssignment[i]);
                    int groupID = -1;
                    
                    if(groupIDA != -1 && groupIDB != -1){ //both are already assigned to a group, so merge the groups. set to groupIDA
                        int thisGroup = mergedGroups[groupIDA];
                        groupID = thisGroup;
                        mergedGroups[groupIDB] = thisGroup;
                    
                    }else if(groupIDA != -1 && groupIDB == -1){ //seqA has a group and seqB is unassigned
                        groupID = mergedGroups[groupIDA]; //assign seqB to seqA's group
                        
                    }else if(groupIDA == -1 && groupIDB != -1) {  //seqB has a group and seqA is unassigned
                        
                        groupID = mergedGroups[groupIDB]; //assign seqA to seqB's group
                        
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
        
        vector<long long> bracketTotals; bracketTotals.resize(22, 0);
        for (itDist = distanceSummary.begin(); itDist != distanceSummary.end(); itDist++) {
            cout << (itDist->first / 100.0) << '\t' << itDist->second << endl;
            int index = 0;
            if (itDist->first <= 10) { index = 0; }
            else if ((itDist->first > 10) && (itDist->first <= 20)) { index = 1; }
            else if ((itDist->first > 20) && (itDist->first <= 30)) { index = 2; }
            else if ((itDist->first > 30) && (itDist->first <= 40)) { index = 3; }
            else if ((itDist->first > 40) && (itDist->first <= 50)) { index = 4; }
            else if ((itDist->first > 50) && (itDist->first <= 60)) { index = 5; }
            else if ((itDist->first > 60) && (itDist->first <= 70)) { index = 6; }
            else if ((itDist->first > 70) && (itDist->first <= 80)) { index = 7; }
            else if ((itDist->first > 80) && (itDist->first <= 90)) { index = 8; }
            else if ((itDist->first > 90) && (itDist->first <= 100)) { index = 9; }
            else if ((itDist->first > 100) && (itDist->first <= 110)) { index = 10; }
            else if ((itDist->first > 110) && (itDist->first <= 120)) { index = 11; }
            else if ((itDist->first > 120) && (itDist->first <= 130)) { index = 12; }
            else if ((itDist->first > 130) && (itDist->first <= 140)) { index = 13; }
            else if ((itDist->first > 140) && (itDist->first <= 150)) { index = 14; }
            else if ((itDist->first > 150) && (itDist->first <= 160)) { index = 15; }
            else if ((itDist->first > 160) && (itDist->first <= 170)) { index = 16; }
            else if ((itDist->first > 170) && (itDist->first <= 180)) { index = 17; }
            else if ((itDist->first > 180) && (itDist->first <= 190)) { index = 18; }
            else if ((itDist->first > 190) && (itDist->first <= 200)) { index = 19; }
            else if ((itDist->first > 200) && (itDist->first <= 210)) { index = 20; }
            else if ((itDist->first > 210) && (itDist->first <= 220)) { index = 21; }
            
           
            bracketTotals[index] += itDist->second;
        }
        
        for (int i = 0; i < bracketTotals.size(); i++) {
            cout << (i+1)*10 << '\t' << bracketTotals[i] << endl;
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
        
        for (int i = 0; i < groups.size(); i++) {
            cout << i << '\t' << groups[i].size() << endl;
        }
        
        return groups;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "split");
        exit(1);
    }
}
/***********************************************************************/
int SplitKmerDistance::fillKmerDB(){
    try {
        
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        
        int maxKmer = power4s[kmerSize];
        
        Kmer kmer(kmerSize);
        
        ifstream inFASTA; util.openInputFile(fastafile, inFASTA);
        numSeqs = 0;
        
        while(!inFASTA.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(inFASTA); util.gobble(inFASTA);
            
            vector<bool> kmerLocations; kmerLocations.resize(maxKmer+1, false);
            
            int numKmers = seq.getNumBases() - kmerSize + 1;
        
            for(int i=0;i<numKmers;i++){
                int kmerNumber = kmer.getKmerNumber(seq.getUnaligned(), i);
                
                kmerLocations[kmerNumber] = true; //ok, we've seen the kmer now
            }
            
            kmerDB.push_back(kmerLocations);
            lengths.push_back(seq.getNumBases());
            numSeqs++;
        }
    
        inFASTA.close();
        
        return kmerDB.size();
  
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "fillKmerDB");
        exit(1);
    }
}
/***********************************************************************/
vector<int> SplitKmerDistance::getUniqueKmers(int i){
    try {
        
        vector<bool> seqsKmers = kmerDB[i];
        vector<int> uniques;
        
        for (int k = 0; k < seqsKmers.size(); k++) {
            if (seqsKmers[k]) {
                uniques.push_back(k);
            }
        }
        
        return uniques;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitKmerDistance", "getUniqueKmers");
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
/***********************************************************************
void SplitKmerDistance::printGroup(vector<long long> group, string outputFileName){
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
