//
//  pam.cpp
//  Mothur
//
//  Created by SarahsWork on 12/10/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "pam.h"
#define DBL_EPSILON 1e-9

/**************************************************************************************************/
Pam::Pam(vector<vector<int> > c, vector<vector<double> > d, int p) : CommunityTypeFinder() {
    try {
        countMatrix = c;
        numSamples = (int)d.size();
        numOTUs = (int)c[0].size();
        numPartitions = p;
        dists = d;
        
        largestDist = 0;
        for (int i = 0; i < dists.size(); i++) {
            for (int j = i; j < dists.size(); j++) {
                if (m->control_pressed) { break; }
                if (dists[i][j] > largestDist) { largestDist = dists[i][j]; } 
            }
        }
        
        buildPhase(); //choosing the medoids
        swapPhase(); //optimize clusters
    }
	catch(exception& e) {
		m->errorOut(e, "Pam", "Pam");
		exit(1);
	}
}
/**************************************************************************************************/
//sets Dp[0] does not set Dp[1]. chooses intial medoids.
int Pam::buildPhase() {
    try {
        
        if (m->debug) { m->mothurOut("[DEBUG]: building medoids\n"); }
        
        vector<double> gains; gains.resize(numSamples);
        
        largestDist *= 1.1 + 1; //make this distance larger than any distance in the matrix
        Dp.resize(numSamples);
        for (int i = 0; i < numSamples; i++) { Dp[i].push_back(largestDist); Dp[i].push_back(largestDist); } //2 smallest dists for this sample in this partition
        
        zMatrix.resize(numPartitions);
        for(int i=0;i<numPartitions;i++){
            zMatrix[i].assign(numSamples, 0);
        }
    
        for (int k = 0; k < numPartitions; k++) {
            
            int medoid = -1;
            double totalGain = 0.0;
            double clusterGain = 0.0;
            
            for (int i = 0; i < numSamples; i++) {  //does this need to be square?? can we do lt?
                if (m->control_pressed) { break; }
        
                if (medoids.count(i) == 0) { //is this sample is NOT a medoid?
                    gains[i] = 0.0;
                
                    for (int j = 0; j < numSamples; j++) {
                        //cout << i << '\t' << j << '\t' <<   Dp[i][0] << '\t' << dists[i][j] << '\t' << totalGain << endl;
                        totalGain = Dp[i][0] - dists[i][j];
                        if (totalGain > 0.0) { gains[i] += totalGain; }
                    }
                    if (m->debug) { m->mothurOut("[DEBUG]: " + toString(i) +  " totalGain = " + toString(totalGain) + "\n"); }
                    
                    if (clusterGain <= gains[i]) {
                        clusterGain = gains[i];
                        medoid = i;
                    }
                }
            }
            
            //save medoid value
            medoids.insert(medoid);
            
            if (m->debug) { m->mothurOut("[DEBUG]: new medoid " + toString(medoid) + "\n"); }
            
            //update dp values
            for (int i = 0; i < numSamples; i++) {
                if (Dp[i][0] > dists[i][medoid]) { Dp[i][0] = dists[i][medoid]; }
            }
        }
        if (m->debug) { m->mothurOut("[DEBUG]: done building medoids\n"); }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "Pam", "buildPhase");
		exit(1);
	}
}
/**************************************************************************************************/
//goal to swap medoids with non-medoids to see if we can reduce the overall cost
int Pam::swapPhase() {
    try {
        if (m->debug) { m->mothurOut("[DEBUG]: swapping  medoids\n"); }
        //calculate cost of initial choice - average distance of samples to their closest medoid
        double sky = 0.0;
        double dzsky = 1.0;
        for (int i = 0; i < numSamples; i++) { sky += Dp[i][0]; }  sky /= (double) numSamples;
        
        bool done = false;
        int hbest, nbest; hbest = -1; nbest = -1;
        while (!done) {
            if (m->control_pressed) { break; }
            
            updateDp();
            
            dzsky = 1;
            
            for (int h = 0; h < numSamples; h++) {
                if (m->control_pressed) { break; }
                if (medoids.count(h) == 0) { //this is NOT a medoid
                    for (int i = 0; i < numSamples; i++) {
                        if (medoids.count(i) != 0) { //this is a medoid
                        
                            double dz = 0.0; //Tih sum of distances between objects and closest medoid caused by swapping i and h. Basically the change in cost. If this < 0 its a "good" swap. When all Tih are > 0, then we stop the algo, because we have the optimal medoids.
                            for (int j = 0; j < numSamples; j++) {
                                if (m->control_pressed) { break; }
                                if (dists[i][j] == Dp[j][0]) {
                                    double small = 0.0;
                                    if (Dp[j][1] > dists[h][j]) {   small = dists[h][j];    }
                                    else                        {   small = Dp[j][1];       }
                                    dz += (small - Dp[j][0]);
                                }else if (dists[h][j] < Dp[j][0]) {
                                    dz += (dists[h][j] - Dp[j][0]);
                                }
                            }
                            if (dzsky > dz) {
                                dzsky = dz;
                                hbest = h; 
                                nbest = i;
                            }
                        }//end if medoid
                    }//end for i
                }//end if NOT medoid
            }//end if h
            
            if (dzsky < -16 *DBL_EPSILON * fabs(sky)) {
                medoids.insert(hbest);
                medoids.erase(nbest);
                if (m->debug) { m->mothurOut("[DEBUG]: swapping " + toString(hbest) + " " + toString(nbest) + "\n"); }
                sky += dzsky;
            }else { done = true; } //stop algo.
        }
        
        
        //fill zmatrix
        int count = 0;
        vector<int> tempMedoids;
        for (set<int>::iterator it = medoids.begin(); it != medoids.end(); it++) {
            medoid2Partition[*it] = count;
            zMatrix[count][*it] = 1; count++; //set medoid in this partition.
            tempMedoids.push_back(*it);
        }
        
        //which partition do you belong to?
        laplace = 0;
        for (int i = 0; i < numSamples; i++) {
            int partition = 0;
            double dist = dists[i][tempMedoids[0]]; //assign to first medoid
            for (int j = 1; j < tempMedoids.size(); j++) {
                if (dists[i][tempMedoids[j]] < dist) { //is this medoid closer?
                    dist = dists[i][tempMedoids[j]];
                    partition = j;
                }
            }
            zMatrix[partition][i] = 1;
            laplace += dist;
        }
        laplace /= (double) numSamples;
        
        if (m->debug) {
            for(int i=0;i<numPartitions;i++){
                m->mothurOut("[DEBUG]: partition 1: "); 
                for (int j = 0; j < numSamples; j++) {
                     m->mothurOut(toString(zMatrix[i][j]) + " ");
                }
                m->mothurOut("\n"); 
            }
            m->mothurOut("[DEBUG]: medoids : ");
            for (set<int>::iterator it = medoids.begin(); it != medoids.end(); it++) { m->mothurOut(toString(*it) + " ");
            }
            m->mothurOut("\n");
            
            m->mothurOut("[DEBUG]: laplace : " + toString(laplace));  m->mothurOut("\n");
        }
        
        if (m->debug) { m->mothurOut("[DEBUG]: done swapping  medoids\n"); }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Pam", "swapPhase");
        exit(1);
    }
}

/**************************************************************************************************/
int Pam::updateDp() {
    try {
        for (int j = 0; j < numSamples; j++) {
            if (m->control_pressed) { break; }
            
            //initialize dp and ep
            Dp[j][0] = largestDist; Dp[j][1] = largestDist;
        
            for (int i = 0; i < numSamples; i++) {
                if (medoids.count(i) != 0) { //is this a medoid? 
                    if (Dp[j][0] > dists[j][i]) {
                        Dp[j][0] = Dp[j][1];
                        Dp[j][0] = dists[j][i];
                    }else if (Dp[j][1] > dists[j][i]) {
                        Dp[j][1] = dists[j][i];
                    }
                }
            }
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Pam", "updateDp");
        exit(1);
    }
}
/**************************************************************************************************/
//The silhouette width S(i)of individual data points i is calculated using the following formula:
/*
 s(i) = b(i) - a(i)
 -----------
 max(b(i),a(i))
 where a(i) is the average dissimilarity (or distance) of sample i to all other samples in the same cluster, while b(i) is the average dissimilarity (or distance) to all objects in the closest other cluster.
 
 The formula implies -1 =< S(i) =< 1 . A sample which is much closer to its own cluster than to any other cluster has a high S(i) value, while S(i) close to 0 implies that the given sample lies somewhere between two clusters. Large negative S(i) values indicate that the sample was assigned to the wrong cluster.
 */

vector<double> Pam::calcSilhouettes(vector<vector<double> > dists) {
    try {
        vector<double> silhouettes; silhouettes.resize(numSamples, 0.0);
        if (numPartitions < 2) { return silhouettes; }
        
        vector<int> whoDp;
        for (int i = 0; i < numSamples; i++) { // assumes at least 2 partitions
            vector<seqDist> tempMeds;
            for (set<int>::iterator it = medoids.begin(); it != medoids.end(); it++) {
                if (m->control_pressed) { break; }
                seqDist temp(*it, *it, dists[i][*it]); //medoid, medoid, distance from medoid to sample
                tempMeds.push_back(temp);
            }
            sort(tempMeds.begin(), tempMeds.end(), compareSequenceDistance); //sort by distance
            whoDp.push_back(tempMeds[1].seq1); //second closest medoid
        }
        
        
        vector<double> a, b; a.resize(numSamples, 0.0); b.resize(numSamples, 0.0);
        
        //calc a - all a[i] are the same in the same partition
        for (int k = 0; k < numPartitions; k++) {
            if (m->control_pressed) { break; }
            
            int count = 0;
            double totalA = 0.0;
            for (int i = 0; i < numSamples; i++) {
                for (int j = 0; j < numSamples; j++) {
                    if (m->control_pressed) { break; }
                    if ((zMatrix[k][i] != 0) && (zMatrix[k][j] != 0)){ //are both samples in the partition, if so add there distance
                        totalA += dists[j][i]; //distance from this sample to the other samples in the partition
                        count++;
                    }
                }
            }
            totalA /= (double) count;
            
            //set a[i] to average for cluster
            for (int i = 0; i < numSamples; i++) {
                if (zMatrix[k][i] != 0) { a[i] = totalA; }
            }
        }
        
        //calc b
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { break; }
            
            int nextClosestMedoid = whoDp[i];
            int thisPartition = medoid2Partition[nextClosestMedoid];
            int count = 0;
            double totalB = 0.0;
            for (int j = 0; j < numSamples; j++) {
                if (zMatrix[thisPartition][j] != 0) { //this sample is in this partition
                    totalB += dists[i][j];
                    count++;
                }
            }
            b[i] = totalB / (double) count;
        }
        
        //calc silhouettes
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { break; }
            
            double denom = a[i];
            if (b[i] > denom) { denom = b[i]; } //max(a[i],b[i])
            
            silhouettes[i] = (b[i] - a[i]) / denom;
            
            //cout << "silhouettes " << i << '\t' << silhouettes[i] << endl;
        }
               
        return silhouettes;
    }
    catch(exception& e) {
        m->errorOut(e, "Pam", "calcSilhouettes");
        exit(1);
    }
}
/**************************************************************************************************/
/*To assess the optimal number of clusters our dataset was most robustly partitioned into, we used the Calinski-Harabasz (CH) Index that has shown good performance in recovering the number of clusters. It is defined as:
 
 CHk=Bk/(k−1)/Wk/(n−k)
 
 where Bk is the between-cluster sum of squares (i.e. the squared distances between all points i and j, for which i and j are not in the same cluster) and Wk is the within-clusters sum of squares (i.e. the squared distances between all points i and j, for which i and j are in the same cluster). This measure implements the idea that the clustering is more robust when between-cluster distances are substantially larger than within-cluster distances. Consequently, we chose the number of clusters k such that CHk was maximal.*/
double Pam::calcCHIndex(vector< vector<double> > dists){ //countMatrix = [numSamples][numOtus]
    try {
        double CH = 0.0;
        
        if (numPartitions < 2) { return CH; }
        
        map<int, int> clusterMap; //map sample to partition
        for (int i = 0; i < numPartitions; i++) {
            for (int j = 0; j < numSamples; j++) {
                if (m->control_pressed) { return 0.0; }
                if (zMatrix[i][j] != 0) { clusterMap[j] = i; }
            }
        }
        
        //make countMatrix a relabund
        vector<vector<double> > relativeAbundance(numSamples);
        //get relative abundance
        for(int i=0;i<numSamples;i++){
            if (m->control_pressed) {  return 0; }
            int groupTotal = 0;
            
            relativeAbundance[i].assign(numOTUs, 0.0);
            
            for(int j=0;j<numOTUs;j++){
                groupTotal += countMatrix[i][j];
            }
            for(int j=0;j<numOTUs;j++){
                relativeAbundance[i][j] = countMatrix[i][j] / (double)groupTotal;
            }
        }
        
        //find centers
        vector<vector<double> > centers; centers.resize(numPartitions);
        int countPartitions = 0;
        for (set<int>::iterator it = medoids.begin(); it != medoids.end(); it++) {
            for (int j = 0; j < numOTUs; j++) {
                centers[countPartitions].push_back(relativeAbundance[*it][j]); //save the relative abundance of the medoid for this partition for this OTU
            }
            countPartitions++;
        }
        
        double sumBetweenCluster = 0.0;
        double sumWithinClusters = 0.0;
        
        for (int i = 0; i < numSamples; i++) { //lt
            for (int j = 0; j < i; j++) {
                if (m->control_pressed) { return 0.0; }
                int partitionI = clusterMap[i];
                int partitionJ = clusterMap[j];
                
                if (partitionI == partitionJ) { //they are from the same cluster so this distance is added to Wk
                    sumWithinClusters += (dists[i][j] * dists[i][j]);
                }else { //they are NOT from the same cluster so this distance is added to Bk
                    sumBetweenCluster += (dists[i][j] * dists[i][j]);
                }
            }
        }
        //cout << numPartitions << '\t' << sumWithinClusters << '\t' << sumBetweenCluster << '\t' <<  (numSamples - numPartitions) << endl;
        
        CH = (sumBetweenCluster / (double)(numPartitions - 1)) / (sumWithinClusters / (double) (numSamples - numPartitions));
        
        return CH;
    }
    catch(exception& e){
        m->errorOut(e, "Pam", "calcCHIndex");
        exit(1);
    }
}

/**************************************************************************************************/



