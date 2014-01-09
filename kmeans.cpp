//
//  kmeans.cpp
//  Mothur
//
//  Created by SarahsWork on 12/4/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "kmeans.h"

/**************************************************************************************************/

KMeans::KMeans(vector<vector<int> > cm, int p) : CommunityTypeFinder() {
    try {
        countMatrix = cm;
        numSamples = (int)countMatrix.size();
        numOTUs = (int)countMatrix[0].size();
        numPartitions = p;
        
        findkMeans();
    }
	catch(exception& e) {
		m->errorOut(e, "KMeans", "KMeans");
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

vector<double> KMeans::calcSilhouettes(vector<vector<double> > dists) {
    try {
        vector<double> silhouettes; silhouettes.resize(numSamples, 0.0);
        if (numPartitions < 2) { return silhouettes; }
        
        map<int, int> clusterMap; //map sample to partition
        for (int j = 0; j < numSamples; j++) {
            double maxValue = 0.0;
            for (int i = 0; i < numPartitions; i++) {
                if (m->control_pressed) { return silhouettes; }
                if (zMatrix[i][j] > maxValue) { //for kmeans zmatrix contains values for each sample in each partition. partition with highest value for that sample is the partition where the sample should be
                    clusterMap[j] = i;
                    maxValue = zMatrix[i][j];
                }
            }
        }
        
        vector<int> nextClosestPartition;
        findSecondClosest(nextClosestPartition, dists, clusterMap);
        
        if (m->control_pressed) { return silhouettes; }
        
        vector<double> a, b; a.resize(numSamples, 0.0); b.resize(numSamples, 0.0);
        
        //calc a - all a[i] are the same in the same partition
        for (int k = 0; k < numPartitions; k++) {
            if (m->control_pressed) { break; }
            
            int count = 0;
            double totalA = 0.0;
            for (int i = 0; i < numSamples; i++) {
                for (int j = 0; j < numSamples; j++) {
                    if (m->control_pressed) { break; }
                    if ((clusterMap[i] == k) && (clusterMap[j] == k)){ //are both samples in the partition, if so add there distance
                        totalA += dists[j][i]; //distance from this sample to the other samples in the partition
                        count++;
                    }
                }
            }
            totalA /= (double) count;
            
            //set a[i] to average for cluster
            for (int i = 0; i < numSamples; i++) {
                if (clusterMap[i] == k) { a[i] = totalA; }
            }
        }
        
        //calc b
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { break; }
            
            int thisPartition = nextClosestPartition[i];
            int count = 0;
            double totalB = 0.0;
            for (int j = 0; j < numSamples; j++) {
                if (clusterMap[j] == thisPartition) { //this sample is in this partition
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
        m->errorOut(e, "KMeans", "calcSilhouettes");
        exit(1);
    }
}
/**************************************************************************************************/
int KMeans::findSecondClosest(vector<int>& nextClosestPartition, vector<vector<double> >& dists, map<int, int> clusterMap) {
    try {
        vector<double> minScores; minScores.resize(numSamples, 1e6);
        nextClosestPartition.resize(numSamples, 0);
        
        
        for (int i = 0; i < numSamples; i++) {
            for (int j = 0; j < numPartitions; j++) {
                if (m->control_pressed) { break; }
                
                //is this the one we are assigned to - ie the "best" cluster. We want second best.
                //if numPartitions = 2, then special case??
                if (clusterMap[i] != j) {
                    double score = 1e6;
                    if (numPartitions == 2) {
                        score = 0.0;  //choose other option, there are only 2.
                    }else{   score = calcScore(i, j, dists, clusterMap); }
                
                    if (m->debug) { m->mothurOut("[DEBUG]: sample = " + toString(i) + " partition = " + toString(j) + " score = " + toString(score) + "\n"); }
                    
                    //is this better than our last find
                    if (score < minScores[i]) {
                        minScores[i] = score;
                        nextClosestPartition[i] = j;
                    }
                }else {} //best case, ignore
            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "KMeans", "findSecondClosest");
        exit(1);
    }
}
/**************************************************************************************************/
double KMeans::calcScore(int sample, int partition, vector<vector<double> >& dists, map<int, int> clusterMap) {
    try {
        //square the distances and then for each pair of clusters, calculate the sum of the squraed distances between the clusters
        //then with the sum of hte squared dsitances take the square root and divide by the number of distances in the sum
        
        double sum = 0.0; int count = 0;
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { break; }
            if (clusterMap[i] == partition) { //samples in this cluster
                sum += (dists[sample][i] * dists[sample][i]);
                count++;
            }
        }
        
        sum = sqrt(sum);
        sum /= (double) count;
                
        return sum;
    }
    catch(exception& e) {
        m->errorOut(e, "KMeans", "calcScore");
        exit(1);
    }
}
/**************************************************************************************************/
/*To assess the optimal number of clusters our dataset was most robustly partitioned into, we used the Calinski-Harabasz (CH) Index that has shown good performance in recovering the number of clusters. It is defined as:
 
 CHk=Bk/(k−1)/Wk/(n−k)
 
 where Bk is the between-cluster sum of squares (i.e. the squared distances between all points i and j, for which i and j are not in the same cluster) and Wk is the within-clusters sum of squares (i.e. the squared distances between all points i and j, for which i and j are in the same cluster). This measure implements the idea that the clustering is more robust when between-cluster distances are substantially larger than within-cluster distances. Consequently, we chose the number of clusters k such that CHk was maximal.*/
double KMeans::calcCHIndex(vector< vector< double> > dists){
    try {
        double CH = 0.0;
        
        if (numPartitions < 2) { return CH; }
        
        map<int, int> clusterMap; //map sample to partition
        for (int j = 0; j < numSamples; j++) {
            double maxValue = 0.0;
            for (int i = 0; i < numPartitions; i++) {
                if (m->control_pressed) { return 0.0; }
                if (zMatrix[i][j] > maxValue) { //for kmeans zmatrix contains values for each sample in each partition. partition with highest value for that sample is the partition where the sample should be
                    clusterMap[j] = i;
                    maxValue = zMatrix[i][j];
                }
            }
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
        m->errorOut(e, "KMeans", "calcCHIndex");
        exit(1);
    }
}
/**************************************************************************************************/




