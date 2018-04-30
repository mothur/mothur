//
//  communitytype.cpp
//  Mothur
//
//  Created by SarahsWork on 12/3/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "communitytype.h"

/**************************************************************************************************/

//can we get these psi/psi1 calculations into their own math class?
//psi calcualtions swiped from gsl library...

static const double psi_cs[23] = {
    -.038057080835217922,
    .491415393029387130,
    -.056815747821244730,
    .008357821225914313,
    -.001333232857994342,
    .000220313287069308,
    -.000037040238178456,
    .000006283793654854,
    -.000001071263908506,
    .000000183128394654,
    -.000000031353509361,
    .000000005372808776,
    -.000000000921168141,
    .000000000157981265,
    -.000000000027098646,
    .000000000004648722,
    -.000000000000797527,
    .000000000000136827,
    -.000000000000023475,
    .000000000000004027,
    -.000000000000000691,
    .000000000000000118,
    -.000000000000000020
};

static double apsi_cs[16] = {
    -.0204749044678185,
    -.0101801271534859,
    .0000559718725387,
    -.0000012917176570,
    .0000000572858606,
    -.0000000038213539,
    .0000000003397434,
    -.0000000000374838,
    .0000000000048990,
    -.0000000000007344,
    .0000000000001233,
    -.0000000000000228,
    .0000000000000045,
    -.0000000000000009,
    .0000000000000002,
    -.0000000000000000
};

/**************************************************************************************************/
/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
    1.00000000000000000000000000000,
    0.083333333333333333333333333333,
    -0.00138888888888888888888888888889,
    0.000033068783068783068783068783069,
    -8.2671957671957671957671957672e-07,
    2.0876756987868098979210090321e-08,
    -5.2841901386874931848476822022e-10,
    1.3382536530684678832826980975e-11,
    -3.3896802963225828668301953912e-13,
    8.5860620562778445641359054504e-15,
    -2.1748686985580618730415164239e-16,
    5.5090028283602295152026526089e-18,
    -1.3954464685812523340707686264e-19,
    3.5347070396294674716932299778e-21,
    -8.9535174270375468504026113181e-23
};

/**************************************************************************************************/
void CommunityTypeFinder::printSilData(ofstream& out, double chi, vector<double> sils){
    try {
        out << setprecision (6) << numPartitions << '\t'  << chi;
        for (int i = 0; i < sils.size(); i++) {
            out << '\t' << sils[i];
        }
        out << endl;
        
        return;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "printSilData");
        exit(1);
    }
}
/**************************************************************************************************/
void CommunityTypeFinder::printSilData(ostream& out, double chi, vector<double> sils){
    try {
        out << setprecision (6) << numPartitions << '\t'  << chi;
        m->mothurOutJustToLog(toString(numPartitions) + '\t' + toString(chi));
        for (int i = 0; i < sils.size(); i++) {
            out << '\t' << sils[i];
            m->mothurOutJustToLog("\t" + toString(sils[i]));
        }
        out << endl;
        m->mothurOutJustToLog("\n");
        
        return;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "printSilData");
        exit(1);
    }
}
/**************************************************************************************************/

void CommunityTypeFinder::printZMatrix(string fileName, vector<string> sampleName){
    try {
        ofstream printMatrix;
        util.openOutputFile(fileName, printMatrix); //(fileName.c_str());
        printMatrix.setf(ios::fixed, ios::floatfield);
        printMatrix.setf(ios::showpoint);
        
        for(int i=0;i<numPartitions;i++){   printMatrix << "\tPartition_" << i+1;   }   printMatrix << endl;
        
        for(int i=0;i<numSamples;i++){
            printMatrix << sampleName[i];
            for(int j=0;j<numPartitions;j++){
                printMatrix << setprecision(4) << '\t' << zMatrix[j][i];
            }
            printMatrix << endl;
        }
        printMatrix.close();
    }
	catch(exception& e) {
		m->errorOut(e, "CommunityTypeFinder", "printZMatrix");
		exit(1);
	}
}

/**************************************************************************************************/

void CommunityTypeFinder::printRelAbund(string fileName, vector<string> otuNames){
    try {
        ofstream printRA;
        util.openOutputFile(fileName, printRA); //(fileName.c_str());
        printRA.setf(ios::fixed, ios::floatfield);
        printRA.setf(ios::showpoint);
        
        vector<double> totals(numPartitions, 0.0000);
        for(int i=0;i<numPartitions;i++){
            for(int j=0;j<numOTUs;j++){
                totals[i] += exp(lambdaMatrix[i][j]);
            }
        }
        
        printRA << "Taxon";
        for(int i=0;i<numPartitions;i++){
            printRA << "\tPartition_" << i+1 << '_' << setprecision(4) << totals[i];
            printRA << "\tPartition_" << i+1 <<"_LCI" << "\tPartition_" << i+1 << "_UCI";
        }
        printRA << endl;
        
        for(int i=0;i<numOTUs;i++){
            
            if (m->getControl_pressed()) { break; }
            
            printRA << otuNames[i];
            for(int j=0;j<numPartitions;j++){
                
                if(error[j][i] >= 0.0000){
                    double std = sqrt(error[j][i]);
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i]) / totals[j];
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i] - 2.0 * std) / totals[j];
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i] + 2.0 * std) / totals[j];
                }
                else{
                    printRA << '\t' << 100 * exp(lambdaMatrix[j][i]) / totals[j];
                    printRA << '\t' << "NA";
                    printRA << '\t' << "NA";
                }
            }
            printRA << endl;
        }
        
        printRA.close();
    }
	catch(exception& e) {
		m->errorOut(e, "CommunityTypeFinder", "printRelAbund");
		exit(1);
	}
}



/**************************************************************************************************/

vector<vector<double> > CommunityTypeFinder::getHessian(){
    try {
        vector<double> alpha(numOTUs, 0.0000);
        double alphaSum = 0.0000;
        
        vector<double> pi = zMatrix[currentPartition];
        vector<double> psi_ajk(numOTUs, 0.0000);
        vector<double> psi_cjk(numOTUs, 0.0000);
        vector<double> psi1_ajk(numOTUs, 0.0000);
        vector<double> psi1_cjk(numOTUs, 0.0000);
        
        for(int j=0;j<numOTUs;j++){
            
            if (m->getControl_pressed()) {  break; }
            
            alpha[j] = exp(lambdaMatrix[currentPartition][j]);
            alphaSum += alpha[j];
            
            for(int i=0;i<numSamples;i++){
                double X = (double) countMatrix[i][j];
                
                psi_ajk[j] += pi[i] * psi(alpha[j]);
                psi1_ajk[j] += pi[i] * psi1(alpha[j]);
                
                psi_cjk[j] += pi[i] * psi(alpha[j] + X);
                psi1_cjk[j] += pi[i] * psi1(alpha[j] + X);
            }
        }
        
        
        double psi_Ck = 0.0000;
        double psi1_Ck = 0.0000;
        
        double weight = 0.0000;
        
        for(int i=0;i<numSamples;i++){
            if (m->getControl_pressed()) {  break; }
            weight += pi[i];
            double sum = 0.0000;
            for(int j=0;j<numOTUs;j++){     sum += alpha[j] + countMatrix[i][j];    }
            
            psi_Ck += pi[i] * psi(sum);
            psi1_Ck += pi[i] * psi1(sum);
        }
        
        double psi_Ak = weight * psi(alphaSum);
        double psi1_Ak = weight * psi1(alphaSum);
        
        vector<vector<double> > hessian(numOTUs);
        for(int i=0;i<numOTUs;i++){ hessian[i].assign(numOTUs, 0.0000); }
        
        for(int i=0;i<numOTUs;i++){
            if (m->getControl_pressed()) {  break; }
            double term1 = -alpha[i] * (- psi_ajk[i] + psi_Ak + psi_cjk[i] - psi_Ck);
            double term2 = -alpha[i] * alpha[i] * (-psi1_ajk[i] + psi1_Ak + psi1_cjk[i] - psi1_Ck);
            double term3 = 0.1 * alpha[i];
            
            hessian[i][i] = term1 + term2 + term3;
            
            for(int j=0;j<i;j++){
                hessian[i][j] = - alpha[i] * alpha[j] * (psi1_Ak - psi1_Ck);
                hessian[j][i] = hessian[i][j];
            }
        }
        
        return hessian;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "getHessian");
        exit(1);
    }
}
/**************************************************************************************************/

double CommunityTypeFinder::psi1(double xx){
    try {
        
        /* Euler-Maclaurin summation formula
         * [Moshier, p. 400, with several typo corrections]
         */
        
        double s = 2.0000;
        const int jmax = 12;
        const int kmax = 10;
        int j, k;
        const double pmax  = pow(kmax + xx, -s);
        double scp = s;
        double pcp = pmax / (kmax + xx);
        double value = pmax*((kmax+xx)/(s-1.0) + 0.5);
        
        for(k=0; k<kmax; k++) {
            if (m->getControl_pressed()) {  return 0; }
            value += pow(k + xx, -s);
        }
        
        for(j=0; j<=jmax; j++) {
            if (m->getControl_pressed()) {  return 0; }
            double delta = hzeta_c[j+1] * scp * pcp;
            value += delta;
            
            if(fabs(delta/value) < 0.5*EPSILON) break;
            
            scp *= (s+2*j+1)*(s+2*j+2);
            pcp /= (kmax + xx)*(kmax + xx);
        }
        
        return value;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "psi1");
        exit(1);
    }
}

/**************************************************************************************************/

double CommunityTypeFinder::psi(double xx){
    try {
        double psiX = 0.0000;
        
        if(xx < 1.0000){
            
            double t1 = 1.0 / xx;
            psiX = cheb_eval(psi_cs, 22, 2.0*xx-1.0);
            psiX = -t1 + psiX;
            
        }
        else if(xx < 2.0000){
            
            const double v = xx - 1.0;
            psiX = cheb_eval(psi_cs, 22, 2.0*v-1.0);
            
        }
        else{
            const double t = 8.0/(xx*xx)-1.0;
            psiX = cheb_eval(apsi_cs, 15, t);
            psiX += log(xx) - 0.5/xx;
        }
        
        return psiX;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "psi");
        exit(1);
    }
}
/**************************************************************************************************/

double CommunityTypeFinder::cheb_eval(const double seriesData[], int order, double xx){
    try {
        double d = 0.0000;
        double dd = 0.0000;
        
        double x2 = xx * 2.0000;
        
        for(int j=order;j>=1;j--){
            if (m->getControl_pressed()) {  return 0; }
            double temp = d;
            d = x2 * d - dd + seriesData[j];
            dd = temp;
        }
        
        d = xx * d - dd + 0.5 * seriesData[0];
        return d;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "cheb_eval");
        exit(1);
    }
}
/**************************************************************************************************/

int CommunityTypeFinder::findkMeans(){
    try {
        error.resize(numPartitions); for (int i = 0; i < numPartitions; i++) { error[i].resize(numOTUs, 0.0); }
        vector<vector<double> > relativeAbundance(numSamples);
        vector<vector<double> > alphaMatrix;
        
        alphaMatrix.resize(numPartitions);
        lambdaMatrix.resize(numPartitions);
        for(int i=0;i<numPartitions;i++){
            alphaMatrix[i].assign(numOTUs, 0);
            lambdaMatrix[i].assign(numOTUs, 0);
        }
        
        //get relative abundance
        for(int i=0;i<numSamples;i++){
            if (m->getControl_pressed()) {  return 0; }
            int groupTotal = 0;
            
            relativeAbundance[i].assign(numOTUs, 0.0);
            
            for(int j=0;j<numOTUs;j++){
                groupTotal += countMatrix[i][j];
            }
            for(int j=0;j<numOTUs;j++){
                relativeAbundance[i][j] = countMatrix[i][j] / (double)groupTotal;
            }
        }
        
        //randomly assign samples into partitions
        zMatrix.resize(numPartitions);
        for(int i=0;i<numPartitions;i++){
            zMatrix[i].assign(numSamples, 0);
        }
        
        //randomize samples
        vector<int> temp;
        for (int i = 0; i < numSamples; i++) { temp.push_back(i); }
        util.mothurRandomShuffle(temp);
        
        //assign each partition at least one random sample
        int numAssignedSamples = 0;
        for (int i = 0; i < numPartitions; i++) {
            zMatrix[i][temp[numAssignedSamples]] = 1;
            numAssignedSamples++;
        }
        
        //assign rest of samples to partitions
        int count = 0;
        for(int i=numAssignedSamples;i<numSamples;i++){
            zMatrix[count%numPartitions][temp[i]] = 1;
            count++;
        }
        
        double maxChange = 1;
        int maxIters = 1000;
        int iteration = 0;
        
        weights.assign(numPartitions, 0);
        
        while(maxChange > 1e-6 && iteration < maxIters){
            
            if (m->getControl_pressed()) {  return 0; }
            //calcualte average relative abundance
            maxChange = 0.0000;
            for(int i=0;i<numPartitions;i++){
                
                double normChange = 0.0;
                
                weights[i] = 0;
                
                for(int j=0;j<numSamples;j++){
                    weights[i] += (double)zMatrix[i][j];
                }
                
                vector<double> averageRelativeAbundance(numOTUs, 0);
                for(int j=0;j<numOTUs;j++){
                    for(int k=0;k<numSamples;k++){
                        averageRelativeAbundance[j] += zMatrix[i][k] * relativeAbundance[k][j];
                    }
                }
                
                for(int j=0;j<numOTUs;j++){
                    averageRelativeAbundance[j] /= weights[i];
                    
                    double difference = averageRelativeAbundance[j] - alphaMatrix[i][j];
                    normChange += difference * difference;
                    alphaMatrix[i][j] = averageRelativeAbundance[j];
                }
                
                normChange = sqrt(normChange);
                
                if(normChange > maxChange){ maxChange = normChange; }
            }
            
            
            //calcualte distance between each sample in partition and the average relative abundance
            for(int i=0;i<numSamples;i++){
                if (m->getControl_pressed()) {  return 0; }
                
                double normalizationFactor = 0;
                vector<double> totalDistToPartition(numPartitions, 0);
                
                for(int j=0;j<numPartitions;j++){
                    for(int k=0;k<numOTUs;k++){
                        double difference = alphaMatrix[j][k] - relativeAbundance[i][k];
                        totalDistToPartition[j] += difference * difference;
                    }
                    totalDistToPartition[j] = sqrt(totalDistToPartition[j]);
                    normalizationFactor += exp(-50.0 * totalDistToPartition[j]);
                }
                
                
                for(int j=0;j<numPartitions;j++){
                    zMatrix[j][i] = exp(-50.0 * totalDistToPartition[j]) / normalizationFactor;
                }
                
            }
            
            iteration++;
            //        cout << "K means: " << iteration << '\t' << maxChange << endl;
            
        }
        
        //    cout << "Iter:-1";
        for(int i=0;i<numPartitions;i++){
            weights[i] = 0.0000;
            
            for(int j=0;j<numSamples;j++){
                weights[i] += zMatrix[i][j];
            }
            //        printf("\tw_%d=%.3f", i, weights[i]);
        }
        //    cout << endl;
        
        
        for(int i=0;i<numOTUs;i++){
            if (m->getControl_pressed()) {  return 0; }
            for(int j=0;j<numPartitions;j++){
                if(alphaMatrix[j][i] > 0){
                    lambdaMatrix[j][i] = log(alphaMatrix[j][i]);
                }
                else{
                    lambdaMatrix[j][i] = -10.0;
                }
            }
        }
        return 0;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "kMeans");
        exit(1);
    }
}

/**************************************************************************************************/
//based on r function .medoid
//results is length numOTUs and holds the distances from x of the sample in d with the min sum of distances to all other samples.
//Basically the "best" medoid.
//returns the sum of the distances squared
double CommunityTypeFinder::rMedoid(vector< vector<double> > x, vector< vector<double> > d){
    try {
        vector<double> results; results.resize(numOTUs, 0.0);
        double minSumDist = 1e6;
        int minGroup = -1;
        
        for (int i = 0; i < d.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            double thisSum = 0.0;
            for (int j = 0; j < d[i].size(); j++) { thisSum += d[i][j];  }
            if (thisSum < minSumDist) {
                minSumDist = thisSum;
                minGroup = i;
            }
        }
        
        if (minGroup != -1) {
            for (int i = 0; i < numOTUs; i++) {  results[i] = x[minGroup][i];  } //save minGroups relativeAbundance for each OTU
        }else { m->mothurOut("[ERROR]: unable to find rMedoid group.\n"); m->setControl_pressed(true); }
        
        
        double allMeanDist = 0.0;
        for (int i = 0; i < x.size(); i++) { //numSamples
            for (int j = 0; j < x[i].size(); j++) { //numOTus
                if (m->getControl_pressed()) { break; }
                allMeanDist += ((x[i][j]-results[j])*(x[i][j]-results[j])); //(otuX sampleY - otuX bestMedoid)^2
                
            }
        }
        return allMeanDist;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "rMedoid");
        exit(1);
    }
}
/**************************************************************************************************/

/*To assess the optimal number of clusters our dataset was most robustly partitioned into, we used the Calinski-Harabasz (CH) Index that has shown good performance in recovering the number of clusters. It is defined as:
 
 CHk=Bk/(k−1)/Wk/(n−k)
 
 where Bk is the between-cluster sum of squares (i.e. the squared distances between all points i and j, for which i and j are not in the same cluster) and Wk is the within-clusters sum of squares (i.e. the squared distances between all points i and j, for which i and j are in the same cluster). This measure implements the idea that the clustering is more robust when between-cluster distances are substantially larger than within-cluster distances. Consequently, we chose the number of clusters k such that CHk was maximal.*/
double CommunityTypeFinder::calcCHIndex(vector< vector< double> > dists){
    try {
        double CH = 0.0;
        
        if (numPartitions < 2) { return CH; }
        
        map<int, int> clusterMap; //map sample to partition
        for (int j = 0; j < numSamples; j++) {
            double maxValue = -1e6;
            for (int i = 0; i < numPartitions; i++) {
                if (m->getControl_pressed()) { return 0.0; }
                if (zMatrix[i][j] > maxValue) { //for kmeans zmatrix contains values for each sample in each partition. partition with highest value for that sample is the partition where the sample should be
                    clusterMap[j] = i;
                    maxValue = zMatrix[i][j];
                }
            }
        }
        
        //make countMatrix a relabund
        vector<vector<double> > relativeAbundance(numSamples); //[numSamples][numOTUs]
        //get relative abundance
        for(int i=0;i<numSamples;i++){
            if (m->getControl_pressed()) {  return 0; }
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
        vector<vector<double> > centers = calcCenters(dists, clusterMap, relativeAbundance);
        
        if (m->getControl_pressed()) { return 0.0; }
        
        double allMeanDist = rMedoid(relativeAbundance, dists);
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: allMeandDist = " + toString(allMeanDist) + "\n"); }
        
        for (int i = 0; i < relativeAbundance.size(); i++) {//numSamples
            for (int j = 0; j < relativeAbundance[i].size(); j++) { //numOtus
                if (m->getControl_pressed()) {  return 0; }
                //x <- (x - centers[cl, ])^2
                relativeAbundance[i][j] = ((relativeAbundance[i][j] - centers[clusterMap[i]][j])*(relativeAbundance[i][j] - centers[clusterMap[i]][j]));
            }
        }
        
        double wgss = 0.0;
        for (int j = 0; j < numOTUs; j++) {
            for(int i=0;i<numSamples;i++){
                if (m->getControl_pressed()) { return 0.0; }
                wgss += relativeAbundance[i][j];
            }
        }
        
        double bgss = allMeanDist - wgss;
        
        CH = (bgss / (double)(numPartitions - 1)) / (wgss / (double) (numSamples - numPartitions));
        
        return CH;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "calcCHIndex");
        exit(1);
    }
}


/**************************************************************************************************/
vector<vector<double> > CommunityTypeFinder::calcCenters(vector<vector<double> >& dists, map<int, int> clusterMap, vector<vector<double> >& relativeAbundance) { //[numsamples][numsamples]
    try {
        //for each partition
        //choose sample with smallest sum of squared dists
        //       cout << "Here" << clusterMap.size() << endl;
        //       for(map<int, int>::iterator it = clusterMap.begin(); it != clusterMap.end(); it++) { cout << it->first << '\t' << it->second <<endl; }
        vector<vector<double> > centers; centers.resize(numPartitions);
        vector<double>  sums;  sums.resize(numSamples, 0.0);
        map<int, vector<int> > partition2Samples; //maps partitions to samples in the partition
        map<int, vector<int> >::iterator it;
        
        for (int i = 0; i < numSamples; i++) {
            int partitionI = clusterMap[i];
            
            //add this sample to list of samples in this partition for access later
            it = partition2Samples.find(partitionI);
            if (it == partition2Samples.end()) {
                vector<int> temp; temp.push_back(i);
                partition2Samples[partitionI] = temp;
            }else {  partition2Samples[partitionI].push_back(i); }
            
            for (int j = 0; j < numSamples; j++) {
                
                int partitionJ = clusterMap[j];
                
                if (partitionI == partitionJ) { //if you are a distance between samples in the same cluster
                    sums[i] += dists[i][j];
                    sums[j] += dists[i][j];
                }else{}//we dont' care about distance between clusters
            }
        }
        
        vector<int> medoidsVector; medoidsVector.resize(numPartitions, -1);
        for (it = partition2Samples.begin(); it != partition2Samples.end(); it++) { //for each partition look for sample with smallest squared
            
            //sum dist to all other samples in cluster
            vector<int> members = it->second;
            double minSumDist = 1e6;
            for (int i = 0; i < members.size(); i++) {
                if (m->getControl_pressed()) { return centers; }
                if (sums[members[i]] < minSumDist) {
                    minSumDist = sums[members[i]];
                    medoidsVector[it->first] = members[i];
                }
            }
            
        }
        
        set<int> medoids;
        for (int i = 0; i < medoidsVector.size(); i++) {
            medoids.insert(medoidsVector[i]);
        }
        
        int countPartitions = 0;
        for (set<int>::iterator it = medoids.begin(); it != medoids.end(); it++) {
            for (int j = 0; j < numOTUs; j++) {
                centers[countPartitions].push_back(relativeAbundance[*it][j]); //save the relative abundance of the medoid for this partition for this OTU
            }
            countPartitions++;
        }
        
        return centers;
    }
    catch(exception& e){
        m->errorOut(e, "CommunityTypeFinder", "calcCenters");
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
//based on silouette.r which calls sildist.c written by Francois Romain
vector<double> CommunityTypeFinder::calcSilhouettes(vector<vector<double> > dists) {
    try {
        vector<double> silhouettes; silhouettes.resize(numSamples, 0.0);
        if (numPartitions < 2) { return silhouettes; }
        
        
        map<int, int> clusterMap; //map sample to partition
        for (int j = 0; j < numSamples; j++) {
            double maxValue = 0.0;
            for (int i = 0; i < numPartitions; i++) {
                if (m->getControl_pressed()) { return silhouettes; }
                if (zMatrix[i][j] > maxValue) { //for kmeans zmatrix contains values for each sample in each partition. partition with highest value for that sample is the partition where the sample should be
                    clusterMap[j] = i;
                    maxValue = zMatrix[i][j];
                }
            }
        }
        
        //count number of samples in each partition
        vector<int> counts; counts.resize(numPartitions, 0);
        vector<double> DiC; DiC.resize((numPartitions*numSamples), 0.0);
        bool computeSi = true;
        
        for (int i = 0; i < numSamples; i++) {
            int partitionI = clusterMap[i];
            counts[partitionI]++;
            
            for (int j = i+1; j < numSamples; j++) {
                if (m->getControl_pressed()) { return silhouettes; }
                int partitionJ = clusterMap[j];
                
                DiC[numPartitions*i+partitionJ] += dists[i][j];
                DiC[numPartitions*j+partitionI] += dists[i][j];
            }
        }
        
        vector<int> neighbor; neighbor.resize(numSamples, -1);
        for (int i = 0; i < numSamples; i++) {
            if (m->getControl_pressed()) { return silhouettes; }
            int ki = numPartitions*i;
            int partitionI = clusterMap[i];
            computeSi = true;
            
            for (int j = 0; j < numPartitions; j++) {
                if (j == partitionI) {
                    if (counts[j] == 1) { //only one sample in cluster
                        computeSi = false;
                    }else { DiC[ki+j] /= (counts[j]-1); }
                }else{
                    DiC[ki+j] /= counts[j];
                }
            }
            
            double ai = DiC[ki+partitionI];
            
            double bi = 0.0;
            if (partitionI == 0) {  bi = DiC[ki+1]; neighbor[i] = 2; }
            else {  bi =  DiC[ki]; neighbor[i] = 1; }
            
            for (int j = 1; j < numPartitions; j++) {
                if (j != partitionI) {
                    if (bi > DiC[ki+j]) {
                        bi = DiC[ki + j];
                        neighbor[i] = j+1;
                    }
                }
            }
            
            silhouettes[i] = 0.0;
            if (computeSi && bi != ai) {
                silhouettes[i] = (bi-ai) / (max(ai, bi));
            }
        }
        
        return silhouettes;
    }
    catch(exception& e) {
        m->errorOut(e, "CommunityTypeFinder", "calcSilhouettes");
        exit(1);
    }
}
/**************************************************************************************************/

