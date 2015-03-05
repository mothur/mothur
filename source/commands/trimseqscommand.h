#ifndef TRIMSEQSCOMMAND_H
#define TRIMSEQSCOMMAND_H

/*
 *  trimseqscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"
#include "qualityscores.h"
#include "trimoligos.h"
#include "counttable.h"


class TrimSeqsCommand : public Command {
public:
	TrimSeqsCommand(string);
	TrimSeqsCommand();
	~TrimSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "trim.seqs";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Trim.seqs"; }
	string getDescription()		{ return "provides the preprocessing features needed to screen and sort pyrosequences"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:    
	bool getOligos(vector<vector<string> >&, vector<vector<string> >&, vector<vector<string> >&);
	bool keepFirstTrim(Sequence&, QualityScores&);
	bool removeLastTrim(Sequence&, QualityScores&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
    string reverseOligo(string);

	bool abort, createGroup;
	string fastaFile, oligoFile, qFileName, groupfile, nameFile, countfile, outputDir;
	
	bool flip, allFiles, qtrim, keepforward, pairedOligos, reorient, logtransform;
	int numFPrimers, numRPrimers, numLinkers, numSpacers, maxAmbig, maxHomoP, minLength, maxLength, processors, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs, comboStarts;
	int qWindowSize, qWindowStep, keepFirst, removeLast;
	double qRollAverage, qThreshold, qWindowAverage, qAverage;
	vector<string> revPrimer, outputNames;
	set<string> filesToRemove;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
	map<string, int> barcodes;
	vector<string> groupVector;
	map<string, int> primers;
    vector<string>  linker;
    vector<string>  spacer;
	map<string, int> combos;
	map<string, int> groupToIndex;
	vector<string> primerNameVector;	//needed here?
	vector<string> barcodeNameVector;	//needed here?
	map<string, int> groupCounts;  
	map<string, string> nameMap;
    map<string, int> nameCount; //for countfile name -> repCount
    map<string, string> groupMap; //for countfile name -> group

	vector<int> processIDS;   //processid
	vector<linePair> lines;
	vector<linePair> qLines;
	
	int driverCreateTrim(string, string, string, string, string, string, string, string, string, string, string, vector<vector<string> >, vector<vector<string> >, vector<vector<string> >, linePair, linePair);	
	int createProcessesCreateTrim(string, string, string, string, string, string, string, string, string, string, string, vector<vector<string> >, vector<vector<string> >, vector<vector<string> >);
	int setLines(string, string);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct trimData {
    unsigned long long start, end;
    MothurOut* m;
    string filename, qFileName, trimFileName, scrapFileName, trimQFileName, scrapQFileName, trimNFileName, scrapNFileName, trimCFileName, scrapCFileName, groupFileName, nameFile, countfile;
	vector<vector<string> > fastaFileNames;
    vector<vector<string> > qualFileNames;
    vector<vector<string> > nameFileNames;
    unsigned long long lineStart, lineEnd, qlineStart, qlineEnd;
    bool flip, allFiles, qtrim, keepforward, createGroup, pairedOligos, reorient, logtransform;
	int numFPrimers, numRPrimers, numLinkers, numSpacers, maxAmbig, maxHomoP, minLength, maxLength, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs;
	int qWindowSize, qWindowStep, keepFirst, removeLast, count;
	double qRollAverage, qThreshold, qWindowAverage, qAverage;
    vector<string> revPrimer;
	map<string, int> barcodes;
	map<string, int> primers;
    map<string, int> nameCount;
    vector<string>  linker;
    vector<string>  spacer;
	map<string, int> combos;
	vector<string> primerNameVector;	
	vector<string> barcodeNameVector;	
	map<string, int> groupCounts;  
	map<string, string> nameMap;
    map<string, string> groupMap;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
    
	trimData(){}
	trimData(string fn, string qn, string nf, string cf, string tn, string sn, string tqn, string sqn, string tnn, string snn, string tcn, string scn,string gn, vector<vector<string> > ffn, vector<vector<string> > qfn, vector<vector<string> > nfn, unsigned long long lstart, unsigned long long lend, unsigned long long qstart, unsigned long long qend,  MothurOut* mout,
                      int pd, int bd, int ld, int sd, int td, map<string, int> pri, map<string, int> bar, vector<string> revP, vector<string> li, vector<string> spa, map<int, oligosPair> pbr, map<int, oligosPair> ppr, bool po,
                      vector<string> priNameVector, vector<string> barNameVector, bool cGroup, bool aFiles, bool keepF, int keepfi, int removeL,
                      int WindowStep, int WindowSize, int WindowAverage, bool trim, double Threshold, double Average, double RollAverage, bool lt,
                      int minL, int maxA, int maxH, int maxL, bool fli, bool reo, map<string, string> nm, map<string, int> ncount) {
        filename = fn;
        qFileName = qn;
        nameFile = nf;
        countfile = cf;
        trimFileName = tn;
        scrapFileName = sn;
        trimQFileName = tqn;
        scrapQFileName = sqn;
        trimNFileName = tnn;
        scrapNFileName = snn;
        trimCFileName = tcn;
        scrapCFileName = scn;
        groupFileName = gn;
        fastaFileNames = ffn;
        qualFileNames = qfn;
        nameFileNames = nfn;
        lineStart = lstart;
        lineEnd = lend;
        qlineStart = qstart;
        qlineEnd = qend;
		m = mout;
        nameCount = ncount;
        
        pdiffs = pd;
        bdiffs = bd;
        ldiffs = ld;
        sdiffs = sd;
        tdiffs = td;
        barcodes = bar;
        pairedPrimers = ppr;
        pairedBarcodes = pbr;
        pairedOligos = po;
        primers = pri;      numFPrimers = primers.size();
        revPrimer = revP;   numRPrimers = revPrimer.size();
        linker = li;        numLinkers = linker.size();
        spacer = spa;       numSpacers = spacer.size();
        primerNameVector = priNameVector;
        barcodeNameVector = barNameVector;
        
        createGroup = cGroup;
        allFiles = aFiles;
        keepforward = keepF;
        keepFirst = keepfi;
        removeLast = removeL;
        qWindowStep = WindowStep;
        qWindowSize = WindowSize;
        qWindowAverage = WindowAverage;
        qtrim = trim;
        qThreshold = Threshold;
        qAverage = Average;
        qRollAverage = RollAverage;
        logtransform = lt;
        minLength = minL;
        maxAmbig = maxA;
        maxHomoP = maxH;
        maxLength = maxL;
        flip = fli;
        reorient = reo;
        nameMap = nm;
        count = 0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyTrimThreadFunction(LPVOID lpParam){ 
	trimData* pDataArray;
	pDataArray = (trimData*)lpParam;
	
	try {
        ofstream trimFASTAFile;
		pDataArray->m->openOutputFile(pDataArray->trimFileName, trimFASTAFile);
		
		ofstream scrapFASTAFile;
		pDataArray->m->openOutputFile(pDataArray->scrapFileName, scrapFASTAFile);
		
		ofstream trimQualFile;
		ofstream scrapQualFile;
		if(pDataArray->qFileName != ""){
			pDataArray->m->openOutputFile(pDataArray->trimQFileName, trimQualFile);
			pDataArray->m->openOutputFile(pDataArray->scrapQFileName, scrapQualFile);
		}
		
		ofstream trimNameFile;
		ofstream scrapNameFile;
		if(pDataArray->nameFile != ""){
			pDataArray->m->openOutputFile(pDataArray->trimNFileName, trimNameFile);
			pDataArray->m->openOutputFile(pDataArray->scrapNFileName, scrapNameFile);
		}
		
		
		ofstream outGroupsFile;
		if ((pDataArray->createGroup) && (pDataArray->countfile == "")){	pDataArray->m->openOutputFile(pDataArray->groupFileName, outGroupsFile);   }
		if(pDataArray->allFiles){
			for (int i = 0; i < pDataArray->fastaFileNames.size(); i++) { //clears old file
				for (int j = 0; j < pDataArray->fastaFileNames[i].size(); j++) { //clears old file
					if (pDataArray->fastaFileNames[i][j] != "") {
						ofstream temp;
						pDataArray->m->openOutputFile(pDataArray->fastaFileNames[i][j], temp);			temp.close();
						if(pDataArray->qFileName != ""){
							pDataArray->m->openOutputFile(pDataArray->qualFileNames[i][j], temp);			temp.close();
						}
						
						if(pDataArray->nameFile != ""){
							pDataArray->m->openOutputFile(pDataArray->nameFileNames[i][j], temp);			temp.close();
						}
					}
				}
			}
		}
		
        ofstream trimCountFile;
		ofstream scrapCountFile;
		if(pDataArray->countfile != ""){
			pDataArray->m->openOutputFile(pDataArray->trimCFileName, trimCountFile);
			pDataArray->m->openOutputFile(pDataArray->scrapCFileName, scrapCountFile);
            if ((pDataArray->lineStart == 0) || (pDataArray->lineStart == 1)) { trimCountFile << "Representative_Sequence\ttotal" << endl; scrapCountFile << "Representative_Sequence\ttotal" << endl; }
		}
        
		ifstream inFASTA;
		pDataArray->m->openInputFile(pDataArray->filename, inFASTA);
		if ((pDataArray->lineStart == 0) || (pDataArray->lineStart == 1)) {
			inFASTA.seekg(0);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->lineStart-1); pDataArray->m->gobble(inFASTA); 
		}
		
		ifstream qFile;
		if(pDataArray->qFileName != "")	{
			pDataArray->m->openInputFile(pDataArray->qFileName, qFile);
			if ((pDataArray->qlineStart == 0) || (pDataArray->qlineStart == 1)) {
                qFile.seekg(0);
            }else { //this accounts for the difference in line endings. 
                qFile.seekg(pDataArray->qlineStart-1); pDataArray->m->gobble(qFile); 
            } 
		}
		
        TrimOligos* trimOligos = NULL;
        int numBarcodes = pDataArray->barcodes.size();
        if (pDataArray->pairedOligos)   {   trimOligos = new TrimOligos(pDataArray->pdiffs, pDataArray->bdiffs, 0, 0, pDataArray->pairedPrimers, pDataArray->pairedBarcodes);   numBarcodes = pDataArray->pairedBarcodes.size(); pDataArray->numFPrimers = pDataArray->pairedPrimers.size(); }
        else                {   trimOligos = new TrimOligos(pDataArray->pdiffs, pDataArray->bdiffs, pDataArray->ldiffs, pDataArray->sdiffs, pDataArray->primers, pDataArray->barcodes, pDataArray->revPrimer, pDataArray->linker, pDataArray->spacer);  }
        
        TrimOligos* rtrimOligos = NULL;
        if (pDataArray->reorient) {
            //create reoriented primer and barcode pairs
            map<int, oligosPair> rpairedPrimers, rpairedBarcodes;
            for (map<int, oligosPair>::iterator it = pDataArray->pairedPrimers.begin(); it != pDataArray->pairedPrimers.end(); it++) {
                oligosPair tempPair(trimOligos->reverseOligo((it->second).reverse), (trimOligos->reverseOligo((it->second).forward))); //reversePrimer, rc ForwardPrimer
                rpairedPrimers[it->first] = tempPair;
            }
            for (map<int, oligosPair>::iterator it = pDataArray->pairedBarcodes.begin(); it != pDataArray->pairedBarcodes.end(); it++) {
                oligosPair tempPair(trimOligos->reverseOligo((it->second).reverse), (trimOligos->reverseOligo((it->second).forward))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[it->first] = tempPair;
            }
            
            int index = rpairedBarcodes.size();
            for (map<string, int>::iterator it = pDataArray->barcodes.begin(); it != pDataArray->barcodes.end(); it++) {
                oligosPair tempPair("", trimOligos->reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[index] = tempPair; index++;
            }
            
            index = rpairedPrimers.size();
            for (map<string, int>::iterator it = pDataArray->primers.begin(); it != pDataArray->primers.end(); it++) {
                oligosPair tempPair("", trimOligos->reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedPrimers[index] = tempPair; index++;
            }

            rtrimOligos = new TrimOligos(pDataArray->pdiffs, pDataArray->bdiffs, 0, 0, rpairedPrimers, rpairedBarcodes); numBarcodes = rpairedBarcodes.size();
        }
        
		pDataArray->count = 0;
		for(int i = 0; i < pDataArray->lineEnd; i++){ //end is the number of sequences to process
			           
			if (pDataArray->m->control_pressed) {
                delete trimOligos; if (pDataArray->reorient) { delete rtrimOligos; }
				inFASTA.close(); trimFASTAFile.close(); scrapFASTAFile.close();
				if ((pDataArray->createGroup) && (pDataArray->countfile == "")) {	 outGroupsFile.close();   }
                if(pDataArray->qFileName != "")	{	qFile.close();	scrapQualFile.close(); trimQualFile.close();	}
                if(pDataArray->nameFile != "")	{	scrapNameFile.close(); trimNameFile.close();	}
                if(pDataArray->countfile != "")	{	scrapCountFile.close(); trimCountFile.close();	}

				if(pDataArray->qFileName != ""){ qFile.close(); }
				return 0;
			}
			
			int success = 1;
			string trashCode = "";
            string commentString = "";
			int currentSeqsDiffs = 0;
            
			Sequence currSeq(inFASTA); pDataArray->m->gobble(inFASTA);
            Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
			
			QualityScores currQual; QualityScores savedQual;
			if(pDataArray->qFileName != ""){
				currQual = QualityScores(qFile);  pDataArray->m->gobble(qFile);
                savedQual.setName(currQual.getName()); savedQual.setScores(currQual.getScores());
			}
              
			
			string origSeq = currSeq.getUnaligned();
			if (origSeq != "") {
                pDataArray->count++;
				
				int barcodeIndex = 0;
				int primerIndex = 0;
				
                if(pDataArray->numLinkers != 0){
					success = trimOligos->stripLinker(currSeq, currQual);
					if(success > pDataArray->ldiffs)		{	trashCode += 'k';	}
					else{ currentSeqsDiffs += success;  }
				}
                
				if(numBarcodes != 0){
					vector<int> results = trimOligos->stripBarcode(currSeq, currQual, barcodeIndex);
                    if (pDataArray->pairedOligos) {
                        success = results[0] + results[2];
                        commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pDataArray->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], pDataArray->bdiffs) + ") ";
                    }
                    else {
                        success = results[0];
                        commentString += "bdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pDataArray->bdiffs) + ") ";
                    }

					if(success > pDataArray->bdiffs)		{	trashCode += 'b';	}
					else{ currentSeqsDiffs += success;  }
				}
                
                if(pDataArray->numSpacers != 0){
					success = trimOligos->stripSpacer(currSeq, currQual);
					if(success > pDataArray->sdiffs)		{	trashCode += 's';	}
					else{ currentSeqsDiffs += success;  }

				}
                
				if(pDataArray->numFPrimers != 0){
					vector<int> results = trimOligos->stripForward(currSeq, currQual, primerIndex, pDataArray->keepforward);
                    if (pDataArray->pairedOligos) {
                        success = results[0] + results[2];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pDataArray->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], pDataArray->pdiffs) + ") ";
                    }
                    else {
                        success = results[0];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pDataArray->pdiffs) + ") ";
                    }

					if(success > pDataArray->pdiffs)		{	trashCode += 'f';	}
					else{ currentSeqsDiffs += success;  }
				}
				
				if (currentSeqsDiffs > pDataArray->tdiffs)	{	trashCode += 't';   }
				
				if(pDataArray->numRPrimers != 0){
					vector<int> results = trimOligos->stripReverse(currSeq, currQual);
                    success = results[0];
                    commentString += "rpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pDataArray->pdiffs) + ") ";
                    if(success > pDataArray->pdiffs)		{	trashCode += 'r';	}
                    else{ currentSeqsDiffs += success;  }
				}
                
                if (pDataArray->reorient && (trashCode != "")) { //if you failed and want to check the reverse
                    int thisSuccess = 0;
                    string thisTrashCode = "";
                    string thiscommentString = "";
                    int thisCurrentSeqsDiffs = 0;
                    
                    int thisBarcodeIndex = 0;
                    int thisPrimerIndex = 0;
                    
                    if(numBarcodes != 0){
                        vector<int> results = rtrimOligos->stripBarcode(savedSeq, savedQual, thisBarcodeIndex);
                        if (pDataArray->pairedOligos) {
                            thisSuccess = results[0] + results[2];
                            thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pDataArray->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], pDataArray->bdiffs) + ") ";
                        }
                        else {
                            thisSuccess = results[0];
                            thiscommentString += "bdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pDataArray->bdiffs) + ") ";
                        }

                        if(thisSuccess > pDataArray->bdiffs)		{	thisTrashCode += 'b';	}
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    if(pDataArray->numFPrimers != 0){
                        vector<int> results = rtrimOligos->stripForward(savedSeq, savedQual, thisPrimerIndex, pDataArray->keepforward);
                        if (pDataArray->pairedOligos) {
                            thisSuccess = results[0] + results[2];
                            thiscommentString += "fpdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pDataArray->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], pDataArray->pdiffs) + ") ";
                        }
                        else {
                            thisSuccess = results[0];
                            thiscommentString += "pdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pDataArray->pdiffs) + ") ";
                        }

                        if(thisSuccess > pDataArray->pdiffs)		{	thisTrashCode += 'f';	}
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    if (thisCurrentSeqsDiffs > pDataArray->tdiffs)	{	thisTrashCode += 't';   }
                    
                    if (thisTrashCode == "") {
                        trashCode = thisTrashCode;
                        success = thisSuccess;
                        commentString = thiscommentString;
                        currentSeqsDiffs = thisCurrentSeqsDiffs;
                        barcodeIndex = thisBarcodeIndex;
                        primerIndex = thisPrimerIndex;
                        savedSeq.reverseComplement();
                        currSeq.setAligned(savedSeq.getAligned());
                        if(pDataArray->qFileName != ""){
                            savedQual.flipQScores();
                            currQual.setScores(savedQual.getScores());
                        }
                    }else { trashCode += "(" + thisTrashCode + ")";  }
                }

                
				if(pDataArray->keepFirst != 0){
					//success = keepFirstTrim(currSeq, currQual);
                    success = 1;
                    if(currQual.getName() != ""){
                        currQual.trimQScores(-1, pDataArray->keepFirst);
                    }
                    currSeq.trim(pDataArray->keepFirst);
				}
				
				if(pDataArray->removeLast != 0){
					//success = removeLastTrim(currSeq, currQual);
                    success = 0;
                    int length = currSeq.getNumBases() - pDataArray->removeLast;
                    
                    if(length > 0){
                        if(currQual.getName() != ""){
                            currQual.trimQScores(-1, length);
                        }
                        currSeq.trim(length);
                        success = 1;
                    }
                    else{ success = 0; }
                    
					if(!success)				{	trashCode += 'l';	}
				}
                
				
				if(pDataArray->qFileName != ""){
					int origLength = currSeq.getNumBases();
					
					if(pDataArray->qThreshold != 0)			{	success = currQual.stripQualThreshold(currSeq, pDataArray->qThreshold);			}
					else if(pDataArray->qAverage != 0)		{	success = currQual.cullQualAverage(currSeq, pDataArray->qAverage, pDataArray->logtransform);				}
					else if(pDataArray->qRollAverage != 0)	{	success = currQual.stripQualRollingAverage(currSeq, pDataArray->qRollAverage, pDataArray->logtransform);	}
					else if(pDataArray->qWindowAverage != 0){	success = currQual.stripQualWindowAverage(currSeq, pDataArray->qWindowStep, pDataArray->qWindowSize, pDataArray->qWindowAverage, pDataArray->logtransform);	}
					else						{	success = 1;				}
					
					//you don't want to trim, if it fails above then scrap it
					if ((!pDataArray->qtrim) && (origLength != currSeq.getNumBases())) {  success = 0; }
					
					if(!success)				{	trashCode += 'q';	}
				}				
                
				if(pDataArray->minLength > 0 || pDataArray->maxLength > 0){
					//success = cullLength(currSeq);
                    int length = currSeq.getNumBases();
                    success = 0;	//guilty until proven innocent
                    if(length >= pDataArray->minLength && pDataArray->maxLength == 0)			{	success = 1;	}
                    else if(length >= pDataArray->minLength && length <= pDataArray->maxLength)	{	success = 1;	}
                    else												{	success = 0;	}
                    
					if(!success)				{	trashCode += 'l';	}
				}
				if(pDataArray->maxHomoP > 0){
					//success = cullHomoP(currSeq);
                    int longHomoP = currSeq.getLongHomoPolymer();
                    success = 0;	//guilty until proven innocent
                    if(longHomoP <= pDataArray->maxHomoP){	success = 1;	}
                    else					{	success = 0;	}
                    
					if(!success)				{	trashCode += 'h';	}
				}
				if(pDataArray->maxAmbig != -1){
					//success = cullAmbigs(currSeq);
                    int numNs = currSeq.getAmbigBases();
                    success = 0;	//guilty until proven innocent
                    if(numNs <= pDataArray->maxAmbig)	{	success = 1;	}
                    else					{	success = 0;	}
					if(!success)				{	trashCode += 'n';	}
				}
				
				if(pDataArray->flip){		// should go last			
					currSeq.reverseComplement();
					if(pDataArray->qFileName != ""){
						currQual.flipQScores();	
					}
				}
				
                string seqComment = currSeq.getComment();
                currSeq.setComment("\t" + commentString + "\t" + seqComment);
                
				if(trashCode.length() == 0){
                    string thisGroup = "";
                    if (pDataArray->createGroup) {
						if(numBarcodes != 0){
							thisGroup = pDataArray->barcodeNameVector[barcodeIndex];
							if (pDataArray->numFPrimers != 0) {
								if (pDataArray->primerNameVector[primerIndex] != "") { 
									if(thisGroup != "") {
										thisGroup += "." + pDataArray->primerNameVector[primerIndex]; 
									}else {
										thisGroup = pDataArray->primerNameVector[primerIndex]; 
									}
								} 
							}
                        }
                    }
                    
                    int pos = thisGroup.find("ignore");
                    if (pos == string::npos) {
                        
                        currSeq.setAligned(currSeq.getUnaligned());
                        currSeq.printSequence(trimFASTAFile);
                        
                        if(pDataArray->qFileName != ""){
                            currQual.printQScores(trimQualFile);
                        }
                        
                        if(pDataArray->nameFile != ""){
                            map<string, string>::iterator itName = pDataArray->nameMap.find(currSeq.getName());
                            if (itName != pDataArray->nameMap.end()) {  trimNameFile << itName->first << '\t' << itName->second << endl; }
                            else { pDataArray->m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); }
                        }
                        
                        int numRedundants = 0;
                        if (pDataArray->countfile != "") {
                            map<string, int>::iterator itCount = pDataArray->nameCount.find(currSeq.getName());
                            if (itCount != pDataArray->nameCount.end()) { 
                                trimCountFile << itCount->first << '\t' << itCount->second << endl;
                                numRedundants = itCount->second-1;
                            }else { pDataArray->m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your count file, please correct."); pDataArray->m->mothurOutEndLine(); }
                        }
                        
                        if (pDataArray->createGroup) {
                            if(numBarcodes != 0){
                                
                                if (pDataArray->countfile == "") { outGroupsFile << currSeq.getName() << '\t' << thisGroup << endl; }
                                else {   pDataArray->groupMap[currSeq.getName()] = thisGroup; }
                                
                                if (pDataArray->nameFile != "") {
                                    map<string, string>::iterator itName = pDataArray->nameMap.find(currSeq.getName());
                                    if (itName != pDataArray->nameMap.end()) { 
                                        vector<string> thisSeqsNames; 
                                        pDataArray->m->splitAtChar(itName->second, thisSeqsNames, ',');
                                        numRedundants = thisSeqsNames.size()-1; //we already include ourselves below
                                        for (int k = 1; k < thisSeqsNames.size(); k++) { //start at 1 to skip self
                                            outGroupsFile << thisSeqsNames[k] << '\t' << thisGroup << endl;
                                        }
                                    }else { pDataArray->m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); }							
                                }
                                
                                map<string, int>::iterator it = pDataArray->groupCounts.find(thisGroup);
                                if (it == pDataArray->groupCounts.end()) {	pDataArray->groupCounts[thisGroup] = 1 + numRedundants; }
                                else { pDataArray->groupCounts[it->first] += (1 + numRedundants); }
                                
                            }
                        }
                        
                        if(pDataArray->allFiles){
                            ofstream output;
                            pDataArray->m->openOutputFileAppend(pDataArray->fastaFileNames[barcodeIndex][primerIndex], output);
                            currSeq.printSequence(output);
                            output.close();
                            
                            if(pDataArray->qFileName != ""){
                                pDataArray->m->openOutputFileAppend(pDataArray->qualFileNames[barcodeIndex][primerIndex], output);
                                currQual.printQScores(output);
                                output.close();							
                            }
                            
                            if(pDataArray->nameFile != ""){
                                map<string, string>::iterator itName = pDataArray->nameMap.find(currSeq.getName());
                                if (itName != pDataArray->nameMap.end()) { 
                                    pDataArray->m->openOutputFileAppend(pDataArray->nameFileNames[barcodeIndex][primerIndex], output);
                                    output << itName->first << '\t' << itName->second << endl; 
                                    output.close();
                                }else { pDataArray->m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); }
                            }
                        }
                    }
				}
				else{
					if(pDataArray->nameFile != ""){ //needs to be before the currSeq name is changed
						map<string, string>::iterator itName = pDataArray->nameMap.find(currSeq.getName());
						if (itName != pDataArray->nameMap.end()) {  scrapNameFile << itName->first << '\t' << itName->second << endl; }
						else { pDataArray->m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); }
					}
                    if (pDataArray->countfile != "") {
                        map<string, int>::iterator itCount = pDataArray->nameCount.find(currSeq.getName());
                        if (itCount != pDataArray->nameCount.end()) { 
                            trimCountFile << itCount->first << '\t' << itCount->second << endl;
                        }else { pDataArray->m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your count file, please correct."); pDataArray->m->mothurOutEndLine(); }
                    }
					currSeq.setName(currSeq.getName() + '|' + trashCode);
					currSeq.setUnaligned(origSeq);
					currSeq.setAligned(origSeq);
					currSeq.printSequence(scrapFASTAFile);
					if(pDataArray->qFileName != ""){
						currQual.printQScores(scrapQualFile);
					}
				}
				
			}
			
			//report progress
			if((pDataArray->count) % 1000 == 0){	pDataArray->m->mothurOut(toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
			
		}
		//report progress
		if((pDataArray->count) % 1000 != 0){	pDataArray->m->mothurOut(toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
		
        if (pDataArray->reorient) { delete rtrimOligos; }
		delete trimOligos;
		inFASTA.close();
		trimFASTAFile.close();
		scrapFASTAFile.close();
		if (pDataArray->createGroup) {	 outGroupsFile.close();   }
		if(pDataArray->qFileName != "")	{	qFile.close();	scrapQualFile.close(); trimQualFile.close();	}
		if(pDataArray->nameFile != "")	{	scrapNameFile.close(); trimNameFile.close();	}
		
        return 0;
            
        }
        catch(exception& e) {
            pDataArray->m->errorOut(e, "TrimSeqsCommand", "MyTrimThreadFunction");
            exit(1);
        }
    } 
#endif
    

/**************************************************************************************************/

#endif

