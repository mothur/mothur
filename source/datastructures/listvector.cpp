/*
 *  list.cpp
 *  
 *
 *  Created by Pat Schloss on 8/8/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


#include "sabundvector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"
#include "listvector.hpp"

//sorts highest to lowest
/***********************************************************************/
inline bool abundNamesSort(string left, string right){
    
    int countLeft = 0;
    if(left != ""){
        countLeft = 1;
        for(int i=0;i<left.size();i++){  if(left[i] == ','){  countLeft++;  }  }
    }
    
    int countRight = 0;
    if(right != ""){
        countRight = 1;
        for(int i=0;i<right.size();i++){  if(right[i] == ','){  countRight++;  }  }
    }
    
	if (countLeft > countRight) {
        return true;
    }
    return false;	
}
//sorts highest to lowest
/***********************************************************************/
inline bool abundNamesSort2(listCt left, listCt right){
    if (left.bin == "") { return false; }
    if (right.bin == "") { return true; }
    if (left.binSize > right.binSize) {
        return true;
    }
    return false;
} 


/***********************************************************************/

ListVector::ListVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0), otuTag("Otu"), printListHeaders(true) {}

/***********************************************************************/

ListVector::ListVector(int n):	DataVector(), data(n, "") , maxRank(0), numBins(0), numSeqs(0), otuTag("Otu"), printListHeaders(true){}

/***********************************************************************/

ListVector::ListVector(string id, vector<string> lv, string& tag) : DataVector(id), data(lv){
	try {
        printListHeaders = true;
        otuTag = tag;
        Utils util;
		for(int i=0;i<data.size();i++){
			if(data[i] != ""){
				int binSize = util.getNumNames(data[i]);
				numBins = i+1;
				if(binSize > maxRank)	{	maxRank = binSize;	}
				numSeqs += binSize;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "ListVector");
		exit(1);
	}
}

/**********************************************************************/

ListVector::ListVector(ifstream& f, string& readHeaders, string& labelTag) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
	try {
        printListHeaders = true;
		int thisNumBins = 0;
        Utils util;
        
        //are we at the beginning of the file??
		if (readHeaders == "") {
			f >> label;
            
			//is this a shared file that has headers
			if (label == "label") {
				
				//gets "numOtus"
				f >> label; util.gobble(f);
				
				//eat rest of line
				label = util.getline(f); util.gobble(f);
				
				//parse labels to save
				istringstream iStringStream(label);
				while(!iStringStream.eof()){
					if (m->getControl_pressed()) { break; }
					string temp;
					iStringStream >> temp;  util.gobble(iStringStream);
					binLabels.push_back(temp);
				}
                if (binLabels.size() != 0) {
                    string binLabelTag = binLabels[0];
                    labelTag = "";
                    for (int i = 0; i < binLabelTag.length(); i++) { if (isalpha(binLabelTag[i])){ labelTag += binLabelTag[i]; } }
                }
				f >> label >> thisNumBins;
			}else {
                //read in first row
                f >> thisNumBins;
                
                //make binlabels because we don't have any
                string snumBins = toString(thisNumBins);
                if (labelTag == "") { labelTag = "Otu"; }
                for (int i = 0; i < thisNumBins; i++) {
                    //if there is a bin label use it otherwise make one
                    string binLabel = labelTag;
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabels.push_back(binLabel);
                }
            }
		}else { f >> label >> thisNumBins; }
		
        util.gobble(f);
		data.assign(thisNumBins, "");
		string inputData = "";
        otuTag = labelTag;
	
        string buffer = util.getline(f); util.gobble(f);
        vector<string> thisListBins = util.splitWhiteSpace(buffer);
        
        if (thisListBins.size() != thisNumBins) {
            m->mothurOut("[ERROR]: Your list file indicates you have " + toString(thisNumBins) + ", but mothur found " + toString(thisListBins.size())+ " bins in your file, please correct.\n"); m->setControl_pressed(true);
        }else {
            for(int i=0;i<thisListBins.size();i++){ set(i, thisListBins[i]); }
        }
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "ListVector");
		exit(1);
	}
}

/***********************************************************************/

void ListVector::set(int binNumber, string seqNames){
	try {
        Utils util;
		int nNames_old = util.getNumNames(data[binNumber]);
		data[binNumber] = seqNames;
		int nNames_new = util.getNumNames(seqNames);
	
		if(nNames_old == 0)			{	numBins++;				}
		if(nNames_new == 0)			{	numBins--;				}
		if(nNames_new > maxRank)	{	maxRank = nNames_new;	}
	
		numSeqs += (nNames_new - nNames_old);
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "set");
		exit(1);
	}
}

/***********************************************************************/

string ListVector::get(int index){
    if (index < data.size()) {  return data[index];  }
    
	return "";
}
/***********************************************************************/

void ListVector::setLabels(vector<string> labels){
	try {
		binLabels = labels;
        getLabels();
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "setLabels");
		exit(1);
	}
}

/***********************************************************************/
//could potentially end up with duplicate binlabel names with code below.
//we don't currently use them in a way that would do that.
//if you had a listfile that had been subsampled and then added to it, dup names would be possible.
vector<string> ListVector::getLabels(){
    try {
        Utils util;
        util.getOTUNames(binLabels, numBins, otuTag);
        return binLabels;
    }
	catch(exception& e) {
		m->errorOut(e, "ListVector", "getLabels");
		exit(1);
	}
}
/***********************************************************************/
string ListVector::getOTUName(int bin){
    try {
        if (binLabels.size() > bin) {  }
        else { getLabels(); }
        return binLabels[bin];
    }
    catch(exception& e) {
        m->errorOut(e, "ListVector", "getOTUName");
        exit(1);
    }
}
/***********************************************************************/
int ListVector::getOTUTotal(string otuLabel){
    try {
        //find bin number
        int binNumber = -1;
        
        getLabels();
        
        for (int i = 0; i < binLabels.size(); i++) {
            if (util.getSimpleLabel(binLabels[i]) == util.getSimpleLabel(otuLabel)) {
                binNumber = i; break;
            }
        }
        
        if (binNumber == -1) { return 0; }
        
        return (util.getNumNames(get(binNumber)));
    }
    catch(exception& e) {
        m->errorOut(e, "ListVector", "getOTUTotal");
        exit(1);
    }
}

/***********************************************************************/

void ListVector::push_back(string seqNames){
	try {
        Utils util;
		data.push_back(seqNames);
		int nNames = util.getNumNames(seqNames);
	
		numBins++;
	
		if(nNames > maxRank)	{	maxRank = nNames;	}
	
		numSeqs += nNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "push_back");
		exit(1);
	}
}/***********************************************************************/
int ListVector::push_back(string bin, int nNames, string binLabel){
    try {
        if (binLabel == "") { //create one
            int otuNum = 1; bool notDone = true;
            
            //find label prefix
            string prefix = "Otu";
            if (binLabels.size() != 0) {
                if (binLabels[binLabels.size()-1][0] == 'P') { prefix = "PhyloType"; }
                
                string tempLabel = binLabels[binLabels.size()-1];
                string simpleLastLabel = util.getSimpleLabel(tempLabel);
                util.mothurConvert(simpleLastLabel, otuNum); otuNum++;
            }
            
            string potentialLabel = toString(otuNum);
            
            while (notDone) {
                if (m->getControl_pressed()) { notDone = false; break; }
                
                potentialLabel = toString(otuNum);
                vector<string>::iterator it = find(binLabels.begin(), binLabels.end(), potentialLabel);
                if (it == binLabels.end()) {
                    potentialLabel = prefix + toString(otuNum);
                    it = find(binLabels.begin(), binLabels.end(), potentialLabel);
                    if (it == binLabels.end()) {
                        notDone = false; break;
                    }
                }
                otuNum++;
            }
            
            binLabel = potentialLabel;
        }
        binLabels.push_back(binLabel);
        
        data.push_back(bin);
        numBins++;
        
        if(nNames > maxRank)	{	maxRank = nNames;	}
        
        numSeqs += nNames;
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ListVector", "push_back");
        exit(1);
    }
}
/***********************************************************************/

void ListVector::resize(int size){
	data.resize(size);		
}

/***********************************************************************/

int ListVector::size(){
	return data.size();
}
/***********************************************************************/

void ListVector::clear(){
	numBins = 0;
	maxRank = 0;
	numSeqs = 0;
	return data.clear();
	
}

/***********************************************************************/
void ListVector::printHeaders(ostream& output, map<string, int>& ct, bool sortPlease){
	try {
        if (printListHeaders) {
            if (binLabels.size() == 0) { sortPlease = false; } //we are creating arbitary otuNames
            
            vector<string> theseLabels = getLabels();
            
            output << "label\tnum" + otuTag + "s";
            
            if (sortPlease) {
                Utils util;
                vector<listCt> hold;
                for (int i = 0; i < data.size(); i++) {
                    if (data[i] != "") {
                        vector<string> binNames;
                        string bin = data[i];
                        util.splitAtComma(bin, binNames);
                        int total = 0;
                        for (int j = 0; j < binNames.size(); j++) {
                            map<string, int>::iterator it = ct.find(binNames[j]);
                            if (it == ct.end()) {
                                m->mothurOut("[ERROR]: " + binNames[j] + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
                            }else { total += ct[binNames[j]]; }
                        }
                        listCt temp(data[i], total, theseLabels[i]);
                        hold.push_back(temp);
                    }
                }
                sort(hold.begin(), hold.end(), abundNamesSort2);

                //print original label for sorted by abundance otu
                for (int i = 0; i < hold.size(); i++) { output  << '\t' << hold[i].label; }
            }else {
                for (int i = 0; i < theseLabels.size(); i++) { output  << '\t' << theseLabels[i]; }
            }
            
            output << endl;
            
            printListHeaders = false;
        }
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "printHeaders");
		exit(1);
	}
}

/***********************************************************************/

void ListVector::print(ostream& output, map<string, int>& ct){
	try {
        printHeaders(output, ct, true);
		output << label << '\t' << numBins;
        
        Utils util;
        vector<listCt> hold;
        for (int i = 0; i < data.size(); i++) {
            if (data[i] != "") {
                vector<string> binNames;
                string bin = data[i];
                util.splitAtComma(bin, binNames);
                int total = 0;
                for (int j = 0; j < binNames.size(); j++) {
                    map<string, int>::iterator it = ct.find(binNames[j]);
                    if (it == ct.end()) {
                        m->mothurOut("[ERROR]: " + binNames[j] + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
                    }else { total += ct[binNames[j]]; }
                }
                listCt temp(data[i], total, "");
                hold.push_back(temp);
            }
        }
        sort(hold.begin(), hold.end(), abundNamesSort2);
        
        for(int i=0;i<hold.size();i++){
            if(hold[i].bin != ""){
                output  << '\t' << hold[i].bin;
            }
        }
        output << endl;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "print");
		exit(1);
	}
}

/***********************************************************************/

void ListVector::print(ostream& output){
    try {
        map<string, int> ct;
        for (int i = 0; i < data.size(); i++) {
            if (data[i] != "") {
                string bin = data[i];
                vector<string> binNames; util.splitAtComma(bin, binNames);
                for (int j = 0; j < binNames.size(); j++) { ct[binNames[j]] = 1; }
            }
        }
        
        print(output, ct);
    }
    catch(exception& e) {
        m->errorOut(e, "ListVector", "print");
        exit(1);
    }
}
/***********************************************************************/
//no sort for subsampling and get.otus and remove.otus
void ListVector::print(ostream& output, bool sortOtus){
    try {
        map<string, int> ct;
        for (int i = 0; i < data.size(); i++) {
            if (data[i] != "") {
                string bin = data[i];
                vector<string> binNames; util.splitAtComma(bin, binNames);
                for (int j = 0; j < binNames.size(); j++) { ct[binNames[j]] = 1; }
            }
        }
        
        printHeaders(output, ct, sortOtus);
        output << label << '\t' << numBins;
        
        vector<string> hold = data;
        if (sortOtus) { sort(hold.begin(), hold.end(), abundNamesSort); }
        
        //find first non blank otu
        int start = 0;
        for(int i=0;i<hold.size();i++){  if(hold[i] != ""){  start = i; break; } }
        
        for(int i=start;i<hold.size();i++){
            if(hold[i] != ""){
                output << '\t' << hold[i];
            }
        }
        output << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "ListVector", "print");
        exit(1);
    }
}

/***********************************************************************/

RAbundVector ListVector::getRAbundVector(){
	try {
		RAbundVector rav;
        Utils util;
		for(int i=0;i<data.size();i++){
			int binSize = util.getNumNames(data[i]);
			rav.push_back(binSize);
		}
		rav.setLabel(label);
	
		return rav;
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "getRAbundVector");
		exit(1);
	}
}

/***********************************************************************/

SAbundVector ListVector::getSAbundVector(){
	try {
		SAbundVector sav(maxRank+1);
        Utils util;
		for(int i=0;i<data.size();i++){
			int binSize = util.getNumNames(data[i]);	
			sav.set(binSize, sav.get(binSize) + 1);	
		}
		sav.set(0, 0);
		sav.setLabel(label);
	
		return sav;
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "getSAbundVector");
		exit(1);
	}
}

/***********************************************************************/

OrderVector ListVector::getOrderVector(map<string,int>* orderMap = NULL){
	try {
        Utils util;
		if(orderMap == NULL){
			
            vector<int> ovData;
			for(int i=0;i<data.size();i++){
				int binSize = util.getNumNames(data[i]);		
				for(int j=0;j<binSize;j++){
					ovData.push_back(i);
				}
			}
			util.mothurRandomShuffle(ovData);
            OrderVector ov;
            for(int i=0;i<ovData.size();i++){ ov.push_back(ovData[i]); }
			ov.setLabel(label);
			ov.getNumBins();
		
			return ov;
		
		}
		else{
			OrderVector ov(numSeqs);
		
			for(int i=0;i<data.size();i++){
				string listOTU = data[i];
				int length = listOTU.size();
				
				string seqName="";
			
				for(int j=0;j<length;j++){
				
					if(listOTU[j] != ','){
						seqName += listOTU[j];
					}
					else{
						if(orderMap->count(seqName) == 0){
							m->mothurOut(seqName + " not found, check *.names file\n");
							exit(1);
						}
					
						ov.set((*orderMap)[seqName], i);
						seqName = "";
					}						
				}
			
				if(orderMap->count(seqName) == 0){
					m->mothurOut(seqName + " not found, check *.names file\n");
					exit(1);
				}
				ov.set((*orderMap)[seqName], i);	
			}
		
			ov.setLabel(label);
			ov.getNumBins();
		
			return ov;		
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ListVector", "getOrderVector");
		exit(1);
	}
}

/***********************************************************************/
