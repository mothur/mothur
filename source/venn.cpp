/*
 *  venn.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "venn.h"
#include "ace.h"
#include "sobs.h"
#include "chao1.h"
#include "nseqs.h"
#include "sharedchao1.h"
#include "sharedsobscollectsummary.h"


//**********************************************************************************************************************
Venn::Venn(string o, bool n, string f, int fs, bool so) : outputDir(o), nseqs(n), inputfile(f), fontSize(fs), sharedOtus(so) {
	try {
		m = MothurOut::getInstance();
	}
	catch(exception& e) {
		m->errorOut(e, "Venn", "Venn");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> Venn::getPic(SAbundVector* sabund, vector<Calculator*> vCalcs) {
	try {
		
		vector<string> outputNames;

		for(int i=0;i<vCalcs.size();i++){
			string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + "." + sabund->getLabel() + "." + vCalcs[i]->getName() + ".svg";
			outputNames.push_back(filenamesvg);
			util.openOutputFile(filenamesvg, outsvg);
			
			if (m->getControl_pressed()) { outsvg.close(); return outputNames; }

			vector<double> data = vCalcs[i]->getValues(sabund);
			
			int width = 1500;
			int height = 1500;
			
			//svg image
			outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " +  toString(width) + " " + toString(height) + " \" >\n";
			outsvg << "<g>\n";
				
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" +  toString(width) +  "\" height=\"" +  toString(height) +  "\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.40 * width)) +  "\" y=\"" +  toString(int(0.05 * height)) +  "\">Venn Diagram at distance " + sabund->getLabel() + "</text>\n";
			outsvg << "<circle fill=\"red\" opacity=\".5\" stroke=\"black\" cx=\"" +  toString(int(0.50 * width)) +  "\" cy=\"" +  toString(int(0.29 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString((width / 2) - ((int)toString(data[0]).length() / 2)) + "\" y=\"" +  toString(int(0.28 * height)) +  "\">" + toString(data[0]) + "</text>\n";  
			
			if (data.size() == 3) { 
				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.35 * width)) +  "\" y=\"" +  toString(int(0.60 * height)) +  "\">The lower bound of the confidence interval is " + toString(data[1]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.35 * width)) +  "\" y=\"" +  toString(int(0.645 * height)) +  "\">The upper bound of the confidence interval is " + toString(data[2]) + "</text>\n";
			}
			
			if (nseqs) {
				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.35 * width)) +  "\" y=\"" +  toString(int(0.70 * height)) +  "\">The number of sequences represented is " + toString(sabund->getNumSeqs()) + "</text>\n";
			}
			
			outsvg << "</g>\n</svg>\n";
			outsvg.close();
		}
		
		return outputNames;
	}
	catch(exception& e) {
		m->errorOut(e, "Venn", "getPic");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> Venn::getPic(vector<SharedRAbundVector*> lookup, vector<Calculator*> vCalcs) {
	try {

		vector<SharedRAbundVector*> subset;
		vector<string> outputNames;
		
		int width = 1500;
		int height = 1500;
		
		/******************* 1 Group **************************/
		if (lookup.size() == 1) {
					
			SAbundVector s;
			s = lookup[0]->getSAbundVector();  SAbundVector* sabund = &s;
			
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
				string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + lookup[0]->getGroup() + ".svg";
				outputNames.push_back(filenamesvg);
				util.openOutputFile(filenamesvg, outsvg);
				
				if (m->getControl_pressed()) { outsvg.close(); return outputNames; }
				
				//in essence you want to run it like a single 
				if (vCalcs[i]->getName() == "sharedsobs") {
					singleCalc = new Sobs();
				}else if (vCalcs[i]->getName() == "sharedchao") {
					singleCalc = new Chao1();
				}else if (vCalcs[i]->getName() == "sharedace") {
					singleCalc = new Ace(10);
				}
				
				vector<double> data = singleCalc->getValues(sabund);
			
				//svg image
				outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " +  toString(width) + " " + toString(height) + " \">\n";
				outsvg << "<g>\n";
				
				outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" +  toString(width) +  "\" height=\"" +  toString(height) +  "\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.40 * width)) +  "\" y=\"" +  toString(int(0.05 * height)) +  "\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
				outsvg << "<circle fill=\"red\" opacity=\".5\" stroke=\"black\" cx=\"" +  toString(int(0.50 * width)) +  "\" cy=\"" +  toString(int(0.29 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.50 * width) - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"" +  toString(int(0.24 * height)) +  "\">" + lookup[0]->getGroup() + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.50 * width) - ((int)toString(data[0]).length() / 2)) + "\" y=\"" +  toString(int(0.28 * height)) +  "\">" + toString(data[0]) + "</text>\n";  
			
				if (data.size() == 3) { 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.35 * width)) +  "\" y=\"" +  toString(int(0.60 * height)) +  "\" >The lower bound of the confidence interval is " + toString(data[1]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.35 * width)) +  "\" y=\"" +  toString(int(0.645 * height)) +  "\">The upper bound of the confidence interval is " + toString(data[2]) + "</text>\n";
				}
				
				if (nseqs) {
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.35 * width)) +  "\" y=\"" +  toString(int(0.70 * height)) +  "\">The number of sequences represented is " + toString(sabund->getNumSeqs()) + "</text>\n";
				}
				
				outsvg << "</g>\n</svg>\n";
				outsvg.close();
				delete singleCalc;
				
			}
		/******************* 2 Groups **************************/	
		
		}else if (lookup.size() == 2) {
			//get sabund vector pointers so you can use the single calculators
			//one for each group
			SAbundVector sA, sB;
			SAbundVector* sabundA; SAbundVector* sabundB;
			sabundA = new SAbundVector(lookup[0]->getSAbundVector());//  sabundA = &sA;
			sabundB = new SAbundVector(lookup[1]->getSAbundVector());//  sabundB = &sB;
			
			subset.clear();
			subset.push_back(lookup[0]); subset.push_back(lookup[1]);
			
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
				string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + lookup[0]->getGroup() + "-" + lookup[1]->getGroup() + ".svg";

				outputNames.push_back(filenamesvg);
				util.openOutputFile(filenamesvg, outsvg);
				
				if (m->getControl_pressed()) { outsvg.close(); return outputNames; }
				
				//get estimates for sharedAB
                vector<string> labels;
				vector<double> shared = vCalcs[i]->getValues(subset, labels);
				
				//in essence you want to run it like a single 
				if (vCalcs[i]->getName() == "sharedsobs") {
					singleCalc = new Sobs();
                    if (sharedOtus &&  (labels.size() != 0)) {
                        string groupsTag = "";
                        for (int h = 0; h < lookup.size()-1; h++) { groupsTag += lookup[h]->getGroup() + "-"; }  groupsTag += lookup[lookup.size()-1]->getGroup();
                        string filenameShared = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + groupsTag + ".sharedotus";
                        
                        outputNames.push_back(filenameShared);
                        ofstream outShared;
                        util.openOutputFile(filenameShared, outShared);
                        outShared << "Groups\tNumShared\tOTULabels\n";
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup() << '\t' << labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared  << endl;
                        outShared.close();
                    }
				}else if (vCalcs[i]->getName() == "sharedchao") {
					singleCalc = new Chao1();
				}
				
				int sharedVal, numSeqsA, numSeqsB, uniqSeqsToA, uniqSeqsToB;
				if (nseqs) {
					NSeqs* nseqsCalc = new NSeqs();
					vector<double> data = nseqsCalc->getValues(lookup);
					sharedVal = data[0] + data[1];
					numSeqsA = sabundA->getNumSeqs();
					numSeqsB = sabundB->getNumSeqs();
					uniqSeqsToA = numSeqsA-data[0];
					uniqSeqsToB = numSeqsB-data[1];
					
					delete nseqsCalc;
				}
				
				
				//get estimates for numA
				vector<double> numA = singleCalc->getValues(sabundA);

				//get estimates for numB
				vector<double> numB = singleCalc->getValues(sabundB);
						
				//image window
				outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " +  toString(width) + " " + toString(height) + " \" >\n";
				outsvg << "<g>\n";

				//draw circles
				outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" +  toString(width) +  "\" height=\"" +  toString(height) +  "\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.40 * width)) +  "\" y=\"" +  toString(int(0.05 * height)) +  "\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
				outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.36 * width)) +  "\" cy=\"" +  toString(int(0.29 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
				outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.62 * width)) +  "\" cy=\"" +  toString(int(0.29 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.29 * width) - ((int)toString(numA[0]).length() / 2)) + "\" y=\"" +  toString(int(0.28 * height)) +  "\">" + toString(numA[0] - shared[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.7 * width) - ((int)toString(numB[0]).length() / 2)) + "\" y=\"" +  toString(int(0.28 * height)) +  "\">" + toString(numB[0] - shared[0]) + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.29 * width) - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"" +  toString(int(0.25 * height)) +  "\">" + lookup[0]->getGroup() + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.7 * width) - ((int)lookup[1]->getGroup().length() / 2)) + "\" y=\"" +  toString(int(0.25 * height)) +  "\">" + lookup[1]->getGroup() + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(shared[0]).length() / 2)) + "\" y=\"" +  toString(int(0.28 * height)) +  "\">" + toString(shared[0]) + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\"   font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.66 * height)) +  "\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA[0]);
				if (numA.size() == 3) { 
					outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]);
				}
				if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsA) + "; " + toString(uniqSeqsToA) + " sequences are not shared";  }  
				outsvg << "</text>\n";
		
				outsvg << "<text fill=\"black\" class=\"seri\"   font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.69 * height)) +  "\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB[0]);
				if (numB.size() == 3) { 
					outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]);
				}
				if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsB) + "; " + toString(uniqSeqsToB) + " sequences are not shared";  }  
				outsvg << "</text>\n";

				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.72 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(shared[0]);
				if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedVal) + "; " + toString((sharedVal / (float)(numSeqsA + numSeqsB))*100) + "% of these sequences are shared";  }  
				outsvg << "</text>\n";
				
				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.75 * height)) +  "\">Percentage of species that are shared in groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString((shared[0] / (float)(numA[0] + numB[0] - shared[0]))*100) + "</text>\n";
				
				outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.78 * height)) +  "\">The total richness for all groups is " + toString((float)(numA[0] + numB[0] - shared[0]))+ "</text>\n";;
				
				
				//close file
				outsvg << "</g>\n</svg>\n";
				outsvg.close();
				delete singleCalc;
			}
		/******************* 3 Groups **************************/
						
		}else if (lookup.size() == 3) {
			
			height = 1600;
			int windowSize = height;
	
			
			//get sabund vector pointers so you can use the single calculators
			//one for each group
			SAbundVector sA, sB, sC;
			SAbundVector* sabundA; SAbundVector* sabundB; SAbundVector* sabundC;
			sA = lookup[0]->getSAbundVector();  sabundA = &sA;
			sB = lookup[1]->getSAbundVector();  sabundB = &sB;
			sC = lookup[2]->getSAbundVector();  sabundC = &sC;
		
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
			
				string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + lookup[0]->getGroup() + "-" + lookup[1]->getGroup() + "-" + lookup[2]->getGroup() + ".svg";

				outputNames.push_back(filenamesvg);
				util.openOutputFile(filenamesvg, outsvg);
				
				if (m->getControl_pressed()) { outsvg.close(); return outputNames; }
				
				int sharedVal, sharedABVal, sharedACVal, sharedBCVal, numSeqsA, numSeqsB, numSeqsC, uniqSeqsToA, uniqSeqsToB, uniqSeqsToC;
				
				if (nseqs) {
					NSeqs* nseqsCalc = new NSeqs();
					vector<double> sharedData = nseqsCalc->getValues(lookup);
						
					vector<SharedRAbundVector*> mysubset; mysubset.push_back(lookup[0]); mysubset.push_back(lookup[1]);
					vector<double> sharedAB = nseqsCalc->getValues(mysubset);
						
					mysubset.clear(); mysubset.push_back(lookup[0]); mysubset.push_back(lookup[2]);
					vector<double> sharedAC = nseqsCalc->getValues(mysubset);
						
					mysubset.clear(); mysubset.push_back(lookup[1]); mysubset.push_back(lookup[2]);
					vector<double> sharedBC = nseqsCalc->getValues(mysubset);
						
					sharedVal = sharedData[0] + sharedData[1] + sharedData[2];
					sharedABVal = sharedAB[0] + sharedAB[1];
					sharedACVal = sharedAC[0] + sharedAC[1];
					sharedBCVal = sharedBC[0] + sharedBC[1];
					numSeqsA = sabundA->getNumSeqs();
					numSeqsB = sabundB->getNumSeqs();
					numSeqsC = sabundC->getNumSeqs();
					uniqSeqsToA = numSeqsA-sharedData[0];
					uniqSeqsToB = numSeqsC-sharedData[1];
					uniqSeqsToC = numSeqsB-sharedData[1];

					delete nseqsCalc;
				}

				
				if (vCalcs[i]->getName() == "sharedace") {
				
					singleCalc = new Ace(10);
					
					//get estimates for numA
					vector<double> numA = singleCalc->getValues(sabundA);
 			
					//get estimates for numB
					vector<double> numB = singleCalc->getValues(sabundB);
 				
					//get estimates for numC
					vector<double> numC = singleCalc->getValues(sabundC);


					//get estimates for sharedAB, sharedAC and sharedBC
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]);
					vector<double> sharedAB = vCalcs[i]->getValues(subset);
					
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]);
					vector<double> sharedAC = vCalcs[i]->getValues(subset);
					
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					vector<double> sharedBC = vCalcs[i]->getValues(subset);
					
					vector<double> sharedAwithBC;
					vector<double> sharedBwithAC;
					vector<double> sharedCwithAB;
					
					//find possible sharedABC values
					float sharedABC1 = 0.0; float sharedABC2 = 0.0; float sharedABC3 = 0.0; float sharedABC = 0.0;

					if (vCalcs[i]->getMultiple() == false) {
						//merge BC and estimate with shared with A
						SharedRAbundVector* merge = new SharedRAbundVector();
						for (int j = 0; j < lookup[1]->size(); j++) {  merge->push_back((lookup[1]->get(j) + lookup[2]->get(j)));  }
					
						subset.clear();
						subset.push_back(lookup[0]); subset.push_back(merge);
						sharedAwithBC = vCalcs[i]->getValues(subset);
				
						delete merge;
						//merge AC and estimate with shared with B
						merge = new SharedRAbundVector();
						for (int j = 0; j < lookup[0]->size(); j++) { merge->push_back((lookup[0]->get(j) + lookup[2]->get(j))); }
					
						subset.clear();
						subset.push_back(merge); subset.push_back(lookup[1]);
						sharedBwithAC = vCalcs[i]->getValues(subset);
				
						delete merge;
						//merge AB and estimate with shared with C
						merge = new SharedRAbundVector();
						for (int j = 0; j < lookup[0]->size(); j++) { merge->push_back((lookup[0]->get(j) + lookup[1]->get(j))); }
					
						subset.clear();
						subset.push_back(lookup[2]); subset.push_back(merge);
						sharedCwithAB = vCalcs[i]->getValues(subset);
						delete merge;
					
						sharedABC1 = sharedAB[0] + sharedAC[0] - sharedAwithBC[0];
						sharedABC2 = sharedAB[0] + sharedBC[0] - sharedBwithAC[0];
						sharedABC3 = sharedAC[0] + sharedBC[0] - sharedCwithAB[0];
	 
						//if any of the possible m's are - throw them out
						if (sharedABC1 < 0.00001) { sharedABC1 = 0; }
						if (sharedABC2 < 0.00001) { sharedABC2 = 0; }
						if (sharedABC3 < 0.00001) { sharedABC3 = 0; }
			
						//sharedABC is the minimum of the 3 possibilities
						if ((sharedABC1 < sharedABC2) && (sharedABC1 < sharedABC3)) { sharedABC = sharedABC1; }
						else if ((sharedABC2 < sharedABC1) && (sharedABC2 < sharedABC3)) { sharedABC = sharedABC2; }
						else if ((sharedABC3 < sharedABC1) && (sharedABC3 < sharedABC2)) { sharedABC = sharedABC3; }	
					}else{
						vector<double> data = vCalcs[i]->getValues(lookup);
						sharedABC = data[0];
						sharedAwithBC.push_back(sharedAB[0] + sharedAC[0] - sharedABC);
						sharedBwithAC.push_back(sharedAB[0] + sharedBC[0] - sharedABC);
						sharedCwithAB.push_back(sharedAC[0] + sharedBC[0] - sharedABC);
					}
					
					//image window
					outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " +  toString(width) + " " + toString(windowSize) + " \" >\n";
					outsvg << "<g>\n";

					//draw circles
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" +  toString(width) +  "\" height=\"" +  toString(height) +  "\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.40 * width)) +  "\" y=\"" +  toString(int(0.44 * height)) +  "\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
					outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.33 * width)) +  "\" cy=\"" +  toString(int(0.25 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
					outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.65 * width)) +  "\" cy=\"" +  toString(int(0.25 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
					outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.5 * width)) +  "\" cy=\"" +  toString(int(0.5 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 

					//place labels within overlaps
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.33 * width) - ((int)toString(numA[0]-sharedAwithBC[0]).length() / 2)) + "\" y=\"" +  toString(int(0.22 * height)) +  "\">" + toString(numA[0]-sharedAwithBC[0]) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.33 * width) - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"" +  toString(int(0.19 * height)) +  "\">" + lookup[0]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedAB[0] - sharedABC).length() / 2)) + "\"  y=\"" +  toString(int(0.22 * height)) +  "\">" + toString(sharedAB[0] - sharedABC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.7 * width) - ((int)toString(numB[0]-sharedBwithAC[0]).length() / 2)) + "\"  y=\"" +  toString(int(0.22 * height)) +  "\">" + toString(numB[0]-sharedBwithAC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.7 * width) - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"" +  toString(int(0.19 * height)) +  "\">" + lookup[1]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.38 * width) - ((int)toString(sharedAC[0] - sharedABC).length() / 2)) + "\"  y=\"" +  toString(int(0.38 * height)) +  "\">" + toString(sharedAC[0] - sharedABC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.5 * width) - ((int)toString(numC[0]-sharedCwithAB[0]).length() / 2)) + "\"   y=\"" +  toString(int(0.54 * height)) +  "\">" + toString(numC[0]-sharedCwithAB[0]) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.5 * width) - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"" +  toString(int(0.52 * height)) +  "\">" + lookup[2]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.58 * width) - ((int)toString(sharedBC[0] - sharedABC).length() / 2)) + "\" y=\"" +  toString(int(0.38 * height)) +  "\">" + toString(sharedBC[0] - sharedABC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"" +  toString(int(0.34 * height)) +  "\">" + toString(sharedABC) + "</text>\n"; 
				
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.825 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(sharedAB[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedABVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.85 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedAC[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedACVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.875 * height)) +  "\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedBC[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedBCVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.9 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and combined groups " + lookup[1]->getGroup() + lookup[2]->getGroup() + " is " + toString(sharedAwithBC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.925 * height)) +  "\">The number of species shared between groups " + lookup[1]->getGroup() + " and combined groups " + lookup[0]->getGroup() + lookup[2]->getGroup() + " is " + toString(sharedBwithAC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.95 * height)) +  "\">The number of species shared between groups " + lookup[2]->getGroup() + " and combined groups " + lookup[0]->getGroup() + lookup[1]->getGroup() + " is " + toString(sharedCwithAB[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.725 * height)) +  "\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA[0]);
					if (numA.size() == 3) { 
						outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]);
					} 
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsA) + "; " + toString(uniqSeqsToA) + " sequences are not shared";  }  
					outsvg << "</text>\n";
			
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.75 * height)) +  "\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB[0]);
					if (numB.size() == 3) { 
						outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]);
					}
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsB) + "; " + toString(uniqSeqsToB) + " sequences are not shared";  }  
					outsvg << "</text>\n";
					
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.775 * height)) +  "\">The number of species in group " + lookup[2]->getGroup() + " is " + toString(numC[0]);
					if (numC.size() == 3) { 
						outsvg << " the lci is " + toString(numC[1]) + " and the hci is " + toString(numC[2]);
					}
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsC) + "; " + toString(uniqSeqsToC) + " sequences are not shared";  }  
					outsvg << "</text>\n";

					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.80 * height)) +  "\">The total richness of all the groups is " + toString(numA[0] + numB[0] + numC[0] - sharedAB[0] - sharedAC[0] - sharedBC[0] + sharedABC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.975 * height)) +  "\">The total shared richness is " + toString(sharedABC);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedVal);  }  
					outsvg << "</text>\n";
					
					delete singleCalc;
					
				}else { //sharedchao and sharedsobs are multigroup
					
                    ofstream outShared;
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs")) {
                        string groupsTag = "";
                        for (int h = 0; h < lookup.size()-1; h++) { groupsTag += lookup[h]->getGroup() + "-"; }  groupsTag += lookup[lookup.size()-1]->getGroup();
                        string filenameShared = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + groupsTag + ".sharedotus";
                        
                        outputNames.push_back(filenameShared);
                       
                        util.openOutputFile(filenameShared, outShared);
                        outShared << "Groups\tNumShared\tOTULabels\n";
                    }
					vector<SharedRAbundVector*> subset;

					//get estimates for numA
					subset.push_back(lookup[0]);
					vector<double> numA = vCalcs[i]->getValues(subset);
 			
					//get estimates for numB
					subset.clear();
					subset.push_back(lookup[1]);
					vector<double> numB = vCalcs[i]->getValues(subset);
 				
					//get estimates for numC
					subset.clear();
					subset.push_back(lookup[2]);
					vector<double> numC = vCalcs[i]->getValues(subset);

					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]);
                    vector<string> labels;
					vector<double> sharedab =  vCalcs[i]->getValues(subset, labels);
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup() << '\t' << labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
					
					subset.clear(); 
					subset.push_back(lookup[0]); subset.push_back(lookup[2]);
					vector<double> sharedac =  vCalcs[i]->getValues(subset, labels);
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[2]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
					
					subset.clear(); 
					subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					vector<double> sharedbc =  vCalcs[i]->getValues(subset, labels);
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[1]->getGroup() + "-" + lookup[2]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }

					
					subset.clear(); 
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					vector<double> sharedabc =  vCalcs[i]->getValues(subset, labels);
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup() + "-" + lookup[2]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                        outShared.close();
                    }
					
					//image window
					outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " +  toString(width) + " " + toString(windowSize) + " \" >\n";
					outsvg << "<g>\n";

					//draw circles
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" +  toString(width) +  "\" height=\"" +  toString(height) +  "\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.40 * width)) +  "\" y=\"" +  toString(int(0.05 * height)) +  "\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
					outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.33 * width)) +  "\" cy=\"" +  toString(int(0.25 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
					outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.65 * width)) +  "\" cy=\"" +  toString(int(0.25 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 
					outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"" +  toString(int(0.5 * width)) +  "\" cy=\"" +  toString(int(0.5 * height)) +  "\" r=\"" +  toString(int(0.22 * width)) +  "\"/>"; 

					//place labels within overlaps
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.29 * width) - ((int)toString(numA[0]-sharedab[0]-sharedac[0]+sharedabc[0]).length() / 2)) + "\" y=\"" +  toString(int(0.22 * height)) +  "\">" + toString(numA[0]-sharedab[0]-sharedac[0]+sharedabc[0]) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.29 * width) - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"" +  toString(int(0.19 * height)) +  "\">" + lookup[0]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedab[0] - sharedabc[0]).length() / 2)) + "\"  y=\"" +  toString(int(0.22 * height)) +  "\">" + toString(sharedab[0] - sharedabc[0]) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.68 * width) - ((int)toString(numB[0]-sharedab[0]-sharedbc[0]+sharedabc[0]).length() / 2)) + "\"  y=\"" +  toString(int(0.22 * height)) +  "\">" + toString(numB[0]-sharedab[0]-sharedbc[0]+sharedabc[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.68 * width) - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"" +  toString(int(0.19 * height)) +  "\">" + lookup[1]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.38 * width) - ((int)toString(sharedac[0] - sharedabc[0]).length() / 2)) + "\"  y=\"" +  toString(int(0.38 * height)) +  "\">" + toString(sharedac[0] - sharedabc[0]) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.5 * width) - ((int)toString(numC[0]-sharedac[0]-sharedbc[0]+sharedabc[0]).length() / 2)) + "\"   y=\"" +  toString(int(0.54 * height)) +  "\">" + toString(numC[0]-sharedac[0]-sharedbc[0]+sharedabc[0]) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.5 * width) - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"" +  toString(int(0.51 * height)) +  "\">" + lookup[2]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.59 * width) - ((int)toString(sharedbc[0] - sharedabc[0]).length() / 2)) + "\" y=\"" +  toString(int(0.38 * height)) +  "\">" + toString(sharedbc[0] - sharedabc[0]) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedabc[0]).length() / 2)) + "\"  y=\"" +  toString(int(0.35 * height)) +  "\">" + toString(sharedabc[0]) + "</text>\n"; 
				
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.825 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(sharedab[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedABVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.85 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedac[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedACVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.875 * height)) +  "\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedbc[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedBCVal);  }  
					outsvg << "</text>\n";
					
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.725 * height)) +  "\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA[0]);
					if (numA.size() == 3) { 
						outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]);
					}
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsA);  }  
					outsvg << "</text>\n";
			
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.75 * height)) +  "\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB[0]);
					if (numB.size() == 3) { 
						outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]);
					}
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsB);  }  
					outsvg << "</text>\n";
										
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.775 * height)) +  "\">The number of species in group " + lookup[2]->getGroup() + " is " + toString(numC[0]);
					if (numC.size() == 3) { 
						outsvg << " the lci is " + toString(numC[1]) + " and the hci is " + toString(numC[2]);
					}
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsC);  }  
					outsvg << "</text>\n";

					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.8 * height)) +  "\">The total richness of all the groups is " + toString(numA[0] + numB[0] + numC[0] - sharedab[0] - sharedac[0] - sharedbc[0] + sharedabc[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.9 * height)) +  "\">The total shared richness is " + toString(sharedabc[0]);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedVal);  }  
					outsvg << "</text>\n";

				}
								
				//close file
				outsvg << "</g>\n</svg>\n";
				outsvg.close();
				

			}
			
		/******************* 4 Groups **************************/
		
		}else if (lookup.size() == 4) {
			
			height = 1600;
			
			int windowSize = height;
		
			//calc the shared otu
			float sharedABCD = 0;
			float numA = 0; float numB = 0; float numC = 0; float numD = 0;
			float sharedAB = 0; float sharedAC = 0; float sharedBC = 0; float sharedAD = 0; float sharedBD = 0; float sharedCD = 0;
			float sharedABC = 0; float sharedACD = 0; float sharedBCD = 0; float sharedABD = 0;
			vector<double> data;
			//get sabund vector pointers so you can use the single calculators
			//one for each group
			SAbundVector sA, sB, sC, sD;
			SAbundVector* sabundA; SAbundVector* sabundB; SAbundVector* sabundC; SAbundVector* sabundD;
			sA = lookup[0]->getSAbundVector();  sabundA = &sA;
			sB = lookup[1]->getSAbundVector();  sabundB = &sB;
			sC = lookup[2]->getSAbundVector();  sabundC = &sC;
			sD = lookup[3]->getSAbundVector();  sabundD = &sD;
			
			//A = red, B = green, C = blue, D = yellow
			
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
				
				if ((vCalcs[i]->getName() != "sharedsobs") && (vCalcs[i]->getName() != "sharedchao")) { m->mothurOut(vCalcs[i]->getName() + " is not a valid calculator with four groups.  It will be disregarded. "); m->mothurOutEndLine(); }
				else{
					string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + lookup[0]->getGroup() + "-" + lookup[1]->getGroup() + "-" + lookup[2]->getGroup() + "-" + lookup[3]->getGroup() + ".svg";
					outputNames.push_back(filenamesvg);
					util.openOutputFile(filenamesvg, outsvg);

					if (m->getControl_pressed()) { outsvg.close(); return outputNames; }
					
					//in essence you want to run it like a single 
					if (vCalcs[i]->getName() == "sharedsobs") {
						singleCalc = new Sobs();
					}else if (vCalcs[i]->getName() == "sharedchao") {
						singleCalc = new Chao1();
					}
				
					//get estimates for numA
					data = singleCalc->getValues(sabundA);
					numA = data[0];
	//cout << "num a = " << numA << endl;	
 			
					//get estimates for numB
					data = singleCalc->getValues(sabundB);
					numB = data[0];
 	//cout << "num b = " << numB << endl;				
					//get estimates for numC
					data = singleCalc->getValues(sabundC);
					numC = data[0];
	//cout << "num c = " << numC << endl;				
					//get estimates for numD
					data = singleCalc->getValues(sabundD);
					numD = data[0];
//cout << "num d = " << numD << endl;	
                    
                    ofstream outShared;
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs")) {
                        string groupsTag = "";
                        for (int h = 0; h < lookup.size()-1; h++) { groupsTag += lookup[h]->getGroup() + "-"; }  groupsTag += lookup[lookup.size()-1]->getGroup();
                        string filenameShared = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + "." + vCalcs[i]->getName() + "." + groupsTag + ".sharedotus";
                        
                        outputNames.push_back(filenameShared);
                        
                        util.openOutputFile(filenameShared, outShared);
                        outShared << "Groups\tNumShared\tOTULabels\n";
                    }

					//get estimates for pairs
					subset.clear();
                    vector<string> labels;
					subset.push_back(lookup[0]); subset.push_back(lookup[1]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedAB = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
	//cout << "num ab = " << sharedAB << endl;			
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedAC = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[2]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
	//cout << "num ac = " << sharedAC << endl;				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedAD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[3]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
	//cout << "num ad = " << sharedAD << endl;			
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedBC = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[1]->getGroup() + "-" + lookup[2]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
	//cout << "num bc = " << sharedBC << endl;				
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedBD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[1]->getGroup() + "-" + lookup[3]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
		//cout << "num bd = " << sharedBD << endl;						
					subset.clear();
					subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedCD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[2]->getGroup() + "-" + lookup[3]->getGroup() << '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
						
	//cout << "num cd = " << sharedCD << endl;				
					//get estimates for combos of 3
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedABC = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup()+ "-" + lookup[2]->getGroup()<< '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
		//cout << "num abc = " << sharedABC << endl;					
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedACD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[2]->getGroup()+ "-" + lookup[3]->getGroup()<< '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                    }
			//cout << "num acd = " << sharedACD << endl;	
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset,labels);
					sharedBCD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[1]->getGroup() + "-" + lookup[2]->getGroup()+ "-" + lookup[3]->getGroup()<< '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        outShared << labels[labels.size()-1]; 
                        outShared << endl;
                    }
		//cout << "num bcd = " << sharedBCD << endl;		
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset, labels);
					sharedABD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup()+ "-" + lookup[3]->getGroup()<< '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        outShared << labels[labels.size()-1]; 
                        outShared << endl;
                    }
//cout << "num abd = " << sharedABD << endl;
					//get estimate for all four
					data = vCalcs[i]->getValues(lookup, labels);
					sharedABCD = data[0];
                    if (sharedOtus && (vCalcs[i]->getName() == "sharedsobs") &&  (labels.size() != 0)) {
                        outShared << lookup[0]->getGroup() + "-" + lookup[1]->getGroup() + "-" + lookup[2]->getGroup()+ "-" + lookup[3]->getGroup()<< '\t'<< labels.size() << '\t';
                        for (int k = 0; k < labels.size()-1; k++) {
                            outShared << labels[k] << ",";
                        }
                        if (labels.size() != 0) { outShared << labels[labels.size()-1]; }
                        outShared << endl;
                        outShared.close();
                    }
		//cout << "num abcd = " << sharedABCD << endl << endl;	
					int sharedVal, sharedABCVal, sharedABDVal, sharedACDVal, sharedBCDVal, sharedABVal, sharedACVal, sharedADVal, sharedBCVal, sharedBDVal, sharedCDVal, numSeqsA, numSeqsB, numSeqsC, numSeqsD;
												
					if (nseqs) {
						NSeqs* nseqsCalc = new NSeqs();
						vector<double> sharedData = nseqsCalc->getValues(lookup);
						
						vector<SharedRAbundVector*> mysubset; mysubset.push_back(lookup[0]); mysubset.push_back(lookup[1]);
						vector<double> sharedAB = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[0]); mysubset.push_back(lookup[2]);
						vector<double> sharedAC = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[0]); mysubset.push_back(lookup[3]);
						vector<double> sharedAD = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[1]); mysubset.push_back(lookup[2]);
						vector<double> sharedBC = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[1]); mysubset.push_back(lookup[3]);
						vector<double> sharedBD = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[2]); mysubset.push_back(lookup[3]);
						vector<double> sharedCD = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[0]);  mysubset.push_back(lookup[1]); mysubset.push_back(lookup[2]);
						vector<double> sharedABC = nseqsCalc->getValues(mysubset);
						
						mysubset.clear(); mysubset.push_back(lookup[0]);  mysubset.push_back(lookup[1]); mysubset.push_back(lookup[3]);
						vector<double> sharedABD = nseqsCalc->getValues(mysubset);

						mysubset.clear(); mysubset.push_back(lookup[0]);  mysubset.push_back(lookup[2]); mysubset.push_back(lookup[3]);
						vector<double> sharedACD = nseqsCalc->getValues(mysubset);

						mysubset.clear(); mysubset.push_back(lookup[1]);  mysubset.push_back(lookup[2]); mysubset.push_back(lookup[3]);
						vector<double> sharedBCD = nseqsCalc->getValues(mysubset);
						
						sharedVal = sharedData[0] + sharedData[1] + sharedData[2] + sharedData[3];
						sharedABCVal = sharedABC[0] + sharedABC[1] + sharedABC[2];
						sharedABDVal = sharedABD[0] + sharedABD[1] + sharedABD[2];
						sharedACDVal = sharedACD[0] + sharedACD[1] + sharedACD[2];
						sharedBCDVal = sharedBCD[0] + sharedBCD[1] + sharedBCD[2];
						sharedABVal = sharedAB[0] + sharedAB[1];
						sharedACVal = sharedAC[0] + sharedAC[1];
						sharedADVal = sharedAD[0] + sharedAD[1];
						sharedBCVal = sharedBC[0] + sharedBC[1];
						sharedBDVal = sharedBD[0] + sharedBD[1];
						sharedCDVal = sharedCD[0] + sharedCD[1];
						numSeqsA = sabundA->getNumSeqs(); 
						numSeqsB = sabundB->getNumSeqs(); 
						numSeqsC = sabundC->getNumSeqs(); 
						numSeqsD = sabundD->getNumSeqs(); 
						
						delete nseqsCalc;
					}

						
					//image window
					outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " +  toString(width) + " " + toString(windowSize) + " \" >\n";
					outsvg << "<g>\n";
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" +  toString(width) +  "\" height=\"" +  toString(windowSize) +  "\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.40 * width)) +  "\" y=\"" +  toString(int(0.05 * height)) +  "\" >Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";

					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\"  y=\"" +  toString(int(0.625 * height)) +  "\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA); 
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsA);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.65 * height)) +  "\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsB);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.675 * height)) +  "\">The number of species in group " + lookup[2]->getGroup() + " is " + toString(numC);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsC);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.7 * height)) +  "\">The number of species in group " + lookup[3]->getGroup() + " is " + toString(numD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(numSeqsD);  }  
					outsvg << "</text>\n";
					
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.725 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(sharedAB);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedABVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.75 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedAC);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedACVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.775 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedAD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedADVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.8 * height)) +  "\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedBC);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedBCVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.825 * height)) +  "\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedBD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedBDVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.85 * height)) +  "\">The number of species shared between groups " + lookup[2]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedCD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedCDVal);  }  
					outsvg << "</text>\n";
					
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.875 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + ", " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedABC);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedABCVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.9 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + ", " + lookup[1]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedABD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedABDVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.925 * height)) +  "\">The number of species shared between groups " + lookup[0]->getGroup() + ", " + lookup[2]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedACD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedACDVal);  }  
					outsvg << "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\"  font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.95 * height)) +  "\">The number of species shared between groups " + lookup[1]->getGroup() + ", " + lookup[2]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedBCD);
					if (nseqs) {  outsvg << ", and the number of squences is " + toString(sharedBCDVal);  }  
					outsvg << "</text>\n";
									
					//make adjustments
					sharedABC = sharedABC - sharedABCD;
			//cout << "num abc = " << sharedABC << endl;		
					sharedABD = sharedABD - sharedABCD;
				//cout << "num abd = " << sharedABD << endl;
					sharedACD = sharedACD - sharedABCD;
				//cout << "num acd = " << sharedACD << endl;
					sharedBCD = sharedBCD - sharedABCD;
				//cout << "num bcd = " << sharedBCD << endl;
					
					sharedAB = sharedAB - sharedABC - sharedABCD - sharedABD;  //cout << "num ab = " << sharedAB << endl;
					sharedAC = sharedAC - sharedABC - sharedABCD - sharedACD;  //cout << "num ac = " << sharedAC << endl;
					sharedAD = sharedAD - sharedABD - sharedABCD - sharedACD;  //cout << "num ad = " << sharedAD << endl;				
					sharedBC = sharedBC - sharedABC - sharedABCD - sharedBCD;  //cout << "num bc = " << sharedBC << endl;
					sharedBD = sharedBD - sharedABD - sharedABCD - sharedBCD; // cout << "num bd = " << sharedBD << endl; 
					sharedCD = sharedCD - sharedACD - sharedABCD - sharedBCD;  //cout << "num cd = " << sharedCD << endl;
					
					numA = numA - sharedAB - sharedAC - sharedAD - sharedABCD - sharedABC - sharedACD - sharedABD;
			//cout << "num a = " << numA << endl;		
					numB = numB - sharedAB - sharedBC - sharedBD - sharedABCD - sharedABC - sharedABD - sharedBCD;
			//cout << "num b = " << numB << endl;		
					numC = numC - sharedAC - sharedBC - sharedCD - sharedABCD - sharedABC - sharedACD - sharedBCD;
			//cout << "num c = " << numC << endl;		
					numD = numD - sharedAD - sharedBD - sharedCD - sharedABCD - sharedBCD - sharedACD - sharedABD;
			//cout << "num d = " << numD << endl;		
					
					//draw circles
					outsvg << "<ellipse fill=\"red\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-45 " +  toString(int(0.51 * width)) +  " " +  toString(int(0.27 * height)) +  ") \" cx=\"" +  toString(int(0.51 * width)) +  "\" cy=\"" +  toString(int(0.27 * height)) +  "\" rx=\"" +  toString(int(0.29 * width)) +  "\" ry=\"" +  toString(int(0.14 * height)) +  "\"/>\n "; 
					outsvg << "<ellipse fill=\"green\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+45 " +  toString(int(0.51 * width)) +  " " +  toString(int(0.27 * height)) +  ") \" cx=\"" +  toString(int(0.51 * width)) +  "\" cy=\"" +  toString(int(0.27 * height)) +  "\" rx=\"" +  toString(int(0.29 * width)) +  "\" ry=\"" +  toString(int(0.14 * height)) +  "\"/>\n ";
					outsvg << "<ellipse fill=\"blue\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-40 " +  toString(int(0.63 * width)) +  " " +  toString(int(0.36 * height)) +  ") \" cx=\"" +  toString(int(0.63 * width)) +  "\" cy=\"" +  toString(int(0.36 * height)) +  "\" rx=\"" +  toString(int(0.29 * width)) +  "\" ry=\"" +  toString(int(0.14 * height)) +  "\"/>\n ";
					outsvg << "<ellipse fill=\"yellow\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+40 " +  toString(int(0.36 * width)) +  " " +  toString(int(0.36 * height)) +  ") \" cx=\"" +  toString(int(0.36 * width)) +  "\" cy=\"" +  toString(int(0.36 * height)) +  "\" rx=\"" +  toString(int(0.29 * width)) +  "\" ry=\"" +  toString(int(0.14 * height)) +  "\"/>\n ";
			
					//A = red, B = green, C = blue, D = yellow
			
					//place labels within overlaps
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.66 * width) - ((int)toString(numA).length() / 2)) + "\" y=\"" +  toString(int(0.14 * height)) +  "\">" + toString(numA)  + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.66 * width) - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"" +  toString(int(0.11 * height)) +  "\">" + lookup[0]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedAB).length() / 2)) + "\"  y=\"" +  toString(int(0.2 * height)) +  "\">" + toString(sharedAB) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.36 * width) - ((int)toString(numB).length() / 2)) + "\"  y=\"" +  toString(int(0.14 * height)) +  "\">" + toString(numB)  + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.36 * width) - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"" +  toString(int(0.11 * height)) +  "\">" + lookup[1]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.7 * width) - ((int)toString(sharedAC).length() / 2)) + "\"  y=\"" +  toString(int(0.24 * height)) +  "\">" + toString(sharedAC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.785 * width) - ((int)toString(numC).length() / 2)) + "\"   y=\"" +  toString(int(0.29 * height)) +  "\">" + toString(numC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.785 * width) - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"" +  toString(int(0.26 * height)) +  "\">" + lookup[2]->getGroup() + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.31 * width) - ((int)toString(sharedBD).length() / 2)) + "\" y=\"" +  toString(int(0.24 * height)) +  "\">" + toString(sharedBD) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.22 * width) - ((int)toString(numD).length() / 2)) + "\"   y=\"" +  toString(int(0.29 * height)) +  "\">" + toString(numD) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\"  x=\"" + toString(int(0.22 * width) - ((int)lookup[3]->getGroup().length() / 2)) + "\"   y=\"" +  toString(int(0.26 * height)) +  "\">" + lookup[3]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.34 * width) - ((int)toString(sharedAD).length() / 2)) + "\" y=\"" +  toString(int(0.41 * height)) +  "\">" + toString(sharedAD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.675 * width) - ((int)toString(sharedBC).length() / 2)) + "\" y=\"" +  toString(int(0.41 * height)) +  "\">" + toString(sharedBC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedCD).length() / 2)) + "\" y=\"" +  toString(int(0.54 * height)) +  "\">" + toString(sharedCD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.39 * width) - ((int)toString(sharedABD).length() / 2)) + "\" y=\"" +  toString(int(0.3 * height)) +  "\">" + toString(sharedABD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.57 * width) - ((int)toString(sharedBCD).length() / 2)) + "\" y=\"" +  toString(int(0.45 * height)) +  "\">" + toString(sharedBCD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.435 * width) - ((int)toString(sharedACD).length() / 2)) + "\" y=\"" +  toString(int(0.45 * height)) +  "\">" + toString(sharedACD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.63 * width) - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"" +  toString(int(0.3 * height)) +  "\">" + toString(sharedABC) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(int(0.5 * width) - ((int)toString(sharedABCD).length() / 2)) + "\"  y=\"" +  toString(int(0.4 * height)) +  "\">" + toString(sharedABCD) + "</text>\n"; 
					
					outsvg << "<text fill=\"black\"  class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" +  toString(int(0.25 * width)) +  "\" y=\"" +  toString(int(0.975 * height)) +  "\">The total richness of all the groups is " + toString((float)(numA + numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD));
					if (nseqs) {  outsvg << ", and the number of squences in the otus shared by all groups is " + toString(sharedVal);  }  
					outsvg << "</text>\n";
					
					outsvg << "</g>\n</svg>\n";
					outsvg.close();
					delete singleCalc;
				}
			}
		}
		
		return outputNames;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Venn", "getPic");
		exit(1);
	}
}


