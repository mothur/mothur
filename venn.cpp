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
#include "sharedchao1.h"
#include "sharedsobscollectsummary.h"


//**********************************************************************************************************************
Venn::Venn(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function Venn. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function Venn. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void Venn::getPic(SAbundVector* sabund, vector<Calculator*> vCalcs) {
	try {
				
		for(int i=0;i<vCalcs.size();i++){
			string filenamesvg = globaldata->inputFileName + ".venn." + sabund->getLabel() + vCalcs[i]->getName() + ".svg";
			openOutputFile(filenamesvg, outsvg);

			vector<double> data = vCalcs[i]->getValues(sabund);
			
			//svg image
			outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
			outsvg << "<g>\n";
				
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + sabund->getLabel() + "</text>\n";
			outsvg << "<circle fill=\"red\" opacity=\".5\" stroke=\"black\" cx=\"350\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(data[0]).length() / 2)) + "\" y=\"195\">" + toString(data[0]) + "</text>\n";  
			
			if (data.size() == 3) { 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"200\" y=\"380\">The lower bound of the confidence interval is " + toString(data[1]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"200\" y=\"410\">The upper bound of the confidence interval is " + toString(data[2]) + "</text>\n";
			}
			
			outsvg << "</g>\n</svg>\n";
			outsvg.close();
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void Venn::getPic(vector<SharedRAbundVector*> lookup, vector<Calculator*> vCalcs) {
	try {
		
		vector<SharedRAbundVector*> subset;
		
		/******************* 1 Group **************************/
		if (lookup.size() == 1) {
					
			SAbundVector s;
			s = lookup[0]->getSAbundVector();  SAbundVector* sabund = &s;
			
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
				string filenamesvg = getRootName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
				openOutputFile(filenamesvg, outsvg);
			
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
				outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
				outsvg << "<g>\n";
				
				outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
				outsvg << "<circle fill=\"red\" opacity=\".5\" stroke=\"black\" cx=\"350\" cy=\"200\" r=\"150\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"165\">" + lookup[0]->getGroup() + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(data[0]).length() / 2)) + "\" y=\"195\">" + toString(data[0]) + "</text>\n";  
			
				if (data.size() == 3) { 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"200\" y=\"380\">The lower bound of the confidence interval is " + toString(data[1]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"200\" y=\"410\">The upper bound of the confidence interval is " + toString(data[2]) + "</text>\n";
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
				string filenamesvg = getRootName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
				openOutputFile(filenamesvg, outsvg);
				
				//get estimates for sharedAB
				vector<double> shared = vCalcs[i]->getValues(subset);
				
				//in essence you want to run it like a single 
				if (vCalcs[i]->getName() == "sharedsobs") {
					singleCalc = new Sobs();
				}else if (vCalcs[i]->getName() == "sharedchao") {
					singleCalc = new Chao1();
				}else if (vCalcs[i]->getName() == "sharedace") {
					singleCalc = new Ace(10);
				}
				
				//get estimates for numA
				vector<double> numA = singleCalc->getValues(sabundA);

				//get estimates for numB
				vector<double> numB = singleCalc->getValues(sabundB);
						
				//image window
				outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
				outsvg << "<g>\n";

				//draw circles
				outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
				outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"250\" cy=\"200\" r=\"150\"/>"; 
				outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"435\" cy=\"200\" r=\"150\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)toString(numA[0]).length() / 2)) + "\" y=\"195\">" + toString(numA[0] - shared[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(490 - ((int)toString(numB[0]).length() / 2)) + "\" y=\"195\">" + toString(numB[0] - shared[0]) + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"175\">" + lookup[0]->getGroup() + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(490 - ((int)lookup[1]->getGroup().length() / 2)) + "\" y=\"175\">" + lookup[1]->getGroup() + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(shared[0]).length() / 2)) + "\" y=\"195\">" + toString(shared[0]) + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"460\">The number of species in group " + globaldata->Groups[0] + " is " + toString(numA[0]);
				if (numA.size() == 3) { 
					outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }
		
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"480\">The number of species in group " + globaldata->Groups[1] + " is " + toString(numB[0]);
				if (numB.size() == 3) { 
					outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }

				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"500\">The number of sepecies shared between groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString(shared[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"520\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString((shared[0] / (float)(numA[0] + numB[0] - shared[0]))) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"540\">The total richness for all groups is " + toString((float)(numA[0] + numB[0] - shared[0])) + "</text>\n";
				
				//close file
				outsvg << "</g>\n</svg>\n";
				outsvg.close();
				delete singleCalc;
			}
		/******************* 3 Groups **************************/
						
		}else if (lookup.size() == 3) {
			//get sabund vector pointers so you can use the single calculators
			//one for each group
			SAbundVector sA, sB, sC;
			SAbundVector* sabundA; SAbundVector* sabundB; SAbundVector* sabundC;
			sA = lookup[0]->getSAbundVector();  sabundA = &sA;
			sB = lookup[1]->getSAbundVector();  sabundB = &sB;
			sC = lookup[2]->getSAbundVector();  sabundC = &sC;
		
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
			
				//in essence you want to run it like a single 
				if (vCalcs[i]->getName() == "sharedsobs") {
					singleCalc = new Sobs();
				}else if (vCalcs[i]->getName() == "sharedchao") {
					singleCalc = new Chao1();
				}else if (vCalcs[i]->getName() == "sharedace") {
					singleCalc = new Ace(10);
				}
				
				//get estimates for numA
				vector<double> numA = singleCalc->getValues(sabundA);
 			
				//get estimates for numB
				vector<double> numB = singleCalc->getValues(sabundB);
 				
				//get estimates for numC
				vector<double> numC = singleCalc->getValues(sabundC);

				
				string filenamesvg = getRootName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
				openOutputFile(filenamesvg, outsvg);

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
				float sharedABC1, sharedABC2, sharedABC3, sharedABC;

				if (vCalcs[i]->getName() != "sharedchao") {
					//merge BC and estimate with shared with A
					SharedRAbundVector* merge = new SharedRAbundVector();
					for (int j = 0; j < lookup[1]->size(); j++) {
						merge->push_back((lookup[1]->getAbundance(j) + lookup[2]->getAbundance(j)), j, "");
					}
				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(merge);
					sharedAwithBC = vCalcs[i]->getValues(subset);
			
					delete merge;
					//merge AC and estimate with shared with B
					merge = new SharedRAbundVector();
					for (int j = 0; j < lookup[0]->size(); j++) {
						merge->push_back((lookup[0]->getAbundance(j) + lookup[2]->getAbundance(j)), j, "");
					}
				
					subset.clear();
					subset.push_back(merge); subset.push_back(lookup[1]);
					sharedBwithAC = vCalcs[i]->getValues(subset);
			
					delete merge;
					//merge AB and estimate with shared with C
					merge = new SharedRAbundVector();
					for (int j = 0; j < lookup[0]->size(); j++) {
						merge->push_back((lookup[0]->getAbundance(j) + lookup[1]->getAbundance(j)), j, "");
					}
				
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
				outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 800 800\" >\n";
				outsvg << "<g>\n";

				//draw circles
				outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"800\" height=\"800\"/>"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
				outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"230\" cy=\"200\" r=\"150\"/>"; 
				outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"455\" cy=\"200\" r=\"150\"/>"; 
				outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"343\" cy=\"400\" r=\"150\"/>"; 
//cout << "numA = " << numA[0] << " numB = " << numB[0] 	<< " numC = " << numC[0] << endl;
//cout << "sharedAB = " << sharedAB[0] << " sharedAC = " << sharedAC[0] << " sharedBC = " << sharedBC[0] << endl;
//cout << "sharedAwithBC = " << sharedAwithBC[0]	<< " sharedBwithAC = " << sharedBwithAC[0] << " sharedCwithAB = " << sharedCwithAB[0] << endl;	
				//place labels within overlaps
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)toString(numA[0]-sharedAwithBC[0]).length() / 2)) + "\" y=\"170\">" + toString(numA[0]-sharedAwithBC[0]) + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"150\">" + lookup[0]->getGroup() + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedAB[0] - sharedABC).length() / 2)) + "\"  y=\"170\">" + toString(sharedAB[0] - sharedABC) + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(485 - ((int)toString(numB[0]-sharedBwithAC[0]).length() / 2)) + "\"  y=\"170\">" + toString(numB[0]-sharedBwithAC[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(485 - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"150\">" + lookup[1]->getGroup() + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(268 - ((int)toString(sharedAC[0] - sharedABC).length() / 2)) + "\"  y=\"305\">" + toString(sharedAC[0] - sharedABC) + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(343 - ((int)toString(numC[0]-sharedCwithAB[0]).length() / 2)) + "\"   y=\"430\">" + toString(numC[0]-sharedCwithAB[0]) + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(343 - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"410\">" + lookup[2]->getGroup() + "</text>\n"; 
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(408 - ((int)toString(sharedBC[0] - sharedABC).length() / 2)) + "\" y=\"305\">" + toString(sharedBC[0] - sharedABC) + "</text>\n";  
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"280\">" + toString(sharedABC) + "</text>\n"; 
			
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"660\">The number of sepecies shared between groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString(sharedAB[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"680\">The number of sepecies shared between groups " + globaldata->Groups[0] + " and " + globaldata->Groups[2] + " is " + toString(sharedAC[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"700\">The number of sepecies shared between groups " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString(sharedBC[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"720\">The number of sepecies shared between groups " + globaldata->Groups[0] + " and combined groups " + globaldata->Groups[1] + globaldata->Groups[2] + " is " + toString(sharedAwithBC[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"740\">The number of sepecies shared between groups " + globaldata->Groups[1] + " and combined groups " + globaldata->Groups[0] + globaldata->Groups[2] + " is " + toString(sharedBwithAC[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"760\">The number of sepecies shared between groups " + globaldata->Groups[2] + " and combined groups " + globaldata->Groups[0] + globaldata->Groups[1] + " is " + toString(sharedCwithAB[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"580\">The number of species in group " + globaldata->Groups[0] + " is " + toString(numA[0]);
				if (numA.size() == 3) { 
					outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }
		
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"600\">The number of species in group " + globaldata->Groups[1] + " is " + toString(numB[0]);
				if (numB.size() == 3) { 
					outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }
				
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"620\">The number of species in group " + globaldata->Groups[2] + " is " + toString(numC[0]);
				if (numC.size() == 3) { 
					outsvg << " the lci is " + toString(numC[1]) + " and the hci is " + toString(numC[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }

				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"640\">The total richness of all the groups is " + toString(numA[0] + numB[0] + numC[0] - sharedAB[0] - sharedAC[0] - sharedBC[0] + sharedABC) + "</text>\n";
				
				//close file
				outsvg << "</g>\n</svg>\n";
				outsvg.close();
				delete singleCalc;

			}
			
		/******************* 4 Groups **************************/
		
		}else if (lookup.size() == 4) {
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
				
				if ((vCalcs[i]->getName() != "sharedsobs") && (vCalcs[i]->getName() != "sharedchao")) { cout << vCalcs[i]->getName() << " is not a valid calculator with four groups.  It will be disregarded. " << endl; }
				else{
					string filenamesvg = getRootName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
					openOutputFile(filenamesvg, outsvg);

				
					//in essence you want to run it like a single 
					if (vCalcs[i]->getName() == "sharedsobs") {
						singleCalc = new Sobs();
					}else if (vCalcs[i]->getName() == "sharedchao") {
						singleCalc = new Chao1();
					}
				
					//get estimates for numA
					data = singleCalc->getValues(sabundA);
					numA = data[0];
 			
					//get estimates for numB
					data = singleCalc->getValues(sabundB);
					numB = data[0];
 				
					//get estimates for numC
					data = singleCalc->getValues(sabundC);
					numC = data[0];
				
					//get estimates for numD
					data = singleCalc->getValues(sabundD);
					numD = data[0];


					//get estimates for pairs
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]);
					data = vCalcs[i]->getValues(subset);
					sharedAB = data[0];
				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset);
					sharedAC = data[0];
				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedAD = data[0];
				
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset);
					sharedBC = data[0];
				
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedBD = data[0];
				
					subset.clear();
					subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedCD = data[0];
				
				
					//get estimates for combos of 3
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset);
					sharedABC = data[0];
				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedACD = data[0];
				
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedBCD = data[0];
				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedABD = data[0];

					//get estimate for all four
					data = vCalcs[i]->getValues(lookup);
					sharedABCD = data[0];

		
					//image window
					outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
					outsvg << "<g>\n";

					//draw circles
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
					outsvg << "<ellipse fill=\"red\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-45 355 215) \" cx=\"355\" cy=\"215\" rx=\"200\" ry=\"115\"/>\n "; 
					outsvg << "<ellipse fill=\"green\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+45 355 215) \" cx=\"355\" cy=\"215\" rx=\"200\" ry=\"115\"/>\n ";
					outsvg << "<ellipse fill=\"blue\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-40 440 315) \" cx=\"440\" cy=\"315\" rx=\"200\" ry=\"115\"/>\n ";
					outsvg << "<ellipse fill=\"yellow\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+40 270 315) \" cx=\"270\" cy=\"315\" rx=\"200\" ry=\"115\"/>\n ";
			
					//A = red, B = green, C = blue, D = yellow
			
					//place labels within overlaps
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(460 - ((int)toString(numA).length() / 2)) + "\" y=\"110\">" + toString(numA) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(460 - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"90\">" + lookup[0]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedAB).length() / 2)) + "\"  y=\"160\">" + toString(sharedAB) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(250 - ((int)toString(numB).length() / 2)) + "\"  y=\"110\">" + toString(numB) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(250 - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"90\">" + lookup[1]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(490 - ((int)toString(sharedAC).length() / 2)) + "\"  y=\"190\">" + toString(sharedAC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(550 - ((int)toString(numC).length() / 2)) + "\"   y=\"230\">" + toString(numC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(550 - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"210\">" + lookup[2]->getGroup() + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(215 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"190\">" + toString(sharedBC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(150 - ((int)toString(numD).length() / 2)) + "\"   y=\"230\">" + toString(numD) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(150 - ((int)lookup[3]->getGroup().length() / 2)) + "\"   y=\"210\">" + lookup[3]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(240 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"325\">" + toString(sharedAD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(470 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"325\">" + toString(sharedBD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedCD).length() / 2)) + "\" y=\"430\">" + toString(sharedCD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(275 - ((int)toString(sharedABD).length() / 2)) + "\" y=\"240\">" + toString(sharedABD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(400 - ((int)toString(sharedBCD).length() / 2)) + "\" y=\"360\">" + toString(sharedBCD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(305 - ((int)toString(sharedACD).length() / 2)) + "\" y=\"360\">" + toString(sharedACD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(440 - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"240\">" + toString(sharedABC) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedABCD).length() / 2)) + "\"  y=\"320\">" + toString(sharedABCD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"250\" y=\"490\">The total richness of all the groups is " + toString((float)(numA + numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD)) + "</text>\n";
			
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"510\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[2] + " is " + toString(((sharedAC + sharedACD + sharedABC + sharedABCD) / (float)(numA + numC + sharedAB + sharedAC + sharedAD + sharedBC + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"530\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[3] + " is " + toString(((sharedAD + sharedACD + sharedABD + sharedABCD) / (float)(numA + numD + sharedAB + sharedAC + sharedAD + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"550\">Percentage of species that are shared in groups " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString(((sharedBC + sharedABC + sharedBCD + sharedABCD) / (float)(numB + numC + sharedAB + sharedAC + sharedCD + sharedBD + sharedBC + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"570\">Percentage of species that are shared in groups " + globaldata->Groups[1] + " and " + globaldata->Groups[3] + " is " + toString(((sharedBD + sharedABD + sharedBCD + sharedABCD) / (float)(numB + numD + sharedAB + sharedAD + sharedCD + sharedBD + sharedBC + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"590\">Percentage of species that are shared in groups " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString(((sharedCD + sharedBCD + sharedACD + sharedABCD) / (float)(numC + numD + sharedAC + sharedAD + sharedCD + sharedBD + sharedBC + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"610\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString(((sharedABC + sharedABCD) / (float)(numA + numB + numC + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"630\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + " and " + globaldata->Groups[3] + " is " + toString(((sharedABD + sharedABCD) / (float)(numA + numB + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"650\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString(((sharedACD + sharedABCD) / (float)(numA + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"670\">Percentage of species that are shared in groups " + globaldata->Groups[1] + ", " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString(((sharedBCD + sharedABCD) / (float)(numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
					//outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"690\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + ", " + globaldata->Groups[2] + " and " + globaldata->Groups[3] + " is " + toString((sharedABCD / (float)(numA + numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD))) + "</text>\n";
		
					outsvg << "</g>\n</svg>\n";
					outsvg.close();
					delete singleCalc;
				}
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}


