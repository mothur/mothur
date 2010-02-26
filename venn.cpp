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
Venn::Venn(string o) : outputDir(o) {
	try {
		globaldata = GlobalData::getInstance();
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
			string filenamesvg = outputDir + getSimpleName(globaldata->inputFileName) + ".venn." + sabund->getLabel() + vCalcs[i]->getName() + ".svg";
			outputNames.push_back(filenamesvg);
			openOutputFile(filenamesvg, outsvg);

			vector<double> data = vCalcs[i]->getValues(sabund);
			
			//svg image
			outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
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
		
		/******************* 1 Group **************************/
		if (lookup.size() == 1) {
					
			SAbundVector s;
			s = lookup[0]->getSAbundVector();  SAbundVector* sabund = &s;
			
			//make a file for each calculator
			for(int i=0;i<vCalcs.size();i++){
				string filenamesvg = outputDir + getSimpleName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
				outputNames.push_back(filenamesvg);
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
				outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
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
				string filenamesvg = outputDir + getSimpleName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
				outputNames.push_back(filenamesvg);
				openOutputFile(filenamesvg, outsvg);
				
				//get estimates for sharedAB
				vector<double> shared = vCalcs[i]->getValues(subset);
				
				//in essence you want to run it like a single 
				if (vCalcs[i]->getName() == "sharedsobs") {
					singleCalc = new Sobs();
				}else if (vCalcs[i]->getName() == "sharedchao") {
					singleCalc = new Chao1();
				}//else if (vCalcs[i]->getName() == "sharedace") {
					//singleCalc = new Ace(10);
				//}
				
				//get estimates for numA
				vector<double> numA = singleCalc->getValues(sabundA);

				//get estimates for numB
				vector<double> numB = singleCalc->getValues(sabundB);
						
				//image window
				outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
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
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"460\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA[0]);
				if (numA.size() == 3) { 
					outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }
		
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"480\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB[0]);
				if (numB.size() == 3) { 
					outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]) + "</text>\n";
				}else { outsvg << "</text>\n"; }

				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"500\">The number of sepecies shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(shared[0]) + "</text>\n";
				outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"520\">Percentage of species that are shared in groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString((shared[0] / (float)(numA[0] + numB[0] - shared[0]))*100) + "</text>\n";
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
			
				string filenamesvg = outputDir + getSimpleName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
				outputNames.push_back(filenamesvg);
				openOutputFile(filenamesvg, outsvg);
				
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
						for (int j = 0; j < lookup[1]->size(); j++) {
							merge->push_back((lookup[1]->getAbundance(j) + lookup[2]->getAbundance(j)), "");
						}
					
						subset.clear();
						subset.push_back(lookup[0]); subset.push_back(merge);
						sharedAwithBC = vCalcs[i]->getValues(subset);
				
						delete merge;
						//merge AC and estimate with shared with B
						merge = new SharedRAbundVector();
						for (int j = 0; j < lookup[0]->size(); j++) {
							merge->push_back((lookup[0]->getAbundance(j) + lookup[2]->getAbundance(j)), "");
						}
					
						subset.clear();
						subset.push_back(merge); subset.push_back(lookup[1]);
						sharedBwithAC = vCalcs[i]->getValues(subset);
				
						delete merge;
						//merge AB and estimate with shared with C
						merge = new SharedRAbundVector();
						for (int j = 0; j < lookup[0]->size(); j++) {
							merge->push_back((lookup[0]->getAbundance(j) + lookup[1]->getAbundance(j)), "");
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
					outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 800 800\" >\n";
					outsvg << "<g>\n";

					//draw circles
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"800\" height=\"800\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
					outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"230\" cy=\"200\" r=\"150\"/>"; 
					outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"455\" cy=\"200\" r=\"150\"/>"; 
					outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"343\" cy=\"400\" r=\"150\"/>"; 

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
				
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"660\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(sharedAB[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"680\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedAC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"700\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedBC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"720\">The number of species shared between groups " + lookup[0]->getGroup() + " and combined groups " + lookup[1]->getGroup() + lookup[2]->getGroup() + " is " + toString(sharedAwithBC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"740\">The number of species shared between groups " + lookup[1]->getGroup() + " and combined groups " + lookup[0]->getGroup() + lookup[2]->getGroup() + " is " + toString(sharedBwithAC[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"760\">The number of species shared between groups " + lookup[2]->getGroup() + " and combined groups " + lookup[0]->getGroup() + lookup[1]->getGroup() + " is " + toString(sharedCwithAB[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"580\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA[0]);
					if (numA.size() == 3) { 
						outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]) + "</text>\n";
					}else { outsvg << "</text>\n"; }
			
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"600\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB[0]);
					if (numB.size() == 3) { 
						outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]) + "</text>\n";
					}else { outsvg << "</text>\n"; }
					
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"620\">The number of species in group " + lookup[2]->getGroup() + " is " + toString(numC[0]);
					if (numC.size() == 3) { 
						outsvg << " the lci is " + toString(numC[1]) + " and the hci is " + toString(numC[2]) + "</text>\n";
					}else { outsvg << "</text>\n"; }

					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"640\">The total richness of all the groups is " + toString(numA[0] + numB[0] + numC[0] - sharedAB[0] - sharedAC[0] - sharedBC[0] + sharedABC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"780\">The total shared richness is " + toString(sharedABC) + "</text>\n";
					
					delete singleCalc;
					
				}else { //sharedchao and sharedsobs are multigroup
					
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
					vector<double> sharedab =  vCalcs[i]->getValues(subset);
					
					subset.clear(); 
					subset.push_back(lookup[0]); subset.push_back(lookup[2]);
					vector<double> sharedac =  vCalcs[i]->getValues(subset);
					
					subset.clear(); 
					subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					vector<double> sharedbc =  vCalcs[i]->getValues(subset);
					
					subset.clear(); 
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					vector<double> sharedabc =  vCalcs[i]->getValues(subset);
					
					//image window
					outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 800 800\" >\n";
					outsvg << "<g>\n";

					//draw circles
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"800\" height=\"800\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";
					outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"230\" cy=\"200\" r=\"150\"/>"; 
					outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"455\" cy=\"200\" r=\"150\"/>"; 
					outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"343\" cy=\"400\" r=\"150\"/>"; 

					//place labels within overlaps
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)toString(numA[0]-sharedab[0]-sharedac[0]+sharedabc[0]).length() / 2)) + "\" y=\"170\">" + toString(numA[0]-sharedab[0]-sharedac[0]+sharedabc[0]) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"150\">" + lookup[0]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedab[0] - sharedabc[0]).length() / 2)) + "\"  y=\"170\">" + toString(sharedab[0] - sharedabc[0]) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(485 - ((int)toString(numB[0]-sharedab[0]-sharedbc[0]+sharedabc[0]).length() / 2)) + "\"  y=\"170\">" + toString(numB[0]-sharedab[0]-sharedbc[0]+sharedabc[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(485 - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"150\">" + lookup[1]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(268 - ((int)toString(sharedac[0] - sharedabc[0]).length() / 2)) + "\"  y=\"305\">" + toString(sharedac[0] - sharedabc[0]) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(343 - ((int)toString(numC[0]-sharedac[0]-sharedbc[0]+sharedabc[0]).length() / 2)) + "\"   y=\"430\">" + toString(numC[0]-sharedac[0]-sharedbc[0]+sharedabc[0]) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(343 - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"410\">" + lookup[2]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(408 - ((int)toString(sharedbc[0] - sharedabc[0]).length() / 2)) + "\" y=\"305\">" + toString(sharedbc[0] - sharedabc[0]) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedabc[0]).length() / 2)) + "\"  y=\"280\">" + toString(sharedabc[0]) + "</text>\n"; 
				
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"660\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(sharedab[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"680\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedac[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"700\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedbc[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"580\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA[0]);
					if (numA.size() == 3) { 
						outsvg << " the lci is " + toString(numA[1]) + " and the hci is " + toString(numA[2]) + "</text>\n";
					}else { outsvg << "</text>\n"; }
			
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"600\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB[0]);
					if (numB.size() == 3) { 
						outsvg << " the lci is " + toString(numB[1]) + " and the hci is " + toString(numB[2]) + "</text>\n";
					}else { outsvg << "</text>\n"; }
					
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"620\">The number of species in group " + lookup[2]->getGroup() + " is " + toString(numC[0]);
					if (numC.size() == 3) { 
						outsvg << " the lci is " + toString(numC[1]) + " and the hci is " + toString(numC[2]) + "</text>\n";
					}else { outsvg << "</text>\n"; }

					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"640\">The total richness of all the groups is " + toString(numA[0] + numB[0] + numC[0] - sharedab[0] - sharedac[0] - sharedbc[0] + sharedabc[0]) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"720\">The total shared richness is " + toString(sharedabc[0]) + "</text>\n";


				}
		
								
				//close file
				outsvg << "</g>\n</svg>\n";
				outsvg.close();
				

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
				
				if ((vCalcs[i]->getName() != "sharedsobs") && (vCalcs[i]->getName() != "sharedchao")) { m->mothurOut(vCalcs[i]->getName() + " is not a valid calculator with four groups.  It will be disregarded. "); m->mothurOutEndLine(); }
				else{
					string filenamesvg = outputDir + getSimpleName(globaldata->inputFileName) + lookup[0]->getLabel() + ".venn." + vCalcs[i]->getName() + ".svg";
					outputNames.push_back(filenamesvg);
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

					//get estimates for pairs
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]);
					data = vCalcs[i]->getValues(subset);
					sharedAB = data[0];
	//cout << "num ab = " << sharedAB << endl;			
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset);
					sharedAC = data[0];
	//cout << "num ac = " << sharedAC << endl;				
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedAD = data[0];
	//cout << "num ad = " << sharedAD << endl;			
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset);
					sharedBC = data[0];
	//cout << "num bc = " << sharedBC << endl;				
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedBD = data[0];
		//cout << "num bd = " << sharedBD << endl;						
					subset.clear();
					subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedCD = data[0];
						
	//cout << "num cd = " << sharedCD << endl;				
					//get estimates for combos of 3
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[2]);
					data = vCalcs[i]->getValues(subset);
					sharedABC = data[0];
		//cout << "num abc = " << sharedABC << endl;					
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedACD = data[0];
			//cout << "num acd = " << sharedACD << endl;	
					subset.clear();
					subset.push_back(lookup[1]); subset.push_back(lookup[2]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedBCD = data[0];
		//cout << "num bcd = " << sharedBCD << endl;		
					subset.clear();
					subset.push_back(lookup[0]); subset.push_back(lookup[1]); subset.push_back(lookup[3]);
					data = vCalcs[i]->getValues(subset);
					sharedABD = data[0];
//cout << "num abd = " << sharedABD << endl;
					//get estimate for all four
					data = vCalcs[i]->getValues(lookup);
					sharedABCD = data[0];
		//cout << "num abcd = " << sharedABCD << endl << endl;	
		
				
						
					//image window
					outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 700 800\" >\n";
					outsvg << "<g>\n";
					outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"800\"/>"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + lookup[0]->getLabel() + "</text>\n";

					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"490\">The number of species in group " + lookup[0]->getGroup() + " is " + toString(numA) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"510\">The number of species in group " + lookup[1]->getGroup() + " is " + toString(numB) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"530\">The number of species in group " + lookup[2]->getGroup() + " is " + toString(numC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"550\">The number of species in group " + lookup[3]->getGroup() + " is " + toString(numD) + "</text>\n";
					
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"570\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[1]->getGroup() + " is " + toString(sharedAB) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"590\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedAC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"610\">The number of species shared between groups " + lookup[0]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedAD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"630\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedBC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"650\">The number of species shared between groups " + lookup[1]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedBD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"670\">The number of species shared between groups " + lookup[2]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedCD) + "</text>\n";
					
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"690\">The number of species shared between groups " + lookup[0]->getGroup() + ", " + lookup[1]->getGroup() + " and " + lookup[2]->getGroup() + " is " + toString(sharedABC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"710\">The number of species shared between groups " + lookup[0]->getGroup() + ", " + lookup[1]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedABD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"730\">The number of species shared between groups " + lookup[0]->getGroup() + ", " + lookup[2]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedACD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"750\">The number of species shared between groups " + lookup[1]->getGroup() + ", " + lookup[2]->getGroup() + " and " + lookup[3]->getGroup() + " is " + toString(sharedBCD) + "</text>\n";
									
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
					outsvg << "<ellipse fill=\"red\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-45 355 215) \" cx=\"355\" cy=\"215\" rx=\"200\" ry=\"115\"/>\n "; 
					outsvg << "<ellipse fill=\"green\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+45 355 215) \" cx=\"355\" cy=\"215\" rx=\"200\" ry=\"115\"/>\n ";
					outsvg << "<ellipse fill=\"blue\" stroke=\"black\" opacity=\".35\" transform=\"rotate(-40 440 315) \" cx=\"440\" cy=\"315\" rx=\"200\" ry=\"115\"/>\n ";
					outsvg << "<ellipse fill=\"yellow\" stroke=\"black\" opacity=\".35\" transform=\"rotate(+40 270 315) \" cx=\"270\" cy=\"315\" rx=\"200\" ry=\"115\"/>\n ";
			
					//A = red, B = green, C = blue, D = yellow
			
					//place labels within overlaps
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(460 - ((int)toString(numA).length() / 2)) + "\" y=\"110\">" + toString(numA)  + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(460 - ((int)lookup[0]->getGroup().length() / 2)) + "\" y=\"90\">" + lookup[0]->getGroup() + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedAB).length() / 2)) + "\"  y=\"160\">" + toString(sharedAB) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(250 - ((int)toString(numB).length() / 2)) + "\"  y=\"110\">" + toString(numB)  + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(250 - ((int)lookup[1]->getGroup().length() / 2)) + "\"  y=\"90\">" + lookup[1]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(490 - ((int)toString(sharedAC).length() / 2)) + "\"  y=\"190\">" + toString(sharedAC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(550 - ((int)toString(numC).length() / 2)) + "\"   y=\"230\">" + toString(numC) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(550 - ((int)lookup[2]->getGroup().length() / 2)) + "\"   y=\"210\">" + lookup[2]->getGroup() + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(215 - ((int)toString(sharedBD).length() / 2)) + "\" y=\"190\">" + toString(sharedBD) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(150 - ((int)toString(numD).length() / 2)) + "\"   y=\"230\">" + toString(numD) + "</text>\n";  
					outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(150 - ((int)lookup[3]->getGroup().length() / 2)) + "\"   y=\"210\">" + lookup[3]->getGroup() + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(240 - ((int)toString(sharedAD).length() / 2)) + "\" y=\"325\">" + toString(sharedAD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(470 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"325\">" + toString(sharedBC) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedCD).length() / 2)) + "\" y=\"430\">" + toString(sharedCD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(275 - ((int)toString(sharedABD).length() / 2)) + "\" y=\"240\">" + toString(sharedABD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(400 - ((int)toString(sharedBCD).length() / 2)) + "\" y=\"360\">" + toString(sharedBCD) + "</text>\n";
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(305 - ((int)toString(sharedACD).length() / 2)) + "\" y=\"360\">" + toString(sharedACD) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(440 - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"240\">" + toString(sharedABC) + "</text>\n"; 
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(350 - ((int)toString(sharedABCD).length() / 2)) + "\"  y=\"320\">" + toString(sharedABCD) + "</text>\n"; 
					
					
										
					outsvg << "<text fill=\"black\" class=\"seri\" x=\"175\" y=\"770\">The total richness of all the groups is " + toString((float)(numA + numB + numC + numD + sharedAB + sharedAC + sharedAD + sharedBC + sharedBD + sharedCD + sharedABC + sharedABD + sharedACD + sharedBCD + sharedABCD)) + "</text>\n";
					

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


