/*
 *  venn.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "venn.h"

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
void Venn::getPic(OrderVector* order) {
	try {
		
		rabund = order->getRAbundVector();
		
		string filenamesvg = globaldata->inputFileName + ".venn." + order->getLabel() + ".svg";
		
		openOutputFile(filenamesvg, outsvg);
		
			
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 200 200\" >\n";
		outsvg << "<g>\n";
				
		outsvg << "<circle fill=\"red\" stroke=\"black\" cx=\"150\" cy=\"150\" r=\"100\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"140\" y=\"150\">" + toString(rabund.getNumBins()) + "</text>\n";  
		outsvg << "</g>\n</svg>\n";
		
		outsvg.close();
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
void Venn::getPic(SharedOrderVector* sharedorder) {
	try {
		
		//fills vector of sharedsabunds - lookup
		getSharedVectors(sharedorder);
				
		string filenamesvg = globaldata->inputFileName + ".venn." + sharedorder->getLabel() + "." + groupComb + ".svg";
		openOutputFile(filenamesvg, outsvg);
		
		//image window
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 700 700\" >\n";
		outsvg << "<g>\n";
				
		if (lookup.size() == 1) {
			outsvg << "<circle fill=\"red\" opacity=\".5\" stroke=\"black\" cx=\"150\" cy=\"150\" r=\"100\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"140\" y=\"150\">" + toString(lookup[0]->getNumBins()) + "</text>\n";  
			outsvg << "</g>\n</svg>\n";
			
		}else if (lookup.size() == 2) {
			//calc the shared otu
			int shared = 0;
			int numA = 0;
			int numB = 0;
			
			float rScaler;
			
			//for each bin
			for (int i = 0; i < lookup[0]->size(); i++) {
				if (lookup[0]->getAbundance(i) != 0) { numA++; }
				if (lookup[1]->getAbundance(i) != 0) { numB++; }
				//are they shared
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0)) { shared++; }
			}
			
			if (numA > numB) { rScaler = 150 / float(numA); }
			else { rScaler = 150 / float(numB); }
			
			//to determine how far over to overlap b onto a.
			float percentOverlap = shared / (float) numA;	
			
			int bx = 200 + (numA * rScaler) - ((2 * (numA * rScaler)) * percentOverlap) + (numB * rScaler);
			int leftedgeB = bx - (numB * rScaler);  //center b - b's radius
			int leftedgeA = 200 - (numA * rScaler); //center a - a's radius
			int rightedgeB = bx + (numB * rScaler);  //center b + b's radius
			int rightedgeA = 200 + (numA * rScaler); //center a + a's radius

			int mida = leftedgeA + ((leftedgeB - leftedgeA) / 2);
			int midb = rightedgeA + ((rightedgeB - rightedgeA) / 2);
			int midab = leftedgeB + ((rightedgeA - leftedgeB) / 2);
			
			//draw circles
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(leftedgeA + ((rightedgeB - leftedgeA) / 2) - 70) + "\" y=\"50\">Venn Diagram at distance " + sharedorder->getLabel() + "</text>\n"; 
			outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".5\" stroke=\"black\" cx=\"200\" cy=\"250\" r=\"" + toString(numA * rScaler) + "\"/>"; 
			outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".5\" stroke=\"black\" cx=\"" + toString(bx)  + "\" cy=\"250\" r=\"" + toString(numB * rScaler) + "\"/>";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(mida) + "\" y=\"250\">" + toString(numA-shared) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(midab) + "\" y=\"250\">" + toString(shared) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(midb) + "\" y=\"250\">" + toString(numB-shared) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(mida) + "\" y=\"500\">Percentage of species that are shared is " + toString((shared / (float)(numA + numB - shared))) + "</text>\n"; 
			
		}else if (lookup.size() == 3) {
			//calc the shared otu
			int sharedABC = 0;
			int numA = 0; int numB = 0; int numC = 0;
			int sharedAB = 0; int sharedAC = 0; int sharedBC = 0;
			
			//float scalerB;
			
			//for each bin
			for (int i = 0; i < lookup[0]->size(); i++) {
				if (lookup[0]->getAbundance(i) != 0) { numA++; }
				if (lookup[1]->getAbundance(i) != 0) { numB++; }
				if (lookup[2]->getAbundance(i) != 0) { numC++; }
				//are they shared by 2
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0)) { sharedAB++; }
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0)) { sharedAC++; }
				if ((lookup[2]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0)) { sharedBC++; }
				
				//are they shared by all
				if ((lookup[0]->getAbundance(i) != 0) && (lookup[1]->getAbundance(i) != 0) && (lookup[2]->getAbundance(i) != 0)) { sharedABC++; }
			}
						
			//draw circles
			outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"700\" height=\"700\"/>"; 
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"265\" y=\"30\">Venn Diagram at distance " + sharedorder->getLabel() + "</text>\n";
			outsvg << "<circle fill=\"rgb(255,0,0)\" opacity=\".3\" stroke=\"black\" cx=\"230\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<circle fill=\"rgb(0,255,0)\" opacity=\".3\" stroke=\"black\" cx=\"455\" cy=\"200\" r=\"150\"/>"; 
			outsvg << "<circle fill=\"rgb(0,0,255)\" opacity=\".3\" stroke=\"black\" cx=\"343\" cy=\"400\" r=\"150\"/>"; 
			
			//place labels within overlaps
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(200 - ((int)toString(numA-sharedAB-sharedAC+sharedABC).length() / 2)) + "\" y=\"170\">" + toString(numA-sharedAB-sharedAC+sharedABC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedAB).length() / 2)) + "\"  y=\"170\">" + toString(sharedAB) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(485 - ((int)toString(numB-sharedAB-sharedBC+sharedABC).length() / 2)) + "\"  y=\"170\">" + toString(numB-sharedAB-sharedBC+sharedABC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(268 - ((int)toString(sharedAC).length() / 2)) + "\"  y=\"305\">" + toString(sharedAC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\"  x=\"" + toString(343 - ((int)toString(numC-sharedAC-sharedBC+sharedABC).length() / 2)) + "\"   y=\"430\">" + toString(numC-sharedAC-sharedBC+sharedABC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(408 - ((int)toString(sharedBC).length() / 2)) + "\" y=\"305\">" + toString(sharedBC) + "</text>\n";  
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(343 - ((int)toString(sharedABC).length() / 2)) + "\"  y=\"280\">" + toString(sharedABC) + "</text>\n"; 
			
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"580\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[1] + " is " + toString((sharedAB / (float)(numA + numB - sharedAB))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"610\">Percentage of species that are shared in groups " + globaldata->Groups[0] + " and " + globaldata->Groups[2] + " is " + toString((sharedAC / (float)(numA + numC - sharedAC))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"640\">Percentage of species that are shared in groups " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString((sharedBC / (float)(numB + numC - sharedBC))) + "</text>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"100\" y=\"670\">Percentage of species that are shared in groups " + globaldata->Groups[0] + ", " + globaldata->Groups[1] + " and " + globaldata->Groups[2] + " is " + toString((sharedABC / (float)(numA + numB + numC - sharedAB - sharedAC - sharedBC - (2 * sharedABC)))) + "</text>\n";
		}
		
		outsvg << "</g>\n</svg>\n";
		outsvg.close();

		
		
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
void Venn::getSharedVectors(SharedOrderVector* order){
	try {
	
		//delete lookup
		for (int j = 0; j < lookup.size(); j++) {
			delete lookup[j];
		}

		lookup.clear();
		
		groupComb = "";
		
		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < globaldata->Groups.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(order->getNumBins());
			temp->setLabel(order->getLabel());
			temp->setGroup(globaldata->Groups[i]);
			groupComb += globaldata->Groups[i];
			lookup.push_back(temp);
		}
		
		int numSeqs = order->size();
		//sample all the members
		for(int i=0;i<numSeqs;i++){
			//get first sample
			individual chosen = order->get(i);
			int abundance; 
					
			//set info for sharedvector in chosens group
			for (int j = 0; j < lookup.size(); j++) { 
				if (chosen.group == lookup[j]->getGroup()) {
					 abundance = lookup[j]->getAbundance(chosen.bin);
					 lookup[j]->set(chosen.bin, (abundance + 1), chosen.group);
					 break;
				}
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Venn class Function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Venn class function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************



