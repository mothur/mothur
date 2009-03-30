/*
 *  heatmap.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/25/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmap.h"

//**********************************************************************************************************************
HeatMap::HeatMap(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		sorted = globaldata->getSorted();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function HeatMap. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function HeatMap. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void HeatMap::getPic(OrderVector* order) {
	try {
		colorScale.clear();
		
		rabund = order->getRAbundVector();
		
		for (int i = 0; i < rabund.size(); i++) {
			colorScale[rabund.get(i)] = "";
		}
		
		float scaler = 255 / (float) colorScale.size();
		
		//go through map and give each score a color value
		for (it = colorScale.begin(); it != colorScale.end(); it++) {
			it->second = toHex(int(float(it->first) * scaler));
			if(it->second.length() == 1) {  it->second = "0" + it->second;  }
		}

		string filenamesvg = globaldata->inputFileName + ".heatmap." + order->getLabel() + ".svg";
		
		openOutputFile(filenamesvg, outsvg);
		
		//scale max rank so the maxrank = bright red
			
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 300 " + toString((rabund.getNumBins()*5 + 15))  + "\">\n";
		outsvg << "<g>\n";
		
		int x = 15;
		int y = 15;
		string color;

		for (int i = 0; i <= rabund.getNumBins(); i++) {
		
			color = colorScale[rabund.get(i)];
			
			outsvg << "<rect fill=\"#" + color + "0000\" stroke=\"#" + color + "0000\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"300\" height=\"5\"/>\n";
			y += 5;
		}
		outsvg << "</g>\n</svg>\n";
		
		outsvg.close();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void HeatMap::getPic(SharedOrderVector* sharedorder) {
	try {
		colorScale.clear();
		
		//fills vector of sharedsabunds - lookup
		getSharedVectors(sharedorder);
		
		//sort lookup so shared bins are on top
		if (sorted == "1") {  sortSharedVectors();  }
		
		//get maxBin
		for (int i = 0; i < lookup.size(); i++) {
			for (int j = 0; j < lookup[i]->size(); j++) {
				colorScale[lookup[i]->getAbundance(j)] = "";
			}
		}
		
		//get scaler
		float scaler = 255 / (float) colorScale.size();
		
		//go through map and give each score a color value
		for (it = colorScale.begin(); it != colorScale.end(); it++) {
			it->second = toHex(int(float(it->first) * scaler));
			if(it->second.length() == 1) {  it->second = "0" + it->second;  }
		}
		
		string filenamesvg = globaldata->inputFileName + ".heatmap." + sharedorder->getLabel() + "." + groupComb + ".svg";
		openOutputFile(filenamesvg, outsvg);
		
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(lookup.size() * 300) + " " + toString((lookup[0]->getNumBins()*5 + 15))  + "\">\n";
		outsvg << "<g>\n";
		
		int x = 15;
		int y = 15;
		string color;

		for (int i = 0; i <= lookup[0]->getNumBins(); i++) {
		
			for (int j = 0; j < lookup.size(); j++) {
			
				color = colorScale[lookup[j]->getAbundance(i)];	
						
				outsvg << "<rect fill=\"#" + color + "0000\" stroke=\"#" + color + "0000\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"300\" height=\"5\"/>\n";
				x += 300;
			}
			x = 15;
			y += 5;
		}
		outsvg << "</g>\n</svg>\n";
		
		outsvg.close();

		
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void HeatMap::getSharedVectors(SharedOrderVector* order){
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
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}


//**********************************************************************************************************************
void HeatMap::sortSharedVectors(){
	try {
		//copy lookup and then clear it to refill with sorted.
		//loop though lookup and determine if they are shared
		//if they are then insert in the front
		//if not push to back
		
		bool shared;
		vector<SharedRAbundVector*> looktemp;
		
		//create and initialize looktemp as a copy of lookup
		for (int i = 0; i < lookup.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(lookup[i]->getNumBins());
			temp->setLabel(lookup[i]->getLabel());
			temp->setGroup(lookup[i]->getGroup());
			//copy lookup i's info
			for (int j = 0; j < lookup[i]->size(); j++) {
				temp->set(j, lookup[i]->getAbundance(j), lookup[i]->getGroup());
			}
			looktemp.push_back(temp);
		}
		
		//clear out lookup to create sorted lookup
		lookup.clear();
		
		//create and initialize lookup to empty vectors
		for (int i = 0; i < looktemp.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector();
			lookup.push_back(temp);
		}
		
		//for each bin
		for (int i = 0; i < looktemp[0]->size(); i++) {
			shared = true;
			//for each group
			for (int j = 0; j < looktemp.size(); j++) {
				if (looktemp[j]->getAbundance(i) == 0) { shared = false; }
			}
			
			//fill lookup
			for (int j = 0; j < looktemp.size(); j++) {
				//if they are not shared then push to back, if they are not insert in front
				if (shared == false)  { lookup[j]->push_back(looktemp[j]->getAbundance(i), i, looktemp[j]->getGroup()); }
				else { lookup[j]->push_front(looktemp[j]->getAbundance(i), i, looktemp[j]->getGroup()); }
			}
		}
		
		//delete looktemp
		for (int j = 0; j < looktemp.size(); j++) {
			delete looktemp[j];
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function sortSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function sortSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

//**********************************************************************************************************************





