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



