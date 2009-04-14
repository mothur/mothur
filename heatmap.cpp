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
		util = new SharedUtil();
		
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

		string filenamesvg = getRootName(globaldata->inputFileName) + order->getLabel() + ".heatmap.svg";
		
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
		util->getSharedVectors(globaldata->Groups, lookup, sharedorder);  //fills group vectors from order vector.
		
		//sort lookup so shared bins are on top
		if (sorted == "1") {  sortSharedVectors();  }
		
		//get users scaling method
		scaler = globaldata->getScaler();
		
		int maxbin = 0;
		for (int i = 0; i < lookup.size(); i++) {
			for (int j = 0; j < lookup[i]->size(); j++) {
				//if (lookup[i]->getAbundance(j) != 0) { //don't want log value of 0.
					//if (scaler == "log10") {
					//	colorScale[-log((log10(lookup[i]->getAbundance(j)) / (float)lookup[i]->getNumSeqs()))] = "";  
	//cout << "abundance  = " << lookup[i]->getAbundance(j) << '\t' << " relative adundance = " << (lookup[i]->getAbundance(j) / (float)lookup[i]->getNumSeqs()) << '\t';
	//cout << -log((log10(lookup[i]->getAbundance(j)) / lookup[i]->getNumSeqs())) << endl;
					//}else if (scaler == "log2") {
						//colorScale[-log((log2(lookup[i]->getAbundance(j)) / (float)lookup[i]->getNumSeqs()))] = "";  //cout << (int)log2(lookup[i]->getAbundance(j)) << endl;
	//cout << "abundance  = " << lookup[i]->getAbundance(j) << '\t' << " relative adundance = " << lookup[i]->getAbundance(j) / (float)lookup[i]->getNumSeqs() << '\t';
	//cout << -log((log2(lookup[i]->getAbundance(j)) / lookup[i]->getNumSeqs())) << endl;
				//	}else if (scaler == "linear") {
						colorScale[lookup[i]->getAbundance(j)] = "";
						if (maxbin < lookup[i]->getAbundance(j)) { maxbin = lookup[i]->getAbundance(j); }
	//cout << "abundance  = " << lookup[i]->getAbundance(j) << '\t' << " relative adundance = " << lookup[i]->getAbundance(j) / (float)lookup[i]->getNumSeqs() << '\t';
	//cout << lookup[i]->getAbundance(j) /(float) lookup[i]->getNumSeqs() << endl;
					//}else {  //if user enters invalid scaler option.
					//	cout << scaler << " is not a valid scaler option. I will use log10." << endl;
					//	colorScale[-log(log10(lookup[i]->getAbundance(j) / lookup[i]->getNumSeqs()))] = "";   
					//} 
				//}else { colorScale[0] = "00";  }
				
			}
		}
		
		//get scaler
		float scalers = 255 / (float) maxbin;
		
		//go through map and give each score a color value
		for (it = colorScale.begin(); it != colorScale.end(); it++) {
			it->second = toHex(int(float(it->first) * scalers));
			if(it->second.length() == 1) {  it->second = "0" + it->second;  }
//cout << it->first << " " << it->second << endl;
		}
		
		string filenamesvg = getRootName(globaldata->inputFileName) + sharedorder->getLabel() + ".heatmap.svg";
		openOutputFile(filenamesvg, outsvg);
		
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(lookup.size() * 300) + " " + toString((lookup[0]->getNumBins()*5 + 120))  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString(lookup.size() * 300) + "\" height=\"" + toString((lookup[0]->getNumBins()*5 + 120))  + "\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString((lookup.size() * 150) - 40) + "\" y=\"25\">Heatmap at distance " + sharedorder->getLabel() + "</text>\n";
		
		//column labels
		for (int h = 0; h < lookup.size(); h++) {
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(((300 * (h+1)) - 150) - ((int)lookup[h]->getGroup().length() / 2)) + "\" y=\"50\">" + lookup[h]->getGroup() + "</text>\n"; 
		}
		
		
		//output legend and color labels
		//go through map and give each score a color value
		string color;
		int x = 0;
		int y = 90 + (lookup[0]->getNumBins()*5);
		for (it = colorScale.begin(); it != colorScale.end(); it++) {
			color = it->second;	
			outsvg << "<rect fill=\"#" + color + "0000\" stroke=\"#" + color + "0000\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"25\" height=\"10\"/>\n";
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\">" + toString(it->first) + "</text>\n";
			x += 25;
		}
		
		x = 0;
		y = 70;
		
		//start at 1 since bin 0 is nothing
		for (int i = 1; i <= lookup[0]->getNumBins(); i++) {
			for (int j = 0; j < lookup.size(); j++) {
				
				//if (lookup[j]->getAbundance(i) != 0) { //don't want log value of 0.
					//if (scaler == "log10") {
					//	color = colorScale[(int)log10(lookup[j]->getAbundance(i))];  
					//}else if (scaler == "log2") {
					//	color = colorScale[(int)log2(lookup[j]->getAbundance(i))];  
				//	}else if (scaler == "linear") {
						color = colorScale[lookup[j]->getAbundance(i)]; 
					//}else {  color = colorScale[(int)log10(lookup[j]->getAbundance(i))];	} 
				//}else { color = "OO";  }

				
				outsvg << "<rect fill=\"#" + color + "0000\" stroke=\"#" + color + "0000\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"300\" height=\"5\"/>\n";
				x += 300;
			}
			x = 0;
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
void HeatMap::sortSharedVectors(){
	try {
		//copy lookup and then clear it to refill with sorted.
		//loop though lookup and determine if they are shared
		//if they are then insert in the front
		//if not push to back
		
		vector<SharedRAbundVector*> looktemp;
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		int count;
		
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
			temp->setLabel(looktemp[i]->getLabel());
			temp->setGroup(looktemp[i]->getGroup());
			lookup.push_back(temp); 
			
			//initialize place map
			place[i] = 0;
		}
		
		
		//for each bin
		for (int i = 0; i < looktemp[0]->size(); i++) {
			count = 0;
			bool updatePlace = false;
			//for each group
			for (int j = 0; j < looktemp.size(); j++) {
				if (looktemp[j]->getAbundance(i) != 0) { count++; }
			}
			
			//fill lookup
			for (int j = 0; j < looktemp.size(); j++) {
				//if they are not shared then push to back, if they are not insert in front
				if (count < 2)  { lookup[j]->push_back(looktemp[j]->getAbundance(i), i, looktemp[j]->getGroup()); }
				//they are shared by some
				else {  lookup[j]->insert(looktemp[j]->getAbundance(i), place[count], looktemp[j]->getGroup());   updatePlace = true; }
			}
			
			if (updatePlace == true) {
				//move place holders below where you entered up to "make space" for you entry
				for (it = place.begin(); it!= place.end(); it++) {  
					if (it->first < count) { it->second++; }
				}
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





