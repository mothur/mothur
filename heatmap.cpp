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

void HeatMap::getPic(RAbundVector* rabund) {
	try {
		//get users scaling method
		scaler = globaldata->getScale();
		
		float maxRelAbund = 0.0;		
		
		for(int i=0;i<rabund->size();i++){				
			float relAbund = rabund->get(i) / (float)rabund->getNumSeqs();
			if(relAbund > maxRelAbund){	maxRelAbund = relAbund;	}
		}
		
		scaler = globaldata->getScale();
		
		vector<string> scaleRelAbund(rabund->size(), "");
		
		for(int i=0;i<rabund->size();i++){
			float relAbund = rabund->get(i) / (float)rabund->getNumSeqs();
			
			if (rabund->get(i) != 0) { //don't want log value of 0.
				if (scaler == "log10") {
					scaleRelAbund[i] = toHex(int(255 * log10(relAbund) / log10(maxRelAbund))) + "0000";  
				}else if (scaler == "log2") {
					scaleRelAbund[i] = toHex(int(255 * log2(relAbund) / log2(maxRelAbund))) + "0000";  
				}else if (scaler == "linear") {
					scaleRelAbund[i] = toHex(int(255 * relAbund / maxRelAbund)) + "0000";  
				}else {  //if user enters invalid scaler option.
					scaleRelAbund[i] = toHex(int(255 * log10(relAbund / log10(maxRelAbund))))  + "0000"; 
				} 
			}
			else { scaleRelAbund[i] = "FFFFFF";  }
		}
		
		
		string filenamesvg = getRootName(globaldata->inputFileName) + rabund->getLabel() + ".heatmap.svg";
		openOutputFile(filenamesvg, outsvg);
		
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(300) + " " + toString((rabund->getNumBins()*5 + 120))  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString(300) + "\" height=\"" + toString((rabund->getNumBins()*5 + 120))  + "\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString((150) - 40) + "\" y=\"25\">Heatmap at distance " + rabund->getLabel() + "</text>\n";
				
		//output legend and color labels
		string color;
		int x = 0;
		int y = 103 + (rabund->getNumBins()*5);
		printLegend(y, maxRelAbund);
		
		y = 70;
		for (int i = 0; i < scaleRelAbund.size(); i++) {
			
			outsvg << "<rect fill=\"#" + scaleRelAbund[i] + "\" stroke=\"#" + scaleRelAbund[i] + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"300\" height=\"5\"/>\n";
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

void HeatMap::getPic(vector<SharedRAbundVector*> lookup) {
	try {
		//sort lookup so shared bins are on top
		if (sorted == "T") {  sortSharedVectors(lookup);  }
		
		vector<vector<string> > scaleRelAbund;
		vector<float> maxRelAbund(lookup.size(), 0.0);		
		float superMaxRelAbund = 0;
		
		for(int i = 0; i < lookup.size(); i++){
			for(int j=0; j<lookup[i]->size(); j++){
				
				float relAbund = lookup[i]->getAbundance(j) / (float)lookup[i]->getNumSeqs();
				if(relAbund > maxRelAbund[i]){	maxRelAbund[i] = relAbund;	}
			}
			if(maxRelAbund[i] > superMaxRelAbund){	superMaxRelAbund = maxRelAbund[i];	}
		}
		
		scaler = globaldata->getScale();
		
		scaleRelAbund.resize(lookup.size());
		for(int i=0;i<lookup.size();i++){
			scaleRelAbund[i].assign(lookup[i]->size(), "");
			for(int j=0;j<lookup[i]->size();j++){
				float relAbund = lookup[i]->getAbundance(j) / (float)lookup[i]->getNumSeqs();
				
				if (lookup[i]->getAbundance(j) != 0) { //don't want log value of 0.
					if (scaler == "log10") {
						scaleRelAbund[i][j] = toHex(int(255 * log10(relAbund) / log10(maxRelAbund[i]))) + "0000";  
					}else if (scaler == "log2") {
						scaleRelAbund[i][j] = toHex(int(255 * log2(relAbund) / log2(maxRelAbund[i]))) + "0000";  
					}else if (scaler == "linear") {
						scaleRelAbund[i][j] = toHex(int(255 * relAbund / maxRelAbund[i])) + "0000";  
					}else {  //if user enters invalid scaler option.
						scaleRelAbund[i][j] = toHex(int(255 * log10(relAbund / log10(maxRelAbund[i]))))  + "0000"; 
					} 
				}else { scaleRelAbund[i][j] = "FFFFFF";  }

			}
		}

		string filenamesvg = getRootName(globaldata->inputFileName) + lookup[0]->getLabel() + ".heatmap.svg";
		openOutputFile(filenamesvg, outsvg);
		
		//svg image
		outsvg << "<svg width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(lookup.size() * 300) + " " + toString((lookup[0]->getNumBins()*5 + 120))  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString(lookup.size() * 300) + "\" height=\"" + toString((lookup[0]->getNumBins()*5 + 120))  + "\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString((lookup.size() * 150) - 40) + "\" y=\"25\">Heatmap at distance " + lookup[0]->getLabel() + "</text>\n";
		
		//column labels
		for (int h = 0; h < lookup.size(); h++) {
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(((300 * (h+1)) - 150) - ((int)lookup[h]->getGroup().length() / 2)) + "\" y=\"50\">" + lookup[h]->getGroup() + "</text>\n"; 
		}
		
		//output legend and color labels
		string color;
		int x = 0;
		int y = 103 + (lookup[0]->getNumBins()*5);
		printLegend(y, superMaxRelAbund);
		
		y = 70;
		for (int i = 0; i < scaleRelAbund[0].size(); i++) {
			for (int j = 0; j < scaleRelAbund.size(); j++) {
				
				outsvg << "<rect fill=\"#" + scaleRelAbund[j][i] + "\" stroke=\"#" + scaleRelAbund[j][i] + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"300\" height=\"5\"/>\n";
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
void HeatMap::sortSharedVectors(vector<SharedRAbundVector*>& lookup){
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
		
		//clear out lookup to create sorted lookup -- Sarah look at - this is causing segmentation faults
		for (int j = 0; j < lookup.size(); j++) {
//			delete lookup[j];
		}
		lookup.clear();  //doesn't this do the job?
		
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
		
		//delete looktemp -- Sarah look at - this is causing segmentation faults
		for (int j = 0; j < looktemp.size(); j++) {
//			delete looktemp[j];
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

void HeatMap::printLegend(int y, float maxbin) {
	try {
		
		//output legend and color labels
		//go through map and give each score a color value
		string color;
		int x = 10;
		
		//prints legend
		for (int i = 1; i < 255; i++) {
			color = toHex(int((float)(i)));
			outsvg << "<rect fill=\"#" + color + "0000\" stroke=\"#" + color + "0000\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"1\" height=\"10\"/>\n";
			x += 1;
		}
		
		//prints legend labels
		x = 10;
		for (int i = 1; i<=5; i++) {
			float label;
			if(scaler== "log10")		{	label = maxbin * log10(51*i) / log10(255);	}
			else if(scaler== "log2")	{	label = maxbin * log2(51*i) / log2(255);	}
			else if(scaler== "linear")	{	label = maxbin * 51 * i / 255;				}
			else						{	label = maxbin * log10(51*i) / log10(255);	}
			file://localhost/Users/westcott/Desktop/c.amazon.fn.0.19.rep.fasta
			label = int(label * 1000 + 0.5);
			label /= 1000.0;
			string text = toString(label, 3);
			
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(x) + "\" y=\"" + toString(y-3) + "\">" + text + "</text>\n";
			x += 60;
		}
	}
	
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function printLegend. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function printLegend. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}




