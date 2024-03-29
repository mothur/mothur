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
HeatMap::HeatMap(string sort, string scale, int num, int fsize, string dir, string i){
	try {
		m = MothurOut::getInstance();
		sorted = sort;
		scaler = scale;
		outputDir = dir;
		numOTU = num;
		fontSize = fsize;
		inputfile = i;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "HeatMap");
		exit(1);
	}
}

//**********************************************************************************************************************

string HeatMap::getPic(RAbundVector* rabund) {
	try {
		
		int numBinsToDisplay = rabund->getNumBins();
		
		if (numOTU != 0) { //user want to display a portion of the otus
			if (numOTU < numBinsToDisplay) {  numBinsToDisplay = numOTU; }
		}
		
		//sort lookup so shared bins are on top
		if (sorted != "none") {  sortRabund(rabund);  }
		
		float maxRelAbund = 0.0;		
		
		for(int i=0;i<rabund->size();i++){				
			float relAbund = rabund->get(i) / (float)rabund->getNumSeqs();
			if(relAbund > maxRelAbund){	maxRelAbund = relAbund;	}
		}
		
		vector<string> scaleRelAbund(numBinsToDisplay, "");
		
		for(int i=0;i<numBinsToDisplay;i++){
			float relAbund = rabund->get(i) / (float)rabund->getNumSeqs();
			
			if (m->getControl_pressed()) { return "control"; }
			
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
		
		
		string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + rabund->getLabel() + ".heatmap.bin.svg";
		util.openOutputFile(filenamesvg, outsvg);
		
		//svg image
		outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(300) + " " + toString((numBinsToDisplay*5 + 120))  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString(300) + "\" height=\"" + toString((numBinsToDisplay*5 + 120))  + "\"/>"; 
		outsvg << "<text fill=\"black\" class=\"seri\" text-anchor=\"middle\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString((150) - 40) + "\" y=\"25\">Heatmap at distance " + rabund->getLabel() + "</text>\n";
				
		//output legend and color labels
		string color;
		int x = 0;
		int y = 103 + (numBinsToDisplay*5);
		printLegend(y, maxRelAbund);
		
		y = 70;

		for (int i = 0; i < scaleRelAbund.size(); i++) {
			if (m->getControl_pressed()) { outsvg.close(); return "control"; }
			
			outsvg << "<rect fill=\"#" + scaleRelAbund[i] + "\" stroke=\"#" + scaleRelAbund[i] + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"300\" height=\"5\"/>\n";
			y += 5;
		}
		
		outsvg << "</g>\n</svg>\n";
		outsvg.close();
		
		return filenamesvg;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "getPic");
		exit(1);
	}
}

//**********************************************************************************************************************

string HeatMap::getPic(SharedRAbundVectors*& data) {
	try {
        vector<SharedRAbundVector*> lookup = data->getSharedRAbundVectors();
        vector<string> groups = data->getNamesGroups();
        
		int numBinsToDisplay = lookup[0]->getNumBins();
		
		if (numOTU != 0) { //user want to display a portion of the otus
			if (numOTU < numBinsToDisplay) {  numBinsToDisplay = numOTU; }
		}
		
		//sort lookup so shared bins are on top
        vector<string> sortedLabels = data->getOTUNames();
		if (sorted != "none") {  sortedLabels = sortSharedVectors(lookup, sortedLabels);  }
		
		vector<vector<string> > scaleRelAbund;
		vector<float> maxRelAbund(lookup.size(), 0.0);
		float superMaxRelAbund = 0;
		
		for(int i = 0; i < lookup.size(); i++){
			for(int j=0; j<lookup[i]->size(); j++){
				
				float relAbund = lookup[i]->get(j) / (float)lookup[i]->getNumSeqs();
				if(relAbund > maxRelAbund[i]){	maxRelAbund[i] = relAbund;	}
			}
			if(maxRelAbund[i] > superMaxRelAbund){	superMaxRelAbund = maxRelAbund[i];	}
		}
		
		scaleRelAbund.resize(lookup.size());
		for(int i=0;i<lookup.size();i++){
			scaleRelAbund[i].assign(numBinsToDisplay, "");
			for(int j=0;j<numBinsToDisplay;j++){
				if (m->getControl_pressed()) {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } return "control"; }
				float relAbund = lookup[i]->get(j) / (float)lookup[i]->getNumSeqs();
				
				if (lookup[i]->get(j) != 0) { //don't want log value of 0.
					if (scaler == "log10") {
                        if (util.isEqual(maxRelAbund[i], 1)) { maxRelAbund[i] -= 0.001; }
                        if (util.isEqual(relAbund, 1)) { relAbund -= 0.001; }
						scaleRelAbund[i][j] = toHex(int(255 * log10(relAbund) / log10(maxRelAbund[i]))) + "0000";
					}else if (scaler == "log2") {
                        if (util.isEqual(maxRelAbund[i], 1)) { maxRelAbund[i] -= 0.001; }
                        if (util.isEqual(relAbund, 1)) { relAbund -= 0.001; }
						scaleRelAbund[i][j] = toHex(int(255 * log2(relAbund) / log2(maxRelAbund[i]))) + "0000";  
					}else if (scaler == "linear") {
						scaleRelAbund[i][j] = toHex(int(255 * relAbund / maxRelAbund[i])) + "0000";
					}else {  //if user enters invalid scaler option.
                        if (util.isEqual(maxRelAbund[i], 1)) { maxRelAbund[i] += 0.001; }
						scaleRelAbund[i][j] = toHex(int(255 * log10(relAbund / log10(maxRelAbund[i]))))  + "0000"; 
					} 
				}else { scaleRelAbund[i][j] = "FFFFFF";  }

			}
		}

		string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + ".heatmap.bin.svg";
		util.openOutputFile(filenamesvg, outsvg);
        int binHeight = 20;
        int labelBump = 100;
        int binWidth = 300;
		
		//svg image
		outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(lookup.size() * binWidth + labelBump) + " " + toString((numBinsToDisplay*binHeight + 120))  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString(lookup.size() * binWidth+labelBump) + "\" height=\"" + toString((numBinsToDisplay*binHeight + 120))  + "\"/>";
		outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" text-anchor=\"middle\" x=\"" + toString((lookup.size() * 150) - 40) + "\" y=\"25\">Heatmap at distance " + lookup[0]->getLabel() + "</text>\n";
		
		//column labels
		for (int h = 0; h < lookup.size()+1; h++) {
            if (h == 0) {
                string tempLabel = "OTU";
                outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(labelBump-labelBump/2+1) + "\" y=\"50\">" + tempLabel + "</text>\n";
            }else {
                outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(((binWidth * h) - 150) - ((int)groups[h-1].length() / 2)+labelBump/2) + "\" y=\"50\">" + groups[h-1] + "</text>\n";
            }
		}

		//output legend and color labels
		string color;
		int x = 0;
		int y = 103 + (numBinsToDisplay*binHeight);
		printLegend(y, superMaxRelAbund);
		
		y = 70;
		for (int i = 0; i < numBinsToDisplay; i++) {
            outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\">" + sortedLabels[i] + "</text>\n";
            x += labelBump;
			for (int j = 0; j < scaleRelAbund.size(); j++) {
				if (m->getControl_pressed()) { outsvg.close(); return "control"; }
				
				outsvg << "<rect fill=\"#" + scaleRelAbund[j][i] + "\" stroke=\"#" + scaleRelAbund[j][i] + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"" + toString(binWidth) +  "\" height=\"" + toString(binHeight) +  "\"/>\n";
				x += binWidth;
			}
			x = 0;
			y += binHeight;
		}
		
		outsvg << "</g>\n</svg>\n";
		outsvg.close();
		
        for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
        
		return filenamesvg;

	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "getPic");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> HeatMap::sortSharedVectors(vector<SharedRAbundVector*> lookup, vector<string> currentLabels){
	try {
        vector<string> sortedLabels; sortedLabels.resize(currentLabels.size(), "");
		vector<SharedRAbundVector*> looktemp;
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
        
		/****************** find order of otus **********************/
		if (sorted == "shared") {
			place = orderShared(lookup);	
		}else if (sorted == "topotu") {
			place = orderTopOtu(lookup);	
		}else if (sorted == "topgroup") {
			place = orderTopGroup(lookup);	
		}else { m->mothurOut("Error: invalid sort option.\n");    }
				
		
		/******************* create copy of lookup *********************/
		//create and initialize looktemp as a copy of lookup
		for (int i = 0; i < lookup.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(*lookup[i]);
			temp->setLabel(lookup[i]->getLabel());
			looktemp.push_back(temp);
		}
	
		/************************ fill lookup in order given by place *********************/
		//for each bin
		for (int i = 0; i < looktemp[0]->getNumBins(); i++) {														//place
			//fill lookup																					// 2 -> 1
			for (int j = 0; j < looktemp.size(); j++) {														// 3 -> 2
				int newAbund = looktemp[j]->get(i);												// 1 -> 3
				lookup[j]->set(place[i], newAbund); //binNumber, abundance, group
                sortedLabels[place[i]] = currentLabels[i];
			}
		}
		
        return sortedLabels;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "sortSharedVectors");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, int> HeatMap::orderShared(vector<SharedRAbundVector*> lookup){
	try {
				
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		vector<int> sharedBins;
		vector<int> uniqueBins;
		
		//for each bin
		for (int i = 0; i < lookup[0]->getNumBins(); i++) {
			int count = 0;												
			
			//is this bin shared
			for (int j = 0; j < lookup.size(); j++) {		if (lookup[j]->get(i) != 0) { count++; }	}
			
			if (count < 2)	{  uniqueBins.push_back(i); }
			else			{  sharedBins.push_back(i); }
		}
		
		//fill place
		for (int i = 0; i < sharedBins.size(); i++) { place[sharedBins[i]] = i; }
		for (int i = 0; i < uniqueBins.size(); i++) { place[uniqueBins[i]] = (sharedBins.size() + i); }
		
		return place;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "orderShared");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, int> HeatMap::orderTopOtu(vector<SharedRAbundVector*> lookup){
	try {
				
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		vector<binCount> totals;
		
		//for each bin
		for (int i = 0; i < lookup[0]->getNumBins(); i++) {
			int total = 0;												
			
			for (int j = 0; j < lookup.size(); j++) {	total += lookup[j]->get(i); 	}
			
			binCount temp(i, total);
			
			totals.push_back(temp);
		}
		
		sort(totals.begin(), totals.end(), comparebinCounts);
		
		//fill place
		for (int i = 0; i < totals.size(); i++) {   place[totals[i].bin] = i;  }
				
		return place;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "orderTopOtu");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, int> HeatMap::orderTopGroup(vector<SharedRAbundVector*> lookup){
	try {
				
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		vector < vector<binCount> > totals; //totals[0] = bin totals for group 0, totals[1] = bin totals for group 1, ...
		totals.resize(lookup.size());
		
		//for each bin
		for (int i = 0; i < lookup[0]->getNumBins(); i++) {
			for (int j = 0; j < lookup.size(); j++) {
				binCount temp(i, (lookup[j]->get(i)));
				totals[j].push_back(temp);
			}
		}
		
		for (int i = 0; i < totals.size(); i++) { sort(totals[i].begin(), totals[i].end(), comparebinCounts);  }
		
		//fill place
		//grab the top otu for each group adding it if its not already added
		int count = 0;
		for (int i = 0; i < totals[0].size(); i++) { 
		
			for (int j = 0; j < totals.size(); j++) {  
				it = place.find(totals[j][i].bin);
				
				if (it == place.end()) { //not added yet
					place[totals[j][i].bin] = count;
					count++;
				}
			}
		}
				
		return place;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "orderTopGroup");
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
			label = int(label * 1000 + 0.5);
			label /= 1000.0;
			string text = toString(label);
			
			outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString(x) + "\" y=\"" + toString(y-3) + "\">" + text + "</text>\n";
			x += 60;
		}
	}
	
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "printLegend");
		exit(1);
	}
}
//**********************************************************************************************************************

string HeatMap::getPic(SharedRAbundFloatVectors*& data) {
	try {
        
        vector<SharedRAbundFloatVector*> lookup = data->getSharedRAbundFloatVectors();
        vector<string> groups = data->getNamesGroups();
        
		int numBinsToDisplay = lookup[0]->getNumBins();
		
		if (numOTU != 0) { //user want to display a portion of the otus
			if (numOTU < numBinsToDisplay) {  numBinsToDisplay = numOTU; }
		}
		
		//sort lookup so shared bins are on top
        vector<string> sortedLabels = data->getOTUNames();
		if (sorted != "none") {  sortedLabels = sortSharedVectors(lookup, sortedLabels);  }
		
		vector<vector<string> > scaleRelAbund;
		vector<float> maxRelAbund(lookup.size(), 0.0);		
		float superMaxRelAbund = 0;
		
		for(int i = 0; i < lookup.size(); i++){
			for(int j=0; j<numBinsToDisplay; j++){
				
				float relAbund = lookup[i]->get(j);
				if(relAbund > maxRelAbund[i]){	maxRelAbund[i] = relAbund;	}
			}
			if(maxRelAbund[i] > superMaxRelAbund){	superMaxRelAbund = maxRelAbund[i];	}
		}
		
		scaleRelAbund.resize(lookup.size());
		for(int i=0;i<lookup.size();i++){
			scaleRelAbund[i].assign(numBinsToDisplay, "");
			for(int j=0;j<numBinsToDisplay;j++){
				if (m->getControl_pressed()) {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } return "control"; }
				float relAbund = lookup[i]->get(j);
				
				if (!util.isEqual(relAbund, 0)) { //don't want log value of 0.
					if (scaler == "log10") {
                        if (util.isEqual(maxRelAbund[i], 1)) { maxRelAbund[i] -= 0.001; }
                        if (util.isEqual(relAbund, 1)) { relAbund -= 0.001; }
						scaleRelAbund[i][j] = toHex(int(255 * log10(relAbund) / log10(maxRelAbund[i]))) + "0000";  
					}else if (scaler == "log2") {
                        if (util.isEqual(maxRelAbund[i], 1)) { maxRelAbund[i] -= 0.001; }
                        if (util.isEqual(relAbund, 1)) { relAbund -= 0.001; }
						scaleRelAbund[i][j] = toHex(int(255 * log2(relAbund) / log2(maxRelAbund[i]))) + "0000";  
					}else if (scaler == "linear") {
						scaleRelAbund[i][j] = toHex(int(255 * relAbund / maxRelAbund[i])) + "0000";
					}else {  //if user enters invalid scaler option.
						scaleRelAbund[i][j] = toHex(int(255 * log10(relAbund / log10(maxRelAbund[i]))))  + "0000"; 
					} 
				}else { scaleRelAbund[i][j] = "FFFFFF";  }

			}
		}

		string filenamesvg = outputDir + util.getRootName(util.getSimpleName(inputfile)) + lookup[0]->getLabel() + ".heatmap.bin.svg";
		util.openOutputFile(filenamesvg, outsvg);
        
        int binHeight = 20;
        int labelBump = 100;
        int binWidth = 300;
		
		//svg image
		outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 " + toString(lookup.size() * binWidth + labelBump) + " " + toString((numBinsToDisplay*binHeight + 120))  + "\">\n";
		outsvg << "<g>\n";
		
		//white backround
		outsvg << "<rect fill=\"white\" stroke=\"white\" x=\"0\" y=\"0\" width=\"" + toString(lookup.size() * binWidth+labelBump) + "\" height=\"" + toString((numBinsToDisplay*binHeight + 120))  + "\"/>";
		outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" text-anchor=\"middle\" x=\"" + toString((lookup.size() * 150) - 40) + "\" y=\"25\">Heatmap at distance " + lookup[0]->getLabel() + "</text>\n";
		
		//column labels
		for (int h = 0; h < lookup.size()+1; h++) {
            if (h == 0) {
                string tempLabel = "OTU";
                outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(labelBump-labelBump/2+1) + "\" y=\"50\">" + tempLabel + "</text>\n";
            }else {
                outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(((binWidth * h) - 150) - ((int)groups[h-1].length() / 2)+labelBump/2) + "\" y=\"50\">" + groups[h-1] + "</text>\n";
            }
		}
        
		//output legend and color labels
		string color;
		int x = 0;
		int y = 103 + (numBinsToDisplay*binHeight);
		printLegend(y, superMaxRelAbund);
		
		y = 70;
		for (int i = 0; i < numBinsToDisplay; i++) {
            outsvg << "<text fill=\"black\" class=\"seri\" font-size=\"" + toString(fontSize) + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\">" + sortedLabels[i] + "</text>\n";
            x += labelBump;
			for (int j = 0; j < scaleRelAbund.size(); j++) {
				if (m->getControl_pressed()) { outsvg.close(); return "control"; }
				
				outsvg << "<rect fill=\"#" + scaleRelAbund[j][i] + "\" stroke=\"#" + scaleRelAbund[j][i] + "\" x=\"" + toString(x) + "\" y=\"" + toString(y) + "\" width=\"" + toString(binWidth) +  "\" height=\"" + toString(binHeight) +  "\"/>\n";
				x += binWidth;
			}
			x = 0;
			y += binHeight;
		}
		
		outsvg << "</g>\n</svg>\n";
		outsvg.close();
        
        for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		
		return filenamesvg;

	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "getPic");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> HeatMap::sortSharedVectors(vector<SharedRAbundFloatVector*> lookup, vector<string> currentLabels){
	try {
				
		vector<SharedRAbundFloatVector*> looktemp;
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
        
        vector<string> sortedLabels; sortedLabels.resize(currentLabels.size(), "");
		
		/****************** find order of otus **********************/
		if (sorted == "shared") {
			place = orderShared(lookup);	
		}else if (sorted == "topotu") {
			place = orderTopOtu(lookup);	
		}else if (sorted == "topgroup") {
			place = orderTopGroup(lookup);	
		}else { m->mothurOut("Error: invalid sort option.\n");   return sortedLabels; }
				
		
		/******************* create copy of lookup *********************/
		//create and initialize looktemp as a copy of lookup
		for (int i = 0; i < lookup.size(); i++) { 
			SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(*lookup[i]);
			temp->setLabel(lookup[i]->getLabel());
            looktemp.push_back(temp);
		}
	
		/************************ fill lookup in order given by place *********************/
		//for each bin
		for (int i = 0; i < looktemp[0]->size(); i++) {														//place
			//fill lookup																					// 2 -> 1
			for (int j = 0; j < looktemp.size(); j++) {														// 3 -> 2
				float newAbund = looktemp[j]->get(i);												// 1 -> 3
				lookup[j]->set(place[i], newAbund); //binNumber, abundance, group
                sortedLabels[place[i]] = currentLabels[i];
			}
		}
		
		return sortedLabels;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "sortSharedVectors");
		exit(1);
	}
}
//**********************************************************************************************************************
int HeatMap::sortRabund(RAbundVector* r){
	try {
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		/****************** find order of otus **********************/
		vector<binCount> totals;
		
		//for each bin
		for (int i = 0; i < r->getNumBins(); i++) {	
			binCount temp(i, r->get(i));
			
			totals.push_back(temp);
		}
		
		sort(totals.begin(), totals.end(), comparebinCounts);
		
		//fill place
		for (int i = 0; i < totals.size(); i++) {   place[totals[i].bin] = i;  }
		
		/******************* create copy of lookup *********************/
		//create and initialize rtemp as a copy of r
		
		RAbundVector* rtemp = new RAbundVector(r->getNumBins());
		for (int i = 0; i < r->size(); i++) {  rtemp->set(i, r->get(i)); }
		rtemp->setLabel(r->getLabel());
			
		/************************ fill lookup in order given by place *********************/
		//for each bin
		for (int i = 0; i < rtemp->size(); i++) {										//place
			//fill lookup																// 2 -> 1
																						// 3 -> 2
			int newAbund = rtemp->get(i);												// 1 -> 3
			r->set(place[i], newAbund); //binNumber, abundance
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "sortRabund");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, int> HeatMap::orderShared(vector<SharedRAbundFloatVector*> lookup){
	try {
				
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		vector<int> sharedBins;
		vector<int> uniqueBins;
		
		//for each bin
		for (int i = 0; i < lookup[0]->getNumBins(); i++) {
			int count = 0;												
			
			//is this bin shared
			for (int j = 0; j < lookup.size(); j++) {
                if (!util.isEqual(lookup[j]->get(i), 0)) { count++; }
            }
			
			if (count < 2)	{  uniqueBins.push_back(i); }
			else			{  sharedBins.push_back(i); }
		}
		
		//fill place
		for (int i = 0; i < sharedBins.size(); i++) { place[sharedBins[i]] = i; }
		for (int i = 0; i < uniqueBins.size(); i++) { place[uniqueBins[i]] = (sharedBins.size() + i); }
		
		return place;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "orderShared");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, int> HeatMap::orderTopOtu(vector<SharedRAbundFloatVector*> lookup){
	try {
				
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		vector<binCountFloat> totals;
		
		//for each bin
		for (int i = 0; i < lookup[0]->size(); i++) {	
			int total = 0;												
			
			for (int j = 0; j < lookup.size(); j++) {	total += lookup[j]->get(i); 	}
			
			binCountFloat temp(i, total);
			
			totals.push_back(temp);
		}
		
		sort(totals.begin(), totals.end(), comparebinFloatCounts);
		
		//fill place
		for (int i = 0; i < totals.size(); i++) {   place[totals[i].bin] = i;  }
				
		return place;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "orderTopOtu");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, int> HeatMap::orderTopGroup(vector<SharedRAbundFloatVector*> lookup){
	try {
				
		map<int, int> place; //spot in lookup where you insert shared by, ie, 3 -> 2 if they are shared by 3 inset into location 2.
		map<int, int>::iterator it;
		
		vector < vector<binCountFloat> > totals; //totals[0] = bin totals for group 0, totals[1] = bin totals for group 1, ...
		totals.resize(lookup.size());
		
		//for each bin
		for (int i = 0; i < lookup[0]->size(); i++) {	
			for (int j = 0; j < lookup.size(); j++) {
				binCountFloat temp(i, (lookup[j]->get(i)));
				totals[j].push_back(temp);
			}
		}
		
		for (int i = 0; i < totals.size(); i++) { sort(totals[i].begin(), totals[i].end(), comparebinFloatCounts);  }
		
		//fill place
		//grab the top otu for each group adding it if its not already added
		int count = 0;
		for (int i = 0; i < totals[0].size(); i++) { 
		
			for (int j = 0; j < totals.size(); j++) {  
				it = place.find(totals[j][i].bin);
				
				if (it == place.end()) { //not added yet
					place[totals[j][i].bin] = count;
					count++;
				}
			}
		}
				
		return place;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMap", "orderTopGroup");
		exit(1);
	}
}
//**********************************************************************************************************************





