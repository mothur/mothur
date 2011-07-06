/*
 *  bayesian.cpp
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "bayesian.h"
#include "kmer.hpp"
#include "phylosummary.h"
#include "referencedb.h"
/**************************************************************************************************/
Bayesian::Bayesian(string tfile, string tempFile, string method, int ksize, int cutoff, int i) : 
Classify(), kmerSize(ksize), confidenceThreshold(cutoff), iters(i)  {
	try {
		ReferenceDB* rdb = ReferenceDB::getInstance();
		
		string baseName = tempFile;
		if (baseName == "saved") { baseName = rdb->getSavedReference(); }
		
		string baseTName = tfile;
		if (baseTName == "saved") { baseTName = rdb->getSavedTaxonomy(); }
		
		/************calculate the probablity that each word will be in a specific taxonomy*************/
		string tfileroot = baseTName.substr(0,baseTName.find_last_of(".")+1);
		string tempfileroot = m->getRootName(m->getSimpleName(baseName));
		string phyloTreeName = tfileroot + "tree.train";
		string phyloTreeSumName = tfileroot + "tree.sum";
		string probFileName = tfileroot + tempfileroot + char('0'+ kmerSize) + "mer.prob";
		string probFileName2 = tfileroot + tempfileroot + char('0'+ kmerSize) + "mer.numNonZero";
		
		ofstream out;
		ofstream out2;
		
		ifstream phyloTreeTest(phyloTreeName.c_str());
		ifstream probFileTest2(probFileName2.c_str());
		ifstream probFileTest(probFileName.c_str());
		ifstream probFileTest3(phyloTreeSumName.c_str());
		
		int start = time(NULL);
		
		//if they are there make sure they were created after this release date
		bool FilesGood = false;
		if(probFileTest && probFileTest2 && phyloTreeTest && probFileTest3){
			FilesGood = checkReleaseDate(probFileTest, probFileTest2, phyloTreeTest, probFileTest3);
		}
		
		//if you want to save, but you dont need to calculate then just read
		if (rdb->save && probFileTest && probFileTest2 && phyloTreeTest && probFileTest3 && FilesGood) {  
			ifstream saveIn;
			m->openInputFile(tempFile, saveIn);
			
			while (!saveIn.eof()) {
				Sequence temp(saveIn);
				m->gobble(saveIn);
				
				rdb->referenceSeqs.push_back(temp); 
			}
			saveIn.close();			
		}
		
		if(probFileTest && probFileTest2 && phyloTreeTest && probFileTest3 && FilesGood){	
			if (tempFile == "saved") { m->mothurOutEndLine();  m->mothurOut("Using sequences from " + rdb->getSavedReference() + " that are saved in memory.");	m->mothurOutEndLine(); }
			
			m->mothurOut("Reading template taxonomy...     "); cout.flush();
			
			phyloTree = new PhyloTree(phyloTreeTest, phyloTreeName);
			
			m->mothurOut("DONE."); m->mothurOutEndLine();
			
			genusNodes = phyloTree->getGenusNodes(); 
			genusTotals = phyloTree->getGenusTotals();
			
			if (tfile == "saved") { 
				m->mothurOutEndLine();  m->mothurOut("Using probabilties from " + rdb->getSavedTaxonomy() + " that are saved in memory...    ");	cout.flush();; 
				wordGenusProb = rdb->wordGenusProb;
			}else {
				m->mothurOut("Reading template probabilities...     "); cout.flush();
				readProbFile(probFileTest, probFileTest2, probFileName, probFileName2);
			}	
			
			//save probabilities
			if (rdb->save) { rdb->wordGenusProb = wordGenusProb; }
		}else{
		
			//create search database and names vector
			generateDatabaseAndNames(tfile, tempFile, method, ksize, 0.0, 0.0, 0.0, 0.0);
			
			//prevents errors caused by creating shortcut files if you had an error in the sanity check.
			if (m->control_pressed) {  remove(phyloTreeName.c_str());  remove(probFileName.c_str()); remove(probFileName2.c_str()); }
			else{ 
				genusNodes = phyloTree->getGenusNodes(); 
				genusTotals = phyloTree->getGenusTotals();
				
				m->mothurOut("Calculating template taxonomy tree...     "); cout.flush();
				
				phyloTree->printTreeNodes(phyloTreeName);
							
				m->mothurOut("DONE."); m->mothurOutEndLine();
				
				m->mothurOut("Calculating template probabilities...     "); cout.flush();
				
				numKmers = database->getMaxKmer() + 1;
			
				//initialze probabilities
				wordGenusProb.resize(numKmers);
			//cout << numKmers << '\t' << genusNodes.size() << endl;
				for (int j = 0; j < wordGenusProb.size(); j++) {	wordGenusProb[j].resize(genusNodes.size());		}
			//cout << numKmers << '\t' << genusNodes.size() << endl;	
				ofstream out;
				ofstream out2;
				
				#ifdef USE_MPI
					int pid;
					MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

					if (pid == 0) {  
				#endif

				
				m->openOutputFile(probFileName, out);
				
				//output mothur version
				out << "#" << m->getVersion() << endl;
				
				out << numKmers << endl;
				
				m->openOutputFile(probFileName2, out2);
				
				//output mothur version
				out2 << "#" << m->getVersion() << endl;
				
				#ifdef USE_MPI
					}
				#endif

				
				//for each word
				for (int i = 0; i < numKmers; i++) {
					if (m->control_pressed) {  break; }
					
					#ifdef USE_MPI
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

						if (pid == 0) {  
					#endif

					out << i << '\t';
					
					#ifdef USE_MPI
						}
					#endif
					
					vector<int> seqsWithWordi = database->getSequencesWithKmer(i);
					
					map<int, int> count;
					for (int k = 0; k < genusNodes.size(); k++) {  count[genusNodes[k]] = 0;  }			
							
					//for each sequence with that word
					for (int j = 0; j < seqsWithWordi.size(); j++) {
						int temp = phyloTree->getIndex(names[seqsWithWordi[j]]);
						count[temp]++;  //increment count of seq in this genus who have this word
					}
					
					//probabilityInTemplate = (# of seqs with that word in template + 0.50) / (total number of seqs in template + 1);
					float probabilityInTemplate = (seqsWithWordi.size() + 0.50) / (float) (names.size() + 1);
					
					int numNotZero = 0;
					for (int k = 0; k < genusNodes.size(); k++) {
						//probabilityInThisTaxonomy = (# of seqs with that word in this taxonomy + probabilityInTemplate) / (total number of seqs in this taxonomy + 1);
						
						
						wordGenusProb[i][k] = log((count[genusNodes[k]] + probabilityInTemplate) / (float) (genusTotals[k] + 1));  
									
						if (count[genusNodes[k]] != 0) { 
							#ifdef USE_MPI
								int pid;
								MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
						
								if (pid == 0) {  
							#endif

							out << k << '\t' << wordGenusProb[i][k] << '\t'; 
							
							#ifdef USE_MPI
								}
							#endif

							numNotZero++;  
						}
					}
					
					#ifdef USE_MPI
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
						if (pid == 0) {  
					#endif
					
					out << endl;
					out2 << probabilityInTemplate << '\t' << numNotZero << endl;
					
					#ifdef USE_MPI
						}
					#endif
				}
				
				#ifdef USE_MPI
					MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
					if (pid == 0) {  
				#endif
				
				out.close();
				out2.close();
				
				#ifdef USE_MPI
					}
				#endif
				
				//read in new phylotree with less info. - its faster
				ifstream phyloTreeTest(phyloTreeName.c_str());
				delete phyloTree;
				
				phyloTree = new PhyloTree(phyloTreeTest, phyloTreeName);
				
				//save probabilities
				if (rdb->save) { rdb->wordGenusProb = wordGenusProb; }
			}
		}
	
		m->mothurOut("DONE."); m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " seconds get probabilities. "); m->mothurOutEndLine();
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "Bayesian");
		exit(1);
	}
}
/**************************************************************************************************/
Bayesian::~Bayesian() {
	try {
		
		 delete phyloTree; 
		 if (database != NULL) {  delete database; }
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "~Bayesian");
		exit(1);
	}
}

/**************************************************************************************************/
string Bayesian::getTaxonomy(Sequence* seq) {
	try {
		string tax = "";
		Kmer kmer(kmerSize);
		
		//get words contained in query
		//getKmerString returns a string where the index in the string is hte kmer number 
		//and the character at that index can be converted to be the number of times that kmer was seen
		
		string queryKmerString = kmer.getKmerString(seq->getUnaligned()); 

		vector<int> queryKmers;
		for (int i = 0; i < queryKmerString.length()-1; i++) {	// the -1 is to ignore any kmer with an N in it
			if (queryKmerString[i] != '!') { //this kmer is in the query
				queryKmers.push_back(i);
			}
		}
		
		if (queryKmers.size() == 0) {  m->mothurOut(seq->getName() + "is bad."); m->mothurOutEndLine(); return "bad seq"; }
		
		
		int index = getMostProbableTaxonomy(queryKmers);
		
		if (m->control_pressed) { return tax; }
					
		//bootstrap - to set confidenceScore
		int numToSelect = queryKmers.size() / 8;
	
		tax = bootstrapResults(queryKmers, index, numToSelect);
				
		return tax;	
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "getTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/
string Bayesian::bootstrapResults(vector<int> kmers, int tax, int numToSelect) {
	try {
				
		map<int, int> confidenceScores; 
		
		//initialize confidences to 0 
		int seqIndex = tax;
		TaxNode seq = phyloTree->get(tax);
		confidenceScores[tax] = 0;
		
		while (seq.level != 0) { //while you are not at the root
			seqIndex = seq.parent;
			confidenceScores[seqIndex] = 0;
			seq = phyloTree->get(seq.parent);
		}
				
		map<int, int>::iterator itBoot;
		map<int, int>::iterator itBoot2;
		map<int, int>::iterator itConvert;
			
		for (int i = 0; i < iters; i++) {
			if (m->control_pressed) { return "control"; }
			
			vector<int> temp;
			for (int j = 0; j < numToSelect; j++) {
				int index = int(rand() % kmers.size());
				
				//add word to temp
				temp.push_back(kmers[index]);
			}
			
			//get taxonomy
			int newTax = getMostProbableTaxonomy(temp);
			//int newTax = 1;
			TaxNode taxonomyTemp = phyloTree->get(newTax);
			
			//add to confidence results
			while (taxonomyTemp.level != 0) { //while you are not at the root
				itBoot2 = confidenceScores.find(newTax); //is this a classification we already have a count on
				
				if (itBoot2 != confidenceScores.end()) { //this is a classification we need a confidence for
					(itBoot2->second)++;
				}
				
				newTax = taxonomyTemp.parent;
				taxonomyTemp = phyloTree->get(newTax);
			}
	
		}
		
		string confidenceTax = "";
		simpleTax = "";
		
		int seqTaxIndex = tax;
		TaxNode seqTax = phyloTree->get(tax);
		
		while (seqTax.level != 0) { //while you are not at the root
					
				itBoot2 = confidenceScores.find(seqTaxIndex); //is this a classification we already have a count on
				
				int confidence = 0;
				if (itBoot2 != confidenceScores.end()) { //already in confidence scores
					confidence = itBoot2->second;
				}
				
				if (((confidence/(float)iters) * 100) >= confidenceThreshold) {
					confidenceTax = seqTax.name + "(" + toString(((confidence/(float)iters) * 100)) + ");" + confidenceTax;
					simpleTax = seqTax.name + ";" + simpleTax;
				}
				
				seqTaxIndex = seqTax.parent;
				seqTax = phyloTree->get(seqTax.parent);
		}
		
		if (confidenceTax == "") { confidenceTax = "unclassified;"; simpleTax = "unclassified;"; }
		return confidenceTax;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "bootstrapResults");
		exit(1);
	}
}
/**************************************************************************************************/
int Bayesian::getMostProbableTaxonomy(vector<int> queryKmer) {
	try {
		int indexofGenus = 0;
		
		double maxProbability = -1000000.0;
		//find taxonomy with highest probability that this sequence is from it
		
		
//		cout << genusNodes.size() << endl;
		
		
		for (int k = 0; k < genusNodes.size(); k++) {
			//for each taxonomy calc its probability
			
			double prob = 0.0000;
			for (int i = 0; i < queryKmer.size(); i++) {
				prob += wordGenusProb[queryKmer[i]][k];
			}
			
//			cout << phyloTree->get(genusNodes[k]).name << '\t' << prob << endl;

			//is this the taxonomy with the greatest probability?
			if (prob > maxProbability) { 
				indexofGenus = genusNodes[k];
				maxProbability = prob;
			}
		}
		
			
		return indexofGenus;
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "getMostProbableTaxonomy");
		exit(1);
	}
}
/*************************************************************************************************
map<string, int> Bayesian::parseTaxMap(string newTax) {
	try{
	
		map<string, int> parsed;
		
		newTax = newTax.substr(0, newTax.length()-1);  //get rid of last ';'
	
		//parse taxonomy
		string individual;
		while (newTax.find_first_of(';') != -1) {
			individual = newTax.substr(0,newTax.find_first_of(';'));
			newTax = newTax.substr(newTax.find_first_of(';')+1, newTax.length());
			parsed[individual] = 1;
		}
		
		//get last one
		parsed[newTax] = 1;

		return parsed;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "parseTax");
		exit(1);
	}
}
/**************************************************************************************************/
void Bayesian::readProbFile(ifstream& in, ifstream& inNum, string inName, string inNumName) {
	try{
		
		#ifdef USE_MPI
			
			int pid, num, num2, processors;
			vector<unsigned long int> positions;
			vector<unsigned long int> positions2;
			
			MPI_Status status; 
			MPI_File inMPI;
			MPI_File inMPI2;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
			MPI_Comm_size(MPI_COMM_WORLD, &processors);
			int tag = 2001;

			char inFileName[1024];
			strcpy(inFileName, inNumName.c_str());
			
			char inFileName2[1024];
			strcpy(inFileName2, inName.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
			MPI_File_open(MPI_COMM_WORLD, inFileName2, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI2);  //comm, filename, mode, info, filepointer

			if (pid == 0) {
				positions = m->setFilePosEachLine(inNumName, num);
				positions2 = m->setFilePosEachLine(inName, num2);
				
				for(int i = 1; i < processors; i++) { 
					MPI_Send(&num, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&positions[0], (num+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
					
					MPI_Send(&num2, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&positions2[0], (num2+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
				}

			}else{
				MPI_Recv(&num, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				positions.resize(num+1);
				MPI_Recv(&positions[0], (num+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
				
				MPI_Recv(&num2, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				positions2.resize(num2+1);
				MPI_Recv(&positions2[0], (num2+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
			}
			
			//read version
			int length = positions2[1] - positions2[0];
			char* buf5 = new char[length];

			MPI_File_read_at(inMPI2, positions2[0], buf5, length, MPI_CHAR, &status);
			delete buf5;

			//read numKmers
			length = positions2[2] - positions2[1];
			char* buf = new char[length];

			MPI_File_read_at(inMPI2, positions2[1], buf, length, MPI_CHAR, &status);

			string tempBuf = buf;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
			delete buf;

			istringstream iss (tempBuf,istringstream::in);
			iss >> numKmers;  
			
			//initialze probabilities
			wordGenusProb.resize(numKmers);
			
			for (int j = 0; j < wordGenusProb.size(); j++) {	wordGenusProb[j].resize(genusNodes.size());		}
			
			int kmer, name;  
			vector<int> numbers; numbers.resize(numKmers);
			float prob;
			vector<float> zeroCountProb; zeroCountProb.resize(numKmers);	
			
			//read version
			length = positions[1] - positions[0];
			char* buf6 = new char[length];

			MPI_File_read_at(inMPI2, positions[0], buf6, length, MPI_CHAR, &status);
			delete buf6;
			
			//read file 
			for(int i=1;i<num;i++){
				//read next sequence
				length = positions[i+1] - positions[i];
				char* buf4 = new char[length];

				MPI_File_read_at(inMPI, positions[i], buf4, length, MPI_CHAR, &status);

				tempBuf = buf4;
				if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
				delete buf4;

				istringstream iss (tempBuf,istringstream::in);
				iss >> zeroCountProb[i] >> numbers[i];  
			}
			
			MPI_File_close(&inMPI);
			
			for(int i=2;i<num2;i++){
				//read next sequence
				length = positions2[i+1] - positions2[i];
				char* buf4 = new char[length];

				MPI_File_read_at(inMPI2, positions2[i], buf4, length, MPI_CHAR, &status);

				tempBuf = buf4;
				if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
				delete buf4;

				istringstream iss (tempBuf,istringstream::in);
				
				iss >> kmer;
				
				//set them all to zero value
				for (int i = 0; i < genusNodes.size(); i++) {
					wordGenusProb[kmer][i] = log(zeroCountProb[kmer] / (float) (genusTotals[i]+1));
				}
				
				//get probs for nonzero values
				for (int i = 0; i < numbers[kmer]; i++) {
					iss >> name >> prob;
					wordGenusProb[kmer][name] = prob;
				}
				
			}
			MPI_File_close(&inMPI2);
			MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
		#else
			//read version
			string line = m->getline(in); m->gobble(in);
			
			in >> numKmers; m->gobble(in);
			
			//initialze probabilities
			wordGenusProb.resize(numKmers);
			
			for (int j = 0; j < wordGenusProb.size(); j++) {	wordGenusProb[j].resize(genusNodes.size());		}
			
			int kmer, name, count;  count = 0;
			vector<int> num; num.resize(numKmers);
			float prob;
			vector<float> zeroCountProb; zeroCountProb.resize(numKmers);		
			
			//read version
			string line2 = m->getline(inNum); m->gobble(inNum);
			
			while (inNum) {
				inNum >> zeroCountProb[count] >> num[count];  
				count++;
				m->gobble(inNum);
			}
			inNum.close();
		
			while(in) {
				in >> kmer;
				
				//set them all to zero value
				for (int i = 0; i < genusNodes.size(); i++) {
					wordGenusProb[kmer][i] = log(zeroCountProb[kmer] / (float) (genusTotals[i]+1));
				}
				
				//get probs for nonzero values
				for (int i = 0; i < num[kmer]; i++) {
					in >> name >> prob;
					wordGenusProb[kmer][name] = prob;
				}
				
				m->gobble(in);
			}
			in.close();
			
		#endif
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "readProbFile");
		exit(1);
	}
}
/**************************************************************************************************/
bool Bayesian::checkReleaseDate(ifstream& file1, ifstream& file2, ifstream& file3, ifstream& file4) {
	try {
		
		bool good = true;
		
		vector<string> lines;
		lines.push_back(m->getline(file1));  
		lines.push_back(m->getline(file2)); 
		lines.push_back(m->getline(file3)); 
		lines.push_back(m->getline(file4)); 

		//before we added this check
		if ((lines[0][0] != '#') || (lines[1][0] != '#') || (lines[2][0] != '#') || (lines[3][0] != '#')) {  good = false;  }
		else {
			//rip off #
			for (int i = 0; i < lines.size(); i++) { lines[i] = lines[i].substr(1);  }
			
			//get mothurs current version
			string version = m->getVersion();
			
			vector<string> versionVector;
			m->splitAtChar(version, versionVector, '.');
			
			//check each files version
			for (int i = 0; i < lines.size(); i++) { 
				vector<string> linesVector;
				m->splitAtChar(lines[i], linesVector, '.');
			
				if (versionVector.size() != linesVector.size()) { good = false; break; }
				else {
					for (int j = 0; j < versionVector.size(); j++) {
						int num1, num2;
						convert(versionVector[j], num1);
						convert(linesVector[j], num2);
						
						//if mothurs version is newer than this files version, then we want to remake it
						if (num1 > num2) {  good = false; break;  }
					}
				}
				
				if (!good) { break; }
			}
		}
		
		if (!good) {  file1.close(); file2.close(); file3.close(); file4.close();  }
		else { file1.seekg(0);  file2.seekg(0);  file3.seekg(0);  file4.seekg(0);  }
		
		return good;
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "checkReleaseDate");
		exit(1);
	}
}
/**************************************************************************************************/






