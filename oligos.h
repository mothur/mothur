//
//  oligos.h
//  Mothur
//
//  Created by Sarah Westcott on 4/4/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#ifndef Mothur_oligos_h
#define Mothur_oligos_h


#include "mothurout.h"


/**************************************************************************************************/

class Oligos {
    
public:
	Oligos(string);
	Oligos();
    ~Oligos() {}
    
    int read(string);
    int read(string, bool); //read without reversing the paired barcodes, for make.contigs.
    bool hasPairedPrimers() { return hasPPrimers; }
    bool hasPairedBarcodes() { return hasPBarcodes; }
    
    //for processing with trimOligos class
    map<int, oligosPair> getPairedPrimers()                 { return pairedPrimers;     }
    map<int, oligosPair> getPairedBarcodes()                { return pairedBarcodes;    }
    map<string, int> getPrimers()                           { return primers;           }
    map<string, int> getBarcodes()                          { return barcodes;          }
    
    map<int, oligosPair> getReorientedPairedPrimers();
    map<int, oligosPair> getReorientedPairedBarcodes();
    map<string, int> getReorientedPrimers();
    map<string, int> getReorientedBarcodes();
    
    
    vector<string> getLinkers()                             { return linker;            }
    vector<string> getSpacers()                             { return spacer;            }
    vector<string> getReversePrimers()                      { return revPrimer;         }
    vector<string> getPrimerNames()                         { return primerNameVector;  }
    vector<string> getBarcodeNames()                        { return barcodeNameVector; }
    vector<string> getGroupNames()                          { return Groups;            }
        
    
    //for printing and other formatting uses
    vector<string> getBarcodes(string); //get barcodes for a group. For paired barcodes will return forward.reverse
    vector<string> getPrimers(string); //get primers for a group. For paired primers will return forward.reverse
    string getBarcodeName(int);
    string getPrimerName(int);
    string getGroupName(int, int);
    
		
private:
    
    set<string> uniqueNames; 
    vector<string> Groups;
	vector<string> revPrimer;
    map<string, vector<string> > Group2Barcode;
    map<string, vector<string> > Group2Primer;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
    map<string, int> primers;
	map<string, int> barcodes;
    vector<string>  linker;
    vector<string>  spacer;
	vector<string> primerNameVector;
	vector<string> barcodeNameVector;
    bool hasPPrimers, hasPBarcodes, pairedOligos, reversePairs;
    string oligosfile;
    int numBarcodes, numFPrimers;
    MothurOut* m;
    
    int indexPrimer;
    int indexBarcode;
    int indexPairedPrimer;
    int indexPairedBarcode;
    set<string> uniquePrimers;
    set<string> uniqueBarcodes;
    
    int readOligos();
    string reverseOligo(string);
};

/**************************************************************************************************/


#endif
