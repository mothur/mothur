#ifndef MOTHUR_FISHER
#define MOTHUR_FISHER

/*
 *  mothurfisher.h
 *  Mothur
 *
 *  Created by westcott on 7/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "mothurout.h"

class MothurFisher {
	
public:
	MothurFisher(){otuLabel = ""; m = MothurOut::getInstance(); }
	~MothurFisher(){}
	
	double fexact(double, double, double, double, string);
	
private:
	MothurOut* m;
	double sleft, sright, sless, slarg;
	double sn11,sn1_,sn_1,sn,sprob;
	double lngamm(double);
	double lnfact(double);
	double lnbico(double, double);
	double hyper_323(double, double, double, double);
	double myhyper(double);
	double hyper0(double, double, double, double);
	double exact(double, double, double, double);
    string otuLabel;
};


#endif

