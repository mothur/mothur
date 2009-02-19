#ifndef PROGRESS_H
#define PROGRESS_H

#include "mothur.h"

using namespace std;

class Progress {
	
public:
	Progress(string, int);
	void update(int);
	void finish();
	
private:
	int nTicks;
	int finalPos;	
};

#endif
