
#include "sparsematrix.hpp"
#include "listvector.hpp"


/***********************************************************************/

SparseMatrix::SparseMatrix() : numNodes(0), minsIndex(0), smallDist(1e6){  m = MothurOut::getInstance();  }

/***********************************************************************/

int SparseMatrix::getNNodes(){
	return numNodes; 
}

/***********************************************************************/

float SparseMatrix::getSmallDist(){
	return smallDist;
}

/***********************************************************************/

MatData SparseMatrix::rmCell(MatData data){
	try {
		if(data->vectorMap != NULL ){
			*(data->vectorMap) = NULL;
			data->vectorMap = NULL;
		}
		data = matrix.erase(data);		
		numNodes--;
		return(data);
	//  seems like i should be updating smallDist here, but the only time we remove cells is when
	//  clustering and the clustering algorithm updates smallDist
	}
	catch(exception& e) {
		m->errorOut(e, "SparseMatrix", "rmCell");
		exit(1);
	}
}

/***********************************************************************/

void SparseMatrix::addCell(PCell value){
	try {
		matrix.push_back(value);
		numNodes++;
		if(value.dist < smallDist){
			smallDist = value.dist;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SparseMatrix", "addCell");
		exit(1);
	}
}

/***********************************************************************/

void SparseMatrix::clear(){
	try {
		matrix.clear();
		mins.clear();
		numNodes = 0;
		minsIndex = 0;
		smallDist = 1e6;
	}
	catch(exception& e) {
		m->errorOut(e, "SparseMatrix", "clear");
		exit(1);
	}
}

/***********************************************************************/

MatData SparseMatrix::begin(){
	return matrix.begin();	
}

/***********************************************************************/

MatData SparseMatrix::end(){
	return matrix.end();	
}

/***********************************************************************/

void SparseMatrix::print(){
	try {
		int index = 0;
	
		cout << endl << "Index\tRow\tColumn\tDistance" << endl;
	
		for(MatData currentCell=matrix.begin();currentCell!=matrix.end();currentCell++){
			cout << index << '\t' << currentCell->row  << '\t' << currentCell->column << '\t' << currentCell->dist << endl;
			index++;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SparseMatrix", "print");
		exit(1);
	}
}

/***********************************************************************/

void SparseMatrix::print(ListVector* list){
	try {
		int index = 0;
	
		m->mothurOutEndLine(); m->mothurOut("Index\tRow\tColumn\tDistance"); m->mothurOutEndLine();
	
		for(MatData currentCell=matrix.begin();currentCell!=matrix.end();currentCell++){
			m->mothurOut(toString(index) + "\t" + toString(list->get(currentCell->row))  + "\t" + toString(list->get(currentCell->column)) + "\t" + toString(currentCell->dist)); m->mothurOutEndLine();
			index++;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SparseMatrix", "print");
		exit(1);
	}
}

/***********************************************************************/

PCell* SparseMatrix::getSmallestCell(){
	try {
	//	this is where I check to see if the next small distance has the correct distance
	//	if it doesn't then I remove the offending Cell -> should also be able to check for
	//	invalid iterator / pointer -- right???
	
		while(!mins.empty() && mins.back() == NULL){
			mins.pop_back();		
		}
	
	//	if the mins vector is empty go here...
		if(mins.empty()){		
			mins.clear();
	
			smallDist = begin()->dist;  //set the first candidate small distance
		
			for(MatData currentCell=begin();currentCell!=end();currentCell++){
			
				float dist = currentCell->dist;
			
				if(dist < smallDist){  //found a new smallest distance
					mins.clear();
					smallDist = dist;
					mins.push_back(&*currentCell);  //this is the address of the data in the list being pointed to by the MatData iterator
				}
				else if(dist == smallDist){  //if a subsequent distance is the same as mins distance add the new iterator to the mins vector
					mins.push_back(&*currentCell); //this is the address of the data in the list being pointed to by the MatData iterator
				}

			}
			m->mothurRandomShuffle(mins);  //randomize the order of the iterators in the mins vector

			for(int i=0;i<mins.size();i++){
				mins[i]->vectorMap = &mins[i];  //assign vectorMap to the address for the container
			}
			
		}
	
		smallCell = mins.back();	//make the smallestCell the last element of the vector

		mins.pop_back();			//remove the last element from the vector

		return smallCell;
	}
	catch(exception& e) {
		m->errorOut(e, "SparseMatrix", "getSmallestCell");
		exit(1);
	}
}

/***********************************************************************/

