#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

/* 2d integer matrix functions*/
int allocateMatrix2di(int ***paanInt, int m, int n);

int reallocateMatrix2di(int ***paanInt, int m1, int n1,int m2,int n2);

int callocateMatrix2di(int ***paanInt, int m, int n);

int recallocateMatrix2di(int ***paanInt,int m1,int n1,int m2,int n2);

int readMatrix2di(FILE *ifp, int **aanMatrix, int m, int n);

void writeMatrix2di(FILE *ofp, int **aanMatrix, int m, int n);

int destroyMatrix2di(int ***paanMatrix,int m);


/*2d float matrix functions*/
int allocateMatrix2df(float ***paafMatrix,int n, int m);

int reallocateMatrix2df(float ***paafMatrix,int m1, int n1,int m2,int n2);

int readMatrix2df(FILE *ifp, float **aafMatrix, int m, int n);

void writeMatrix2df(FILE *ofp, float  **aafMatrix, int m, int n);

int destroyMatrix2df(float ***paafMatrix, int m);


/*2d double matrix functions*/
int allocateMatrix2dd(double  ***paadMatrix,int n, int m);

int callocateMatrix2dd(double  ***paadMatrix,int n, int m);

int reallocateMatrix2dd(double ***paadMatrix,int m1, int n1,int m2,int n2);

int readMatrix2dd(FILE *ifp, double **aadMatrix, int m, int n);

void writeMatrix2dd(FILE *ofp, double **aadMatrix, int m, int n);

int destroyMatrix2dd(double ***paadMatrix, int m);


/*1d float matrix functions*/
int allocateMatrix1df(float **pafMatrix,int m);

int reallocateMatrix1df(float **pafMatrix,int m);

int destroyMatrix1df(float **pafMatrix);

/*1d double matrix functions*/
int allocateMatrix1dd(double **padMatrix,int m);

int callocateMatrix1dd(double **padMatrix,int m);

int reallocateMatrix1dd(double **padMatrix,int m);

/*1d integer matrix functions*/
int allocateMatrix1di(int **panMatrix,int m);

int callocateMatrix1di(int **panMatrix, int m);

int recallocateMatrix1di(int **panMatrix, int m1, int m2);

#endif





