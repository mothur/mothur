#ifndef FILE_UTILS_H
#define FILE_UTILS_H 

char *strpret(char *szString);

char *parseFromFile(const char *szFile, char *szParam);

int extractFloatFromFile(const char *szFile,char *szName,float *pFloat);

int extractDoubleFromFile(const char *szFile,char *szName,double *pDouble);

int extractIntFromFile(const char *szFile,char *szName,int *pInt);

int extractLongFromFile(const char *szFile,char *szName,long *pLong);

int extractUShortsFromFile(const char *szFile, char *szName, unsigned short* aShorts, int nShorts);

int extractIntMatrixFromFile(const char *szFile, char *szName, int** aanInts, int nRows, int nCols);

int extractFloatMatrixFromFile(const char *szFile, char *szName, float** aafFloats, int nRows, int nCols);

int extractDoubleMatrixFromFile(const char *szFile, char *szName, double** aadDoubles, int nRows, int nCols);

int extractDoublesFromFile(const char *szFile, char *szName, double* aDoubles, int nDoubles);

#endif


