//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#ifndef evalue_h
#define evalue_h

#include <float.h>

void SetKarlin(double GappedLambda, double UngappedLambda,
  double GappedK, double UngappedK, double DBLength);\

double GetKarlinDBLength();
void SetKarlinDBLength(double DBLength);
void LogKarlin();
void SetKarlinAmino(double DBLength);
void SetKarlinNucleo(double DBLength);
void SetKarlin(double DBLength, bool Nucleo);
double ComputeBitScoreGapped(double Score);
double ComputeBitScoreUngapped(double Score);
double ComputeEvalueGapped(double Score, unsigned QueryLength);
double ComputeEvalueUngapped(double Score, unsigned QueryLength);
double ComputeMinScoreGivenEvalueAGapped(double Evalue, unsigned Area);
double ComputeMinScoreGivenEvalueAUngapped(double Evalue, unsigned Area);
double ComputeMinScoreGivenEvalueQGapped(double Evalue, unsigned QueryLength);
double ComputeMinScoreGivenEvalueQUngapped(double Evalue, unsigned QueryLength);
double ComputeEvalueGappedFromBitScore(double BitScore, unsigned QueryLength);

#endif // evalue_h
