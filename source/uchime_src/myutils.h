#ifndef myutils_h
#define myutils_h

#define RCE_MALLOC	0

#include <stdio.h>
#include <sys/types.h>
#include <string>
#include <string.h>
#include <memory.h>
#include <vector>
#include <math.h>
#include <stdarg.h>
#include <cstdlib>
#include <climits>

#ifndef _MSC_VER
#include <inttypes.h>
#endif

using namespace std;

#ifdef _MSC_VER
#include <crtdbg.h>
#pragma warning(disable: 4996)	// deprecated functions
#define _CRT_SECURE_NO_DEPRECATE	1
#endif

#if defined(_DEBUG) && !defined(DEBUG)
#define DEBUG	1
#endif

#if defined(DEBUG) && !defined(_DEBUG)
#define _DEBUG	1
#endif

#ifndef NDEBUG
#define	DEBUG	1
#define	_DEBUG	1
#endif

typedef unsigned char byte;
typedef unsigned short uint16;
typedef unsigned uint32;
typedef int int32;
typedef double float32;
typedef signed char int8;
typedef unsigned char uint8;

#ifdef _MSC_VER

typedef __int64 int64;
typedef unsigned __int64 uint64;

#define INT64_PRINTF		"lld"
#define UINT64_PRINTF		"llu"

#define SIZE_T_PRINTF		"u"
#define OFF64_T_PRINTF		"lld"

#define INT64_PRINTFX		"llx"
#define UINT64_PRINTFX		"llx"

#define SIZE_T_PRINTFX		"x"
#define OFF64_T_PRINTFX		"llx"

#elif defined(__x86_64__)

typedef long int64;
typedef unsigned long uint64;

#define INT64_PRINTF		"ld"
#define UINT64_PRINTF		"lu"

#define SIZE_T_PRINTF		"lu"
#define OFF64_T_PRINTF		"ld"

#define INT64_PRINTFX		"lx"
#define UINT64_PRINTFX		"lx"

#define SIZE_T_PRINTFX		"lx"
#define OFF64_T_PRINTFX		"lx"

#else

typedef long long int64;
typedef unsigned long long uint64;

#define INT64_PRINTF		"lld"
#define UINT64_PRINTF		"llu"

#define SIZE_T_PRINTF		"u"
#define OFF64_T_PRINTF		"lld"

#define INT64_PRINTFX		"llx"
#define UINT64_PRINTFX		"llx"

#define SIZE_T_PRINTFX		"x"
#define OFF64_T_PRINTFX		"llx"
#endif

#define d64		INT64_PRINTF
#define	u64		UINT64_PRINTF
#define	x64		UINT64_PRINTFX

// const uint64 UINT64_MAX			= (~((uint64) 0));

void myassertfail(const char *Exp, const char *File, unsigned Line);
#undef  assert
#ifdef  NDEBUG
#define assert(exp)     ((void)0)
#define myassert(exp)     ((void)0)
#else
#define assert(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )
#define myassert(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )
#endif
#define asserta(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )

#define ureturn(x)	return (x)

#define NotUsed(v)	((void *) &v)

// pom=plus or minus, tof=true or false
static inline char pom(bool Plus)	{ return Plus ? '+' : '-'; }
static inline char tof(bool x)		{ return x ? 'T' : 'F';	}
static inline char yon(bool x)		{ return x ? 'Y' : 'N';	}
unsigned GetElapsedSecs();

#if	RCE_MALLOC

void *rce_malloc(unsigned bytes, const char *FileName, int Line);
void rce_free(void *p, const char *FileName, int LineNr);
void rce_chkmem();

void rce_dumpmem_(const char *FileName, int LineNr);
#define rce_dumpmem()		rce_dumpmem_(__FILE__, __LINE__)

void rce_assertvalidptr_(void *p, const char *FileName, int LineNr);
#define rce_assertvalidptr(p)	rce_assertvalidptr_(p, __FILE__, __LINE__)

void rce_dumpptr_(void *p, const char *FileName, int LineNr);
#define rce_dumpptr(p)	rce_dumpptr_(p, __FILE__, __LINE__)

#define mymalloc(n)		rce_malloc((n), __FILE__, __LINE__)
#define myfree(p)		rce_free(p, __FILE__, __LINE__)
#define myfree2(p,n)	rce_free(p, __FILE__, __LINE__)
#define myalloc(t, n)	(t *) rce_malloc((n)*sizeof(t), __FILE__, __LINE__)

#else // RCE_MALLOC
void *mymalloc(unsigned bytes);
void myfree2(void *p, unsigned Bytes);
void myfree(void *p);
#define rce_chkmem()	/* empty */
#define myalloc(t, n)	(t *) mymalloc((n)*sizeof(t))
#endif // RCE_MALLOC

#define SIZE(c)	unsigned((c).size())

bool myisatty(int fd);

#ifdef _MSC_VER
#define off_t	__int64
#endif

FILE *OpenStdioFile(const string &FileName);
FILE *CreateStdioFile(const string &FileName);
bool CanSetStdioFilePos(FILE *f);
void CloseStdioFile(FILE *f);
void SetStdioFilePos(FILE *f, off_t Pos);
void ReadStdioFile(FILE *f, off_t Pos, void *Buffer, unsigned Bytes);
void ReadStdioFile(FILE *f, void *Buffer, unsigned Bytes);
void WriteStdioFile(FILE *f, off_t Pos, const void *Buffer, unsigned Bytes);
void WriteStdioFile(FILE *f, const void *Buffer, unsigned Bytes);
bool ReadLineStdioFile(FILE *f, char *Line, unsigned Bytes);
bool ReadLineStdioFile(FILE *f, string &Line);
byte *ReadAllStdioFile(FILE *f, off_t &FileSize);
byte *ReadAllStdioFile(const string &FileName, off_t &FileSize);
void AppendStdioFileToFile(FILE *fFrom, FILE *fTo);
void FlushStdioFile(FILE *f);
bool StdioFileExists(const string &FileName);
off_t GetStdioFilePos(FILE *f);
off_t GetStdioFileSize(FILE *f);
void LogStdioFileState(FILE *f);
void RenameStdioFile(const string &FileNameFrom, const string &FileNameTo);
void DeleteStdioFile(const string &FileName);

void myvstrprintf(string &Str, const char *szFormat, va_list ArgList);
void myvstrprintf(string &Str, const char *szFormat, ...);

void SetLogFileName(const string &FileName);
void Log(const char *szFormat, ...);

void Die(const char *szFormat, ...);
void Warning(const char *szFormat, ...);

void ProgressStep(unsigned i, unsigned N, const char *Format, ...);
void Progress(const char *szFormat, ...);
void Progress(const string &Str);
void ProgressLog(const char *szFormat, ...);
void ProgressExit();

char *mystrsave(const char *s);

double GetPeakMemUseBytes();

// Are two floats equal to within epsilon?
const double epsilon = 0.01;
inline bool feq(double x, double y, double epsilon)
	{
	if (fabs(x) > 10000)
		epsilon = fabs(x)/10000;
	if (fabs(x - y) > epsilon)
		return false;
	return true;
	}

inline bool feq(double x, double y)
	{
	if (x < -1e6 && y < -1e6)
		return true;
	double e = epsilon;
	if (fabs(x) > 10000)
		e = fabs(x)/10000;
	if (fabs(x - y) > e)
		return false;
	return true;
	}

#define asserteq(x, y)	assert(feq(x, y))
#define assertaeq(x, y)	asserta(feq(x, y))

#define	zero(a, n)	memset(a, 0, n*sizeof(a[0]))

void InitRand();
unsigned randu32();
void Split(const string &Str, vector<string> &Fields, char Sep = 0);
double Pct(double x, double y);
double GetMemUseBytes();
const char *MemBytesToStr(double Bytes);
const char *IntToStr(unsigned i);
const char *FloatToStr(double d);
const char *SecsToStr(double Secs);
void Logu(unsigned u, unsigned w, unsigned prefixspaces = 2);
void Logf(float x, unsigned w, unsigned prefixspaces = 2);
const char *SecsToHHMMSS(int Secs);

void MyCmdLine(int argc, char **argv);
void CmdLineErr(const char *Format, ...);
void Help();
void GetCmdLine(string &s);

#define FLAG_OPT(LongName)						extern bool opt_##LongName; extern bool optset_##LongName;
#define TOG_OPT(LongName, Default)				extern bool opt_##LongName; extern bool optset_##LongName;
#define INT_OPT(LongName, Default, Min, Max)	extern int opt_##LongName; extern bool optset_##LongName;
#define UNS_OPT(LongName, Default, Min, Max)	extern unsigned opt_##LongName; extern bool optset_##LongName;
#define FLT_OPT(LongName, Default, Min, Max)	extern double opt_##LongName; extern bool optset_##LongName;
#define STR_OPT(LongName, Default)				extern string opt_##LongName; extern bool optset_##LongName;
#define ENUM_OPT(LongName, Default, Values)		extern int opt_##LongName; extern bool optset_##LongName;
#include "myopts.h"
#undef FLAG_OPT
#undef TOG_OPT
#undef INT_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
#undef ENUM_OPT

extern const char *SVN_VERSION;
extern const char *SVN_MODS;
extern bool opt_quiet;
extern bool opt_version;
extern FILE *g_fLog;

#endif	// myutils_h
