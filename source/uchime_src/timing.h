#define TIMING 0
#ifndef timing_h
#define timing_h

#define BG_TIMING	0

#if !TIMING
#undef BG_TIMING
#define BG_TIMING	0
#endif

#if	UCHIMES
#undef TIMING
#define TIMING	0
#endif

#if TIMING

enum TIMER
	{
	TIMER_None,
#define T(x)	TIMER_##x,
#include "timers.h"
#undef T
	};

const unsigned TimerCount =
	1	// TIMER_None
#define T(x)	+1
#include "timers.h"
#undef T
	;

enum COUNTER
	{
#define C(x)	COUNTER_##x,
#include "counters.h"
#undef C
	};

enum ALLOCER
	{
#define A(x)	ALLOCER_##x,
#include "allocs.h"
#undef A
	};

const unsigned CounterCount =
#define C(x)	+1
#include "counters.h"
#undef C
	;

const unsigned AllocerCount =
#define A(x)	+1
#include "allocs.h"
#undef A
	;

#ifdef _MSC_VER

typedef unsigned __int64 TICKS;

#pragma warning(disable:4035)
inline TICKS GetClockTicks()
	{
	_asm
		{
		_emit	0x0f
		_emit	0x31
		}
	}

#else	// ifdef _MSC_VER

typedef uint64_t TICKS;
__inline__ uint64_t GetClockTicks()
	{
	uint32_t lo, hi;
	/* We cannot use "=A", since this would use %rax on x86_64 */
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return (uint64_t)hi << 32 | lo;
	}

#endif	// ifdef _MSC_VER

//void AddTicks(const string &Name, TICKS Ticks1, TICKS Ticks2);
//void AddBytes(const string &Name, double Bytes);
//#define SubBytes(Name, Bytes)	AddBytes(Name, -double(Bytes))

const char *TimerToStr(TIMER t);

extern TICKS g_BeginTicks[TimerCount];
extern double g_TotalTicks[TimerCount];
extern double g_TotalCounts[TimerCount];
extern double g_Counters[CounterCount];
extern unsigned g_AllocNewCount[AllocerCount];
extern unsigned g_AllocFreeCount[AllocerCount];
extern double g_AllocNewBytes[AllocerCount];
extern double g_AllocFreeBytes[AllocerCount];
extern double g_AllocNetBytes[AllocerCount];
extern double g_AllocPeakBytes[AllocerCount];
extern bool g_Timer2[TimerCount];
extern TIMER g_CurrTimer;
#if	BG_TIMING
extern TIMER g_BackgroundTimer;
#endif

#define MYALLOC(Type, N, Name)		(Type *) MyAlloc_((N)*sizeof(Type), ALLOCER_##Name, __FILE__, __LINE__)
#define MYFREE(Array, N, Name)		MyFree_(Array, N*sizeof(Array[0]), ALLOCER_##Name, __FILE__, __LINE__)

inline void *MyAlloc_(unsigned Bytes, unsigned a, const char *FileName, int Line)
	{
	++g_AllocNewCount[a];
	g_AllocNewBytes[a] += Bytes;
	g_AllocNetBytes[a] += Bytes;
	if (g_AllocNetBytes[a] > g_AllocPeakBytes[a])
		g_AllocPeakBytes[a] = g_AllocNetBytes[a];
	return mymalloc(Bytes);
	}

inline void MyFree_(void *p, unsigned Bytes, unsigned a, const char *FileName, int Line)
	{
	++g_AllocFreeCount[a];
	g_AllocFreeBytes[a] += Bytes;
	g_AllocNetBytes[a] -= Bytes;
	myfree2(p, Bytes);
	}

#if	BG_TIMING
inline void SetBackgroundTimer_(TIMER Timer)
	{
	TICKS Now = GetClockTicks();
	if (g_BeginTicks[g_BackgroundTimer] != 0)
		{
		++g_TotalCounts[g_BackgroundTimer];
		g_TotalTicks[g_BackgroundTimer] += double(Now - g_BeginTicks[g_BackgroundTimer]);
		}
	g_BackgroundTimer = Timer;
	g_BeginTicks[Timer] = Now;
	}
#else
#define SetBackgroundTimer_(Timer)	/* empty */
#endif

inline void StartTimer_(TIMER Timer)
	{
	if (g_CurrTimer != TIMER_None)
		Die("StartTimer(%s), curr=%s", TimerToStr(Timer), TimerToStr(g_CurrTimer));

	TICKS Now = GetClockTicks();
#if	BG_TIMING
	if (g_BeginTicks[g_BackgroundTimer] != 0)
		{
		++g_TotalCounts[g_BackgroundTimer];
		g_TotalTicks[g_BackgroundTimer] += double(Now - g_BeginTicks[g_BackgroundTimer]);
		}
#endif
	g_BeginTicks[Timer] = Now;
	g_CurrTimer = Timer;
	}

inline void PauseTimer_(TIMER Timer)
	{
	if (Timer != g_CurrTimer)
		Die("PauseTimer(%s), curr=%s", TimerToStr(Timer), TimerToStr(g_CurrTimer));

	TICKS Now = GetClockTicks();
	g_TotalTicks[Timer] += double(Now - g_BeginTicks[Timer]);
	g_BeginTicks[Timer] = Now;
	g_CurrTimer = TIMER_None;
	}

inline void EndTimer_(TIMER Timer)
	{
	if (Timer != g_CurrTimer)
		Die("EndTimer(%s), curr=%s", TimerToStr(Timer), TimerToStr(g_CurrTimer));

	TICKS Now = GetClockTicks();
#if	BG_TIMING
	g_BeginTicks[g_BackgroundTimer] = Now;
#endif
	g_TotalTicks[Timer] += double(Now - g_BeginTicks[Timer]);
	++g_TotalCounts[Timer];
	g_CurrTimer = TIMER_None;
	}

inline void StartTimer2_(TIMER Timer)
	{
	g_Timer2[Timer] = true;
	g_BeginTicks[Timer] = GetClockTicks();
	}

inline void EndTimer2_(TIMER Timer)
	{
	g_TotalTicks[Timer] += double(GetClockTicks() - g_BeginTicks[Timer]);
	++g_TotalCounts[Timer];
	}

#define AddCounter(x, N)	g_Counters[COUNTER_##x] += N
#define IncCounter(x)		++(g_Counters[COUNTER_##x])
#define StartTimer(x)		StartTimer_(TIMER_##x)
#define PauseTimer(x)		PauseTimer_(TIMER_##x)
#define EndTimer(x)			EndTimer_(TIMER_##x)
#define StartTimer2(x)		StartTimer2_(TIMER_##x)
#define EndTimer2(x)		EndTimer2_(TIMER_##x)

#if	BG_TIMING
#define SetBackgroundTimer(x)	SetBackgroundTimer_(TIMER_##x)
#else
#define SetBackgroundTimer(x)	/* empty */
#endif

#else	// if TIMING

#define AddCounter(x, N)	/* empty */
#define IncCounter(x)		/* empty */
#define StartTimer(x)		/* empty */
#define PauseTimer(x)		/* empty */
#define EndTimer(x)			/* empty */
#define StartTimer2(x)		/* empty */
#define PauseTimer2(x)		/* empty */
#define EndTimer2(x)		/* empty */
#define SetBackgroundTimer(x)	/* empty */
#define MYALLOC(Type, N, Name)		myalloc(Type, N)
#define MYFREE(Array, N, Name)		myfree(Array)

#endif	// if TIMING

void LogMemStats();
void LogTickStats();
void LogStats();
void LogAllocs();

#define AddBytes(x, n)	/* empty */
#define SubBytes(x, n)	/* empty */

#endif	// if timing_h
