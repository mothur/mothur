//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <signal.h>
#include <float.h>

#ifdef _MSC_VER
#include <crtdbg.h>
#include <process.h>
#include <windows.h>
#include <psapi.h>
#include <io.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#endif

#include "myutils.h"

const char *SVN_VERSION =
#include "svnversion.h"
;

#define	TEST_UTILS			0

using namespace std;

const unsigned MY_IO_BUFSIZ = 32000;
const unsigned MAX_FORMATTED_STRING_LENGTH = 64000;

static char *g_IOBuffers[256];
static time_t g_StartTime = time(0);
static vector<string> g_Argv;
static double g_PeakMemUseBytes;

#if	TEST_UTILS
void TestUtils()
	{
	const int C = 100000000;
	for (int i = 0; i < C; ++i)
		ProgressStep(i, C, "something or other");

	Progress("\n");
	Progress("Longer message\r");
	Sleep(1000);
	Progress("Short\r");
	Sleep(1000);
	Progress("And longer again\r");
	Sleep(1000);
	Progress("Shrt\n");
	Sleep(1000);
	const unsigned N = 10;
	unsigned M = 10;
	for (unsigned i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Allocating 1MB blocks");
		for (unsigned j = 0; j < M; ++j)
			{
			ProgressStep(j, M, "Inner loop"); 
			malloc(100000);
			Sleep(500);
			}
		}
	}
#endif // TEST_UTILS

static void AllocBuffer(FILE *f)
	{
	int fd = fileno(f);
	if (fd < 0 || fd >= 256)
		return;
	if (g_IOBuffers[fd] == 0)
		g_IOBuffers[fd] = myalloc(char, MY_IO_BUFSIZ);
	setvbuf(f, g_IOBuffers[fd], _IOFBF, MY_IO_BUFSIZ);
	}

static void FreeBuffer(FILE *f)
	{
	int fd = fileno(f);
	if (fd < 0 || fd >= 256)
		return;
	if (g_IOBuffers[fd] == 0)
		return;
	myfree(g_IOBuffers[fd]);
	g_IOBuffers[fd] = 0;
	}

unsigned GetElapsedSecs()
	{
	return (unsigned) (time(0) - g_StartTime);
	}

static unsigned g_NewCalls;
static unsigned g_FreeCalls;
static double g_InitialMemUseBytes;
static double g_TotalAllocBytes;
static double g_TotalFreeBytes;
static double g_NetBytes;
static double g_MaxNetBytes;

void LogAllocStats()
	{
	Log("\n");
	Log("       Allocs  %u\n", g_NewCalls);
	Log("        Frees  %u\n", g_FreeCalls);
	Log("Initial alloc  %s\n", MemBytesToStr(g_InitialMemUseBytes));
	Log("  Total alloc  %s\n", MemBytesToStr(g_TotalAllocBytes));
	Log("   Total free  %s\n", MemBytesToStr(g_TotalFreeBytes));
	Log("    Net bytes  %s\n", MemBytesToStr(g_NetBytes));
	Log("Max net bytes  %s\n", MemBytesToStr(g_MaxNetBytes));
	Log("   Peak total  %s\n", MemBytesToStr(g_MaxNetBytes + g_InitialMemUseBytes));
	}

bool StdioFileExists(const string &FileName)
	{
	struct stat SD;
	int i = stat(FileName.c_str(), &SD);
	return i == 0;
	}

void myassertfail(const char *Exp, const char *File, unsigned Line)
	{
	Die("%s(%u) assert failed: %s", File, Line, Exp);
	}

bool myisatty(int fd)
	{
	return isatty(fd) != 0;
	}

#ifdef _MSC_VER
#include <io.h>
int fseeko(FILE *stream, off_t offset, int whence)
	{
	off_t FilePos = _fseeki64(stream, offset, whence);
	return (FilePos == -1L) ? -1 : 0;
	}
#define ftello(fm) (off_t) _ftelli64(fm)
#endif

void LogStdioFileState(FILE *f)
	{
	unsigned long tellpos = (unsigned long) ftello(f);
	long fseek_pos = fseek(f, 0, SEEK_CUR);
	int fd = fileno(f);
	Log("FILE *     %p\n", f);
	Log("fileno     %d\n", fd);
	Log("feof       %d\n", feof(f));
	Log("ferror     %d\n", ferror(f));
	Log("ftell      %ld\n", tellpos);
	Log("fseek      %ld\n", fseek_pos);
#if	!defined(_GNU_SOURCE) && !defined(__APPLE_CC__)
	fpos_t fpos;
	int fgetpos_retval = fgetpos(f, &fpos);
	Log("fpos       %ld (retval %d)\n", (long) fpos, fgetpos_retval);
//	Log("eof        %d\n", _eof(fd));
#endif
#ifdef _MSC_VER
	__int64 pos64 = _ftelli64(f);
	Log("_ftelli64  %lld\n", pos64);
#endif
	}

FILE *OpenStdioFile(const string &FileName)
	{
	const char *Mode = "rb";
	FILE *f = fopen(FileName.c_str(), Mode);
	if (f == 0)
		{
		if (errno == EFBIG)
			{
			if (sizeof(off_t) == 4)
				Die("File too big, off_t is 32 bits, recompile needed");
			else
				Die("Cannot open '%s', file too big (off_t=%u bits)",
				  FileName.c_str(), sizeof(off_t)*8);
			}
		Die("Cannot open %s, errno=%d %s",
		  FileName.c_str(), errno, strerror(errno));
		}
	AllocBuffer(f);
	return f;
	}

FILE *CreateStdioFile(const string &FileName)
	{
	FILE *f = fopen(FileName.c_str(), "wb+");
	if (0 == f)
		Die("Cannot create %s, errno=%d %s",
		  FileName.c_str(), errno, strerror(errno));
	AllocBuffer(f);
	return f;
	}

void SetStdioFilePos(FILE *f, off_t Pos)
	{
	if (0 == f)
		Die("SetStdioFilePos failed, f=NULL");
	int Ok = fseeko(f, Pos, SEEK_SET);
	off_t NewPos = ftello(f);
	if (Ok != 0 || Pos != NewPos)
		{
		LogStdioFileState(f);
		Die("SetStdioFilePos(%d) failed, Ok=%d NewPos=%d",
		  (int) Pos, Ok, (int) NewPos);
		}
	}

void ReadStdioFile(FILE *f, off_t Pos, void *Buffer, unsigned Bytes)
	{
	if (0 == f)
		Die("ReadStdioFile failed, f=NULL");
	SetStdioFilePos(f, Pos);
	unsigned BytesRead = fread(Buffer, 1, Bytes, f);
	if (BytesRead != Bytes)
		{
		LogStdioFileState(f);
		Die("ReadStdioFile failed, attempted %d bytes, read %d bytes, errno=%d",
		  (int) Bytes, (int) BytesRead, errno);
		}
	}

void ReadStdioFile(FILE *f, void *Buffer, unsigned Bytes)
	{
	if (0 == f)
		Die("ReadStdioFile failed, f=NULL");
	unsigned BytesRead = fread(Buffer, 1, Bytes, f);
	if (BytesRead != Bytes)
		{
		LogStdioFileState(f);
		Die("ReadStdioFile failed, attempted %d bytes, read %d bytes, errno=%d",
		  (int) Bytes, (int) BytesRead, errno);
		}
	}

// Return values from functions like lseek, ftell, fgetpos are
// "undefined" for files that cannot seek. Attempt to detect
// whether a file can seek by checking for error returns.
bool CanSetStdioFilePos(FILE *f)
	{
// Common special cases
	if (f == stdin || f == stdout || f == stderr)
		return false;

	fpos_t CurrPos;
	int ok1 = fgetpos(f, &CurrPos);
	if (ok1 < 0)
		return false;
	int ok2 = fseek(f, 0, SEEK_END);
	if (ok2 < 0)
		return false;
	fpos_t EndPos;
	int ok3 = fgetpos(f, &EndPos);
	int ok4 = fsetpos(f, &CurrPos);
	if (!ok3 || !ok4)
		return false;
	return true;
	}

byte *ReadAllStdioFile(FILE *f, unsigned &FileSize)
	{
	const unsigned BUFF_SIZE = 1024*1024;

	if (CanSetStdioFilePos(f))
		{
		off_t Pos = GetStdioFilePos(f);
		off_t FileSize = GetStdioFileSize(f);
		if (FileSize > UINT_MAX)
			Die("ReadAllStdioFile: file size > UINT_MAX");
		SetStdioFilePos(f, 0);
		byte *Buffer = myalloc(byte, unsigned(FileSize));
		ReadStdioFile(f, Buffer, unsigned(FileSize));
		SetStdioFilePos(f, Pos);
		FileSize = unsigned(FileSize);
		return Buffer;
		}

// Can't seek, read one buffer at a time.
	FileSize = 0;

// Just to initialize so that first call to realloc works.
	byte *Buffer = (byte *) malloc(4);
	if (Buffer == 0)
		Die("ReadAllStdioFile, out of memory");
	for (;;)
		{
		Buffer = (byte *) realloc(Buffer, FileSize + BUFF_SIZE);
		unsigned BytesRead = fread(Buffer + FileSize, 1, BUFF_SIZE, f);
		FileSize += BytesRead;
		if (BytesRead < BUFF_SIZE)
			{
			Buffer = (byte *) realloc(Buffer, FileSize);
			return Buffer;
			}
		}
	}

byte *ReadAllStdioFile(const std::string &FileName, off_t &FileSize)
	{
#if	WIN32
	FILE *f = OpenStdioFile(FileName);
	FileSize = GetStdioFileSize(f);
	CloseStdioFile(f);

	HANDLE h = CreateFile(FileName.c_str(), GENERIC_READ, FILE_SHARE_READ,
	  NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (h == INVALID_HANDLE_VALUE)
		Die("ReadAllStdioFile:Open(%s) failed", FileName.c_str());

	unsigned uFileSize = (unsigned) FileSize;
	if ((off_t) uFileSize != FileSize)
		Die("File too big (%.1f Gb): %s", double(FileSize)/1e9, FileName.c_str());

	byte *Buffer = myalloc(byte, uFileSize);
	DWORD BytesRead;
	ReadFile(h, Buffer, uFileSize, &BytesRead, NULL);
	if (FileSize != BytesRead)
		Die("ReadAllStdioFile:Error reading %s, attempted %u got %u",
		  FileName.c_str(), FileSize, (unsigned) BytesRead);

	CloseHandle(h);
	return Buffer;
#else
	int h = open(FileName.c_str(), O_RDONLY);
	if (h < 0)
		Die("ReadAllStdioFile:Cannot open %s", FileName.c_str());
	FileSize = lseek(h, 0, SEEK_END);
	if (FileSize == (off_t) (-1))
		Die("ReadAllStdioFile:Error seeking %s", FileName.c_str());
	// byte *Buffer = myalloc<byte>(FileSize);
	size_t stBytes = (size_t) FileSize;
	if ((off_t) stBytes != FileSize)
		Die("ReadAllStdioFile: off_t overflow");
	byte *Buffer = (byte *) malloc(stBytes);
	if (Buffer == 0)
		Die("ReadAllStdioFile: failed to allocate %s", MemBytesToStr(stBytes));
	lseek(h, 0, SEEK_SET);
	size_t n = read(h, Buffer, stBytes);
	if (n != FileSize)
		Die("ReadAllStdioFile, Error reading %s, attempted %g got %g",
		  FileName.c_str(), (double) FileSize, (double) n);
	close(h);
	return Buffer;
#endif
	}

void WriteStdioFile(FILE *f, off_t Pos, const void *Buffer, unsigned Bytes)
	{
	if (0 == f)
		Die("WriteStdioFile failed, f=NULL");
	SetStdioFilePos(f, Pos);
	unsigned BytesWritten = fwrite(Buffer, 1, Bytes, f);
	if (BytesWritten != Bytes)
		{
		LogStdioFileState(f);
		Die("WriteStdioFile failed, attempted %d bytes, wrote %d bytes, errno=%d",
		  (int) Bytes, (int) BytesWritten, errno);
		}
	}

void WriteStdioFile(FILE *f, const void *Buffer, unsigned Bytes)
	{
	if (0 == f)
		Die("WriteStdioFile failed, f=NULL");
	unsigned BytesWritten = fwrite(Buffer, 1, Bytes, f);
	if (BytesWritten != Bytes)
		{
		LogStdioFileState(f);
		Die("WriteStdioFile failed, attempted %d bytes, wrote %d bytes, errno=%d",
		  (int) Bytes, (int) BytesWritten, errno);
		}
	}

// Return false on EOF, true if line successfully read.
bool ReadLineStdioFile(FILE *f, char *Line, unsigned Bytes)
	{
	if (feof(f))
		return false;
	if ((int) Bytes < 0)
		Die("ReadLineStdioFile: Bytes < 0");
	char *RetVal = fgets(Line, (int) Bytes, f);
	if (NULL == RetVal)
		{
		if (feof(f))
			return false;
		if (ferror(f))
			Die("ReadLineStdioFile: errno=%d", errno);
		Die("ReadLineStdioFile: fgets=0, feof=0, ferror=0");
		}

	if (RetVal != Line)
		Die("ReadLineStdioFile: fgets != Buffer");
	unsigned n = strlen(Line);
	if (n < 1 || Line[n-1] != '\n')
		Die("ReadLineStdioFile: line too long or missing end-of-line");
	if (n > 0 && (Line[n-1] == '\r' || Line[n-1] == '\n'))
		Line[n-1] = 0;
	if (n > 1 && (Line[n-2] == '\r' || Line[n-2] == '\n'))
		Line[n-2] = 0;
	return true;
	}

// Return false on EOF, true if line successfully read.
bool ReadLineStdioFile(FILE *f, string &Line)
	{
	Line.clear();
	for (;;)
		{
		int c = fgetc(f);
		if (c == -1)
			{
			if (feof(f))
				{
				if (!Line.empty())
					return true;
				return false;
				}
			Die("ReadLineStdioFile, errno=%d", errno);
			}
		if (c == '\r')
			continue;
		if (c == '\n')
			return true;
		Line.push_back((char) c);
		}
	}

// Copies all of fFrom regardless of current
// file position, appends to fTo.
void AppendStdioFileToFile(FILE *fFrom, FILE *fTo)
	{
	off_t SavedFromPos = GetStdioFilePos(fFrom);
	off_t FileSize = GetStdioFileSize(fFrom);
	const off_t BUFF_SIZE = 1024*1024;
	char *Buffer = myalloc(char, BUFF_SIZE);
	SetStdioFilePos(fFrom, 0);
	off_t BytesRemaining = FileSize;
	while (BytesRemaining > 0)
		{
		off_t BytesToRead = BytesRemaining;
		if (BytesToRead > BUFF_SIZE)
			BytesToRead = BUFF_SIZE;
		ReadStdioFile(fFrom, Buffer, (unsigned) BytesToRead);
		WriteStdioFile(fTo, Buffer, (unsigned) BytesToRead);
		BytesRemaining -= BytesToRead;
		}
	SetStdioFilePos(fFrom, SavedFromPos);
	}

void RenameStdioFile(const string &FileNameFrom, const string &FileNameTo)
	{
	int Ok = rename(FileNameFrom.c_str(), FileNameTo.c_str());
	if (Ok != 0)
		Die("RenameStdioFile(%s,%s) failed, errno=%d %s",
		  FileNameFrom.c_str(), FileNameTo.c_str(), errno, strerror(errno));
	}

void FlushStdioFile(FILE *f)
	{
	int Ok = fflush(f);
	if (Ok != 0)
		Die("fflush(%p)=%d,", f, Ok);
	}

void CloseStdioFile(FILE *f)
	{
	if (f == 0)
		return;
	int Ok = fclose(f);
	if (Ok != 0)
		Die("fclose(%p)=%d", f, Ok);
	FreeBuffer(f);
	}

off_t GetStdioFilePos(FILE *f)
	{
	off_t FilePos = ftello(f);
	if (FilePos < 0)
		Die("ftello=%d", (int) FilePos);
	return FilePos;
	}

off_t GetStdioFileSize(FILE *f)
	{
	off_t CurrentPos = GetStdioFilePos(f);
	int Ok = fseeko(f, 0, SEEK_END);
	if (Ok < 0)
		Die("fseek in GetFileSize");

	off_t Length = ftello(f);
	if (Length < 0)
		Die("ftello in GetFileSize");
	SetStdioFilePos(f, CurrentPos);
	return Length;
	}

void DeleteStdioFile(const string &FileName)
	{
	int Ok = remove(FileName.c_str());
	if (Ok != 0)
		Die("remove(%s) failed, errno=%d %s", FileName.c_str(), errno, strerror(errno));
	}

void myvstrprintf(string &Str, const char *Format, va_list ArgList)
	{
	static char szStr[MAX_FORMATTED_STRING_LENGTH];
	vsnprintf(szStr, MAX_FORMATTED_STRING_LENGTH-1, Format, ArgList);
	szStr[MAX_FORMATTED_STRING_LENGTH - 1] = '\0';
	Str.assign(szStr);
	}

void myvstrprintf(string &Str, const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);
	}

FILE *g_fLog = 0;

void SetLogFileName(const string &FileName)
	{
	if (g_fLog != 0)
		CloseStdioFile(g_fLog);
	g_fLog = 0;
	if (FileName.empty())
		return;
	g_fLog = CreateStdioFile(FileName);
	}

void Log(const char *Format, ...)
	{
	if (g_fLog == 0)
		return;

	static bool InLog = false;
	if (InLog)
		return;

	InLog = true;
	va_list ArgList;
	va_start(ArgList, Format);
	vfprintf(g_fLog, Format, ArgList);
	va_end(ArgList);
	fflush(g_fLog);
	InLog = false;
	}

void Die(const char *Format, ...)
	{
	static bool InDie = false;
	if (InDie)
		exit(1);
	InDie = true;
	string Msg;

	if (g_fLog != 0)
		setbuf(g_fLog, 0);
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Msg, Format, ArgList);
	va_end(ArgList);

	fprintf(stderr, "\n\n");
	Log("\n");
	time_t t = time(0);
	Log("%s", asctime(localtime(&t)));
	for (unsigned i = 0; i < g_Argv.size(); i++)
		{
		fprintf(stderr, (i == 0) ? "%s" : " %s", g_Argv[i].c_str());
		Log((i == 0) ? "%s" : " %s", g_Argv[i].c_str());
		}
	fprintf(stderr, "\n");
	Log("\n");

	time_t CurrentTime = time(0);
	unsigned ElapsedSeconds = unsigned(CurrentTime - g_StartTime);
	const char *sstr = SecsToStr(ElapsedSeconds);
	Log("Elapsed time: %s\n", sstr);

	const char *szStr = Msg.c_str();
	fprintf(stderr, "\n---Fatal error---\n%s\n", szStr);
	Log("\n---Fatal error---\n%s\n", szStr);

#ifdef _MSC_VER
	if (IsDebuggerPresent())
 		__debugbreak();
	_CrtSetDbgFlag(0);
#endif

	exit(1);
	}

void Warning(const char *Format, ...)
	{
	string Msg;

	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Msg, Format, ArgList);
	va_end(ArgList);

	const char *szStr = Msg.c_str();

	fprintf(stderr, "\nWARNING: %s\n", szStr);
	if (g_fLog != stdout)
		{
		Log("\nWARNING: %s\n", szStr);
		fflush(g_fLog);
		}
	}

#ifdef _MSC_VER
double GetMemUseBytes()
	{
	HANDLE hProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS PMC;
	BOOL bOk = GetProcessMemoryInfo(hProc, &PMC, sizeof(PMC));
	if (!bOk)
		return 1000000;
	double Bytes = (double) PMC.WorkingSetSize;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
	}
#elif	linux || __linux__
double GetMemUseBytes()
	{
	static char statm[64];
	static int PageSize = 1;
	if (0 == statm[0])
		{
		PageSize = sysconf(_SC_PAGESIZE);
		pid_t pid = getpid();
		sprintf(statm, "/proc/%d/statm", (int) pid);
		}

	int fd = open(statm, O_RDONLY);
	if (-1 == fd)
		return 1000000;
	char Buffer[64];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		return 1000000;

	Buffer[n] = 0;
	double Pages = atof(Buffer);

	double Bytes = Pages*PageSize;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
	}
#elif defined(__MACH__)
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/socket.h>
#include <sys/gmon.h>
#include <mach/vm_param.h>
#include <netinet/in.h>
#include <netinet/icmp6.h>
#include <sys/vmmeter.h>
#include <sys/proc.h>
#include <mach/task_info.h>
#include <mach/task.h>
#include <mach/mach_init.h>
#include <mach/vm_statistics.h>

#define DEFAULT_MEM_USE	100000000.0

double GetMemUseBytes()
	{
	task_t mytask = mach_task_self();
	struct task_basic_info ti;
	memset((void *) &ti, 0, sizeof(ti));
	mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;
	kern_return_t ok = task_info(mytask, TASK_BASIC_INFO, (task_info_t) &ti, &count);
	if (ok == KERN_INVALID_ARGUMENT)
		return DEFAULT_MEM_USE;

	if (ok != KERN_SUCCESS)
		return DEFAULT_MEM_USE;

	double Bytes = (double ) ti.resident_size;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
	}
#else
double GetMemUseBytes()
	{
	return 0;
	}
#endif

double GetPeakMemUseBytes()
	{
	return g_PeakMemUseBytes;
	}

const char *SecsToHHMMSS(int Secs)
	{
	int HH = Secs/3600;
	int MM = (Secs - HH*3600)/60;
	int SS = Secs%60;
	static char Str[16];
	if (HH == 0)
		sprintf(Str, "%02d:%02d", MM, SS);
	else
		sprintf(Str, "%02d:%02d:%02d", HH, MM, SS);
	return Str;
	}

const char *SecsToStr(double Secs)
	{
	if (Secs >= 10.0)
		return SecsToHHMMSS((int) Secs);

	static char Str[16];
	if (Secs < 1e-6)
		sprintf(Str, "%.2gs", Secs);
	else if (Secs < 1e-3)
		sprintf(Str, "%.2fms", Secs*1e3);
	else
		sprintf(Str, "%.3fs", Secs);
	return Str;
	}

const char *MemBytesToStr(double Bytes)
	{
	static char Str[32];

	if (Bytes < 1e6)
		sprintf(Str, "%.1fkb", Bytes/1e3);
	else if (Bytes < 10e6)
		sprintf(Str, "%.1fMb", Bytes/1e6);
	else if (Bytes < 1e9)
		sprintf(Str, "%.0fMb", Bytes/1e6);
	else if (Bytes < 10e9)
		sprintf(Str, "%.1fGb", Bytes/1e9);
	else if (Bytes < 100e9)
		sprintf(Str, "%.0fGb", Bytes/1e9);
	else
		sprintf(Str, "%.3gb", Bytes);
	return Str;
	}

const char *IntToStr(unsigned i)
	{
	static char Str[32];

	double d = (double) i;
	if (i < 10000)
		sprintf(Str, "%u", i);
	else if (i < 1e6)
		sprintf(Str, "%.1fk", d/1e3);
	else if (i < 10e6)
		sprintf(Str, "%.1fM", d/1e6);
	else if (i < 1e9)
		sprintf(Str, "%.0fM", d/1e6);
	else if (i < 10e9)
		sprintf(Str, "%.1fG", d/1e9);
	else if (i < 100e9)
		sprintf(Str, "%.0fG", d/1e9);
	else
		sprintf(Str, "%.3g", d);
	return Str;
	}

const char *FloatToStr(double d)
	{
	static char Str[32];

	double a = fabs(d);
	if (a < 0.01)
		sprintf(Str, "%.3g", a);
	else if (a >= 0.01 && a < 1)
		sprintf(Str, "%.3f", a);
	else if (a <= 10 && a >= 1)
		{
		double intpart;
		if (modf(a, &intpart) < 0.05)
			sprintf(Str, "%.0f", d);
		else
			sprintf(Str, "%.1f", d);
		}
	else if (a > 10 && a < 10000)
		sprintf(Str, "%.0f", d);
	else if (a < 1e6)
		sprintf(Str, "%.1fk", d/1e3);
	else if (a < 10e6)
		sprintf(Str, "%.1fM", d/1e6);
	else if (a < 1e9)
		sprintf(Str, "%.0fM", d/1e6);
	else if (a < 10e9)
		sprintf(Str, "%.1fG", d/1e9);
	else if (a < 100e9)
		sprintf(Str, "%.0fG", d/1e9);
	else
		sprintf(Str, "%.3g", d);
	return Str;
	}

bool opt_quiet = false;
bool opt_version = false;
bool opt_logopts = false;
bool opt_compilerinfo = false;
bool opt_help = false;
string opt_log = "";

bool optset_quiet = false;
bool optset_version = false;
bool optset_logopts = false;
bool optset_compilerinfo = false;
bool optset_help = false;
bool optset_log = false;

static string g_CurrentProgressLine;
static string g_ProgressDesc;
static unsigned g_ProgressIndex;
static unsigned g_ProgressCount;

static unsigned g_CurrProgressLineLength;
static unsigned g_LastProgressLineLength;
static unsigned g_CountsInterval;
static unsigned g_StepCalls;
static time_t g_TimeLastOutputStep;

static string &GetProgressPrefixStr(string &s)
	{
	double Bytes = GetMemUseBytes();
	unsigned Secs = GetElapsedSecs();
	s = string(SecsToHHMMSS(Secs));
	if (Bytes > 0)
		{
		s.push_back(' ');
		char Str[32];
		sprintf(Str, "%5.5s", MemBytesToStr(Bytes));
		s += string(Str);
		}
	s.push_back(' ');
	return s;
	}

void ProgressLog(const char *Format, ...)
	{
	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

	Log("%s", Str.c_str());
	Progress("%s", Str.c_str());
	}

void Progress(const char *Format, ...)
	{
	if (opt_quiet)
		return;

	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

#if	0
	Log("Progress(");
	for (unsigned i = 0; i < Str.size(); ++i)
		{
		char c = Str[i];
		if (c == '\r')
			Log("\\r");
		else if (c == '\n')
			Log("\\n");
		else
			Log("%c", c);
		}
	Log(")\n");
#endif //0

	for (unsigned i = 0; i < Str.size(); ++i)
		{
		if (g_CurrProgressLineLength == 0)
			{
			string s;
			GetProgressPrefixStr(s);
			for (unsigned j = 0; j < s.size(); ++j)
				{
				fputc(s[j], stderr);
				++g_CurrProgressLineLength;
				}
			}

		char c = Str[i];
		if (c == '\n' || c == '\r')
			{
			for (unsigned j = g_CurrProgressLineLength; j < g_LastProgressLineLength; ++j)
				fputc(' ', stderr);
			if (c == '\n')
				g_LastProgressLineLength = 0;
			else
				g_LastProgressLineLength = g_CurrProgressLineLength;
			g_CurrProgressLineLength = 0;
			fputc(c, stderr);
			}
		else
			{
			fputc(c, stderr);
			++g_CurrProgressLineLength;
			}
		}
	}

void ProgressExit()
	{
	time_t Now = time(0);
	struct tm *t = localtime(&Now);
	const char *s = asctime(t);
	unsigned Secs = GetElapsedSecs();

	Log("\n");
	Log("Finished %s", s); // there is a newline in s
	Log("Elapsed time %s\n", SecsToHHMMSS((int) Secs));
	Log("Max memory %s\n", MemBytesToStr(g_PeakMemUseBytes));
#if	WIN32 && DEBUG
// Skip exit(), which can be very slow in DEBUG build
// VERY DANGEROUS practice, because it skips global destructors.
// But if you know the rules, you can break 'em, right?
	ExitProcess(0);
#endif
	}

const char *PctStr(double x, double y)
	{
	if (y == 0)
		{
		if (x == 0)
			return "100%";
		else
			return "inf%";
		}
	static char Str[16];
	double p = x*100.0/y;
	sprintf(Str, "%5.1f%%", p);
	return Str;
	}

string &GetProgressLevelStr(string &s)
	{
	unsigned Index = g_ProgressIndex;
	unsigned Count = g_ProgressCount;
	if (Count == UINT_MAX)
		{
		if (Index == UINT_MAX)
			s = "100%";
		else
			{
			char Tmp[16];
			sprintf(Tmp, "%u", Index); 
			s = Tmp;
			}
		}
	else
		s = string(PctStr(Index+1, Count));
	s += string(" ") + g_ProgressDesc;
	return s;
	}

void ProgressStep(unsigned i, unsigned N, const char *Format, ...)
	{
	if (opt_quiet)
		return;

	if (i == 0)
		{
		string Str;
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(Str, Format, ArgList);
		va_end(ArgList);
		g_ProgressDesc = Str;
		g_ProgressIndex = 0;
		g_ProgressCount = N;
		g_CountsInterval = 1;
		g_StepCalls = 0;
		g_TimeLastOutputStep = 0;
		if (g_CurrProgressLineLength > 0)
			Progress("\n");
		}

	if (i >= N && i != UINT_MAX)
		Die("ProgressStep(%u,%u)", i, N);
	bool IsLastStep = (i == UINT_MAX || i + 1 == N);
	if (!IsLastStep)
		{
		++g_StepCalls;
		if (g_StepCalls%g_CountsInterval != 0)
			return;

		time_t Now = time(0);
		if (Now == g_TimeLastOutputStep)
			{
			if (g_CountsInterval < 128)
				g_CountsInterval = (g_CountsInterval*3)/2;
			else
				g_CountsInterval += 64;
			return;
			}
		else
			{
			time_t Secs = Now - g_TimeLastOutputStep;
			if (Secs > 1)
				g_CountsInterval = unsigned(g_CountsInterval/(Secs*8));
			}

		if (g_CountsInterval < 1)
			g_CountsInterval = 1;

		g_TimeLastOutputStep = Now;
		}

	g_ProgressIndex = i;

	if (i > 0)
		{
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(g_ProgressDesc, Format, ArgList);
		}

	string LevelStr;
	GetProgressLevelStr(LevelStr);
	Progress(" %s\r", LevelStr.c_str());

	if (IsLastStep)
		{
		g_CountsInterval = 1;
		fputc('\n', stderr);
		}
	}

enum OptType
	{
	OT_Flag,
	OT_Tog,
	OT_Int,
	OT_Uns,
	OT_Str,
	OT_Float,
	OT_Enum
	};

struct OptInfo
	{
	void *Value;
	bool *OptSet;
	string LongName;
	OptType Type;
	int iMin;
	int iMax;
	unsigned uMin;
	unsigned uMax;
	double dMin;
	double dMax;
	map<string, unsigned> EnumValues;

	bool bDefault;
	int iDefault;
	unsigned uDefault;
	double dDefault;
	string strDefault;

	string Help;

	bool operator<(const OptInfo &rhs) const
		{
		return LongName < rhs.LongName;
		}
	};

static set<OptInfo> g_Opts;

void Help()
	{
	printf("\n");

	void Usage();
	Usage();

	for (set<OptInfo>::const_iterator p = g_Opts.begin(); p != g_Opts.end(); ++p)
		{
		const OptInfo &Opt = *p;

		printf("\n");
		string LongName = Opt.LongName.c_str();
		if (Opt.Type == OT_Tog)
			LongName = string("[no]") + LongName;
		printf("  --%s ", LongName.c_str());

		switch (Opt.Type)
			{
		case OT_Flag:
			break;
		case OT_Tog:
			break;
		case OT_Int:
			printf("<int>");
			break;
		case OT_Uns:
			printf("<uint>");
			break;
		case OT_Str:
			printf("<str>");
			break;
		case OT_Float:
			printf("<float>");
			break;
		case OT_Enum:
			printf("<enum>");
			break;
		default:
			printf("??type");
			break;
			}

		printf("  ");
		const string &s = Opt.Help;
		for (string::const_iterator q = s.begin(); q != s.end(); ++q)
			{
			char c = *q;
			if (c == '\n')
				printf("\n   ");
			else
				printf("%c", c);
			}
		printf("\n");
		}
	printf("\n");
	exit(0);
	}

void CmdLineErr(const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	string Str;
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);
	fprintf(stderr, "\n");
	fprintf(stderr, "Invalid command line\n");
	fprintf(stderr, "%s\n", Str.c_str());
	fprintf(stderr, "For list of command-line options use --help.\n");
	fprintf(stderr, "\n");
	exit(1);
	}

static set<OptInfo>::iterator GetOptInfo(const string &LongName,
  bool ErrIfNotFound)
	{
	for (set<OptInfo>::iterator p = g_Opts.begin();
	  p != g_Opts.end(); ++p)
		{
		const OptInfo &Opt = *p;
		if (Opt.LongName == LongName)
			return p;
		if (Opt.Type == OT_Tog && "no" + Opt.LongName == LongName)
			return p;
		}
	if (ErrIfNotFound)
		CmdLineErr("Option --%s is invalid", LongName.c_str());
	return g_Opts.end();
	}

static void AddOpt(const OptInfo &Opt)
	{
	if (GetOptInfo(Opt.LongName, false) != g_Opts.end())
		Die("Option --%s defined twice", Opt.LongName.c_str());
	g_Opts.insert(Opt);
	}

#ifdef _MSC_VER
#pragma warning(disable: 4505) // unreferenced local function
#endif

static void DefineFlagOpt(const string &LongName, const string &Help,
  void *Value, bool *OptSet)
	{
	*(bool *) Value = false;

	OptInfo Opt;
	Opt.Value = Value;
	Opt.OptSet = OptSet;
	Opt.LongName = LongName;
	Opt.bDefault = false;
	Opt.Help = Help;
	Opt.Type = OT_Flag;
	AddOpt(Opt);
	}

static void DefineTogOpt(const string &LongName, bool Default, const string &Help,
  void *Value, bool *OptSet)
	{
	*(bool *) Value = Default;

	OptInfo Opt;
	Opt.Value = Value;
	Opt.OptSet = OptSet;
	Opt.LongName = LongName;
	Opt.bDefault = Default;
	Opt.Help = Help;
	Opt.Type = OT_Tog;
	AddOpt(Opt);
	}

static void DefineIntOpt(const string &LongName, int Default, int Min, int Max,
  const string &Help, void *Value, bool *OptSet)
	{
	*(int *) Value = Default;

	OptInfo Opt;
	Opt.Value = Value;
	Opt.OptSet = OptSet;
	Opt.LongName = LongName;
	Opt.iDefault = Default;
	Opt.iMin = Min;
	Opt.iMax = Max;
	Opt.Help = Help;
	Opt.Type = OT_Int;
	AddOpt(Opt);
	}

static void DefineUnsOpt(const string &LongName, unsigned Default, unsigned Min,
  unsigned Max, const string &Help, void *Value, bool *OptSet)
	{
	*(unsigned *) Value = Default;

	OptInfo Opt;
	Opt.Value = Value;
	Opt.OptSet = OptSet;
	Opt.LongName = LongName;
	Opt.uDefault = Default;
	Opt.uMin = Min;
	Opt.uMax = Max;
	Opt.Help = Help;
	Opt.Type = OT_Uns;
	AddOpt(Opt);
	}

static void DefineFloatOpt(const string &LongName, double Default, double Min,
  double Max, const string &Help, void *Value, bool *OptSet)
	{
	*(double *) Value = Default;

	OptInfo Opt;
	Opt.Value = Value;
	Opt.OptSet = OptSet;
	Opt.LongName = LongName;
	Opt.dDefault = Default;
	Opt.dMin = Min;
	Opt.dMax = Max;
	Opt.Help = Help;
	Opt.Type = OT_Float;
	AddOpt(Opt);
	}

static void DefineStrOpt(const string &LongName, const char *Default,
  const string &Help, void *Value, bool *OptSet)
	{
	*(string *) Value = (Default == 0 ? "" : string(Default));

	OptInfo Opt;
	Opt.Value = Value;
	Opt.OptSet = OptSet;
	Opt.LongName = LongName;
	Opt.strDefault = (Default == 0 ? "" : string(Default));
	Opt.Help = Help;
	Opt.Type = OT_Str;
	AddOpt(Opt);
	}

static void ParseEnumValues(const string &Values, map<string, unsigned> &EnumValues)
	{
	EnumValues.clear();
	
	string Name;
	string Value;
	bool Eq = false;
	for (string::const_iterator p = Values.begin(); ; ++p)
		{
		char c = (p == Values.end() ? '|' : *p);
		if (isspace(c))
			;
		else if (c == '|')
			{
			if (EnumValues.find(Name) != EnumValues.end())
				Die("Invalid enum values, '%s' defined twice: '%s'",
				  Name.c_str(), Values.c_str());
			if (Name.empty() || Value.empty())
				Die("Invalid enum values, empty name or value: '%s'",
				  Values.c_str());

			EnumValues[Name] = atoi(Value.c_str());
			Name.clear();
			Value.clear();
			Eq = false;
			}
		else if (c == '=')
			Eq = true;
		else if (Eq)
			Value.push_back(c);
		else
			Name.push_back(c);
		if (p == Values.end())
			return;
		}
	}

static void DefineEnumOpt(const string &LongName, const string &ShortName,
  int Default, const string &Values, const string &Help, void *Value)
	{
	*(int *) Value = Default;

	OptInfo Opt;
	Opt.Value = Value;
	Opt.LongName = LongName;
	Opt.iDefault = Default;
	Opt.Help = Help;
	Opt.Type = OT_Enum;
	ParseEnumValues(Values, Opt.EnumValues);
	AddOpt(Opt);
	}
#undef FLAG_OPT
#undef TOG_OPT
#undef INT_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
#undef ENUM_OPT
#define FLAG_OPT(LongName)							bool opt_##LongName; bool optset_##LongName;
#define TOG_OPT(LongName, Default)					bool opt_##LongName; bool optset_##LongName;
#define INT_OPT(LongName, Default, Min, Max)		int opt_##LongName; bool optset_##LongName;
#define UNS_OPT(LongName, Default, Min, Max)		unsigned opt_##LongName; bool optset_##LongName;
#define FLT_OPT(LongName, Default, Min, Max)		double opt_##LongName; bool optset_##LongName;
#define STR_OPT(LongName, Default)					string opt_##LongName; bool optset_##LongName;
#define ENUM_OPT(LongName, Values, Default)			int opt_##LongName; bool optset_##LongName;
#include "myopts.h"

static int EnumStrToInt(const OptInfo &Opt, const string &Value)
	{
	const map<string, unsigned> &e = Opt.EnumValues;
	string s;
	for (map<string, unsigned>::const_iterator p = e.begin(); p != e.end(); ++p)
		{
		if (Value == p->first)
			return p->second;
		s += " " + p->first;
		}
	CmdLineErr("--%s %s not recognized, valid are: %s",
	  Opt.LongName.c_str(), Value.c_str(), s.c_str());
	ureturn(-1);
	}

static void SetOpt(OptInfo &Opt, const string &Value)
	{
	*Opt.OptSet = true;
	switch (Opt.Type)
		{
	case OT_Int:
		{
		*(int *) Opt.Value = atoi(Value.c_str());
		break;
		}
	case OT_Uns:
		{
		unsigned uValue = 0;
		int n = sscanf(Value.c_str(), "%u", &uValue);
		if (n != 1)
			CmdLineErr("Invalid value '%s' for --%s",
			  Value.c_str(), Opt.LongName.c_str());
		*(unsigned *) Opt.Value = uValue;
		break;
		}
	case OT_Float:
		{
		*(double *) Opt.Value = atof(Value.c_str());
		break;
		}
	case OT_Str:
		{
		*(string *) Opt.Value = Value;
		break;
		}
	case OT_Enum:
		{
		*(int *) Opt.Value = EnumStrToInt(Opt, Value);
		break;
		}
	default:
		asserta(false);
		}
	}

void LogOpts()
	{
	for (set<OptInfo>::const_iterator p = g_Opts.begin(); p != g_Opts.end(); ++p)
		{
		const OptInfo &Opt = *p;
		Log("%s = ", Opt.LongName.c_str());
		switch (Opt.Type)
			{
		case OT_Flag:
			Log("%s", (*(bool *) Opt.Value) ? "yes" : "no");
			break;
		case OT_Tog:
			Log("%s", (*(bool *) Opt.Value) ? "on" : "off");
			break;
		case OT_Int:
			Log("%d", *(int *) Opt.Value);
			break;
		case OT_Uns:
			Log("%u", *(unsigned *) Opt.Value);
			break;
		case OT_Float:
			{
			double Value = *(double *) Opt.Value;
			if (Value == FLT_MAX)
				Log("*");
			else
				Log("%g", Value);
			break;
			}
		case OT_Str:
			Log("%s", (*(string *) Opt.Value).c_str());
			break;
		case OT_Enum:
			Log("%d", *(int *) Opt.Value);
			break;
		default:
			asserta(false);
			}
		Log("\n");
		}
	}

static void CompilerInfo()
	{
#ifdef _FILE_OFFSET_BITS
    printf("_FILE_OFFSET_BITS=%d\n", _FILE_OFFSET_BITS);
#else
    printf("_FILE_OFFSET_BITS not defined\n");
#endif

#define x(t)	printf("sizeof(" #t ") = %d\n", (int) sizeof(t));
	x(int)
	x(long)
	x(float)
	x(double)
	x(void *)
	x(off_t)
#undef x
	exit(0);
	}

void Split(const string &Str, vector<string> &Fields, char Sep)
	{
	Fields.clear();
	const unsigned Length = (unsigned) Str.size();
	string s;
	for (unsigned i = 0; i < Length; ++i)
		{
		char c = Str[i];
		if ((Sep == 0 && isspace(c)) || c == Sep)
			{
			if (!s.empty() || Sep != 0)
				Fields.push_back(s);
			s.clear();
			}
		else
			s.push_back(c);
		}
	if (!s.empty())
		Fields.push_back(s);
	}

static void GetArgsFromFile(const string &FileName, vector<string> &Args)
	{
	Args.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		size_t n = Line.find('#');
		if (n != string::npos)
			Line = Line.substr(0, n);
		vector<string> Fields;
		Split(Line, Fields);
		Args.insert(Args.end(), Fields.begin(), Fields.end());
		}
	CloseStdioFile(f);
	}

void MyCmdLine(int argc, char **argv)
	{
	static unsigned RecurseDepth = 0;
	++RecurseDepth;

	DefineFlagOpt("compilerinfo", "Write info about compiler types and #defines to stdout.",
	  (void *) &opt_compilerinfo, &optset_compilerinfo);
	DefineFlagOpt("quiet", "Turn off progress messages.", (void *) &opt_quiet, &optset_quiet);
	DefineFlagOpt("version", "Show version and exit.", (void *) &opt_version, &optset_version);
	DefineFlagOpt("logopts", "Log options.", (void *) &opt_logopts, &optset_logopts);
	DefineFlagOpt("help", "Display command-line options.", (void *) &opt_help, &optset_help);
	DefineStrOpt("log", "", "Log file name.", (void *) &opt_log, &optset_log);

#undef FLAG_OPT
#undef TOG_OPT
#undef INT_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
#undef ENUM_OPT
#define FLAG_OPT(LongName)						DefineFlagOpt(#LongName, "help", (void *) &opt_##LongName, &optset_##LongName);
#define TOG_OPT(LongName, Default)				DefineTogOpt(#LongName, Default, "help", (void *) &opt_##LongName, &optset_##LongName);
#define INT_OPT(LongName, Default, Min, Max)	DefineIntOpt(#LongName, Default, Min, Max, "help", (void *) &opt_##LongName, &optset_##LongName);
#define UNS_OPT(LongName, Default, Min, Max)	DefineUnsOpt(#LongName, Default, Min, Max, "help", (void *) &opt_##LongName, &optset_##LongName);
#define FLT_OPT(LongName, Default, Min, Max)	DefineFloatOpt(#LongName, Default, Min, Max, "help", (void *) &opt_##LongName, &optset_##LongName);
#define STR_OPT(LongName, Default)				DefineStrOpt(#LongName, Default, "help", (void *) &opt_##LongName, &optset_##LongName);
#define ENUM_OPT(LongName, Values, Default)		DefineEnumOpt(#LongName, Values, Default, "help", (void *) &opt_##LongName, &optset_##LongName);
#include "myopts.h"

	if (RecurseDepth == 0)
		g_Argv.clear();

	for (int i = 0; i < argc; ++i) 
		g_Argv.push_back(string(argv[i]));
	

	int i = 1;
	for (;;)
		{
		if (i >= argc)
			break;
		const string &Arg = g_Argv[i];
		
		if (Arg.empty())
			continue;
		else if (Arg == "file:" && i + 1 < argc)
			{
			const string &FileName = g_Argv[i+1];
			vector<string> Args;
			GetArgsFromFile(FileName, Args);
			for (vector<string>::const_iterator p = Args.begin();
			  p != Args.end(); ++p)
				{
				g_Argv.push_back(*p);
				++argc;
				}
			i += 2;
			continue;
			}
		else if (Arg.size() > 1 && Arg[0] == '-')
			{
			string LongName = (Arg.size() > 2 && Arg[1] == '-' ? Arg.substr(2) : Arg.substr(1));
			OptInfo Opt = *GetOptInfo(LongName, true);
			*Opt.OptSet = true;
			if (Opt.Type == OT_Flag)
				{
				g_Opts.erase(Opt);
				*(bool *) Opt.Value = true;
				g_Opts.insert(Opt);
				++i;
				continue;
				}
			else if (Opt.Type == OT_Tog)
				{
				g_Opts.erase(Opt);
				if (string("no") + Opt.LongName == LongName)
					*(bool *) Opt.Value = false;
				else
					{
					asserta(Opt.LongName == LongName);
					*(bool *) Opt.Value = true;
					}
				g_Opts.insert(Opt);
				++i;
				continue;
				}

			++i;
			if (i >= argc)
				CmdLineErr("Missing value for option --%s", LongName.c_str());

			string Value = g_Argv[i];
			SetOpt(Opt, Value);

			++i;
			continue;
			}
		else
			CmdLineErr("Expected -option_name or --option_name, got '%s'", Arg.c_str());
		}

	--RecurseDepth;
	if (RecurseDepth > 0)
		return;

	if (opt_help)
		Help();

	if (opt_compilerinfo)
		CompilerInfo();

	SetLogFileName(opt_log);

	if (opt_log != "")
		{
		for (int i = 0; i < argc; ++i)
			Log("%s%s", i == 0 ? "" : " ", g_Argv[i].c_str());
		Log("\n");
		time_t Now = time(0);
		struct tm *t = localtime(&Now);
		const char *s = asctime(t);
		Log("Started %s", s); // there is a newline in s
		Log("Version " MY_VERSION ".%s\n", SVN_VERSION);
		Log("\n");
		}

	if (opt_logopts)
		LogOpts();
	}

double Pct(double x, double y)
	{
	if (y == 0.0f)
		return 0.0f;
	return (x*100.0f)/y;
	}

void GetCmdLine(string &s)
	{
	s.clear();
	for (unsigned i = 0; i < SIZE(g_Argv); ++i)
		{
		if (i > 0)
			s += " ";
		s += g_Argv[i];
		}
	}

char *mystrsave(const char *s)
	{
	unsigned n = unsigned(strlen(s));
	char *t = myalloc(char, n+1);
	memcpy(t, s, n+1);
	return t;
	}

void Logu(unsigned u, unsigned w, unsigned prefixspaces)
	{
	for (unsigned i = 0; i < prefixspaces; ++i)
		Log(" ");
	if (u == UINT_MAX)
		Log("%*.*s", w, w, "*");
	else
		Log("%*u", w, u);
	}

void Logf(float x, unsigned w, unsigned prefixspaces)
	{
	for (unsigned i = 0; i < prefixspaces; ++i)
		Log(" ");
	if (x == FLT_MAX)
		Log("%*.*s", w, w, "*");
	else
		Log("%*.2f", w, x);
	}

static uint32 g_SLCG_state = 1;

// Numerical values used by Microsoft C, according to wikipedia:
// http://en.wikipedia.org/wiki/Linear_congruential_generator
static uint32 g_SLCG_a = 214013;
static uint32 g_SLCG_c = 2531011;

// Simple Linear Congruential Generator
// Bad properties; used just to initialize the better generator.
static uint32 SLCG_rand()
	{
	g_SLCG_state = g_SLCG_state*g_SLCG_a + g_SLCG_c;
	return g_SLCG_state;
	}

static void SLCG_srand(uint32 Seed)
	{
	g_SLCG_state = Seed;
	for (int i = 0; i < 10; ++i)
		SLCG_rand();
	}

/***
A multiply-with-carry random number generator, see:
http://en.wikipedia.org/wiki/Multiply-with-carry

The particular multipliers used here were found on
the web where they are attributed to George Marsaglia.
***/

static bool g_InitRandDone = false;
static uint32 g_X[5];

uint32 RandInt32()
	{
	InitRand();

	uint64 Sum = 2111111111*(uint64) g_X[3] + 1492*(uint64) g_X[2] +
	  1776*(uint64) g_X[1] + 5115*(uint64) g_X[0] + g_X[4];
	g_X[3] = g_X[2];
	g_X[2] = g_X[1];
	g_X[1] = g_X[0];
	g_X[4] = (uint32) (Sum >> 32);
	g_X[0] = (uint32) Sum;
	return g_X[0];
	}

unsigned randu32()
	{
	return (unsigned) RandInt32();
	}

void InitRand()
	{
	if (g_InitRandDone)
		return;
// Do this first to avoid recursion
	g_InitRandDone = true;

	unsigned Seed = (optset_randseed ? opt_randseed : (unsigned) (time(0)*getpid()));
	Log("RandSeed=%u\n", Seed);
	SLCG_srand(Seed);

	for (unsigned i = 0; i < 5; i++)
		g_X[i] = SLCG_rand();

	for (unsigned i = 0; i < 100; i++)
		RandInt32();
	}

// MUST COME AT END BECAUSE OF #undef
#if	RCE_MALLOC
#undef mymalloc
#undef myfree
#undef myfree2
void *mymalloc(unsigned bytes, const char *FileName, int Line)
	{
	void *rce_malloc(unsigned bytes, const char *FileName, int Line);
	return rce_malloc(bytes, FileName, Line);
	}

void myfree(void *p, const char *FileName, int Line)
	{
	void rce_free(void *p, const char *FileName, int Line);
	rce_free(p, FileName, Line);
	}

void myfree2(void *p, unsigned bytes, const char *FileName, int Line)
	{
	void rce_free(void *p, const char *FileName, int Line);
	rce_free(p, FileName, Line);
	}

#else // RCE_MALLOC
void *mymalloc(unsigned bytes)
	{
	++g_NewCalls;
	if (g_InitialMemUseBytes == 0)
		g_InitialMemUseBytes = GetMemUseBytes();

	g_TotalAllocBytes += bytes;
	g_NetBytes += bytes;
	if (g_NetBytes > g_MaxNetBytes)
		{
		if (g_NetBytes > g_MaxNetBytes + 10000000)
			GetMemUseBytes();//to force update of peak
		g_MaxNetBytes = g_NetBytes;
		}
	void *p = malloc(bytes);
	//void *p = _malloc_dbg(bytes, _NORMAL_BLOCK, __FILE__, __LINE__);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\nOut of memory mymalloc(%u), curr %.3g bytes",
		  (unsigned) bytes, b);
		void LogAllocs();
		LogAllocs();
#if DEBUG && defined(_MSC_VER)
		asserta(_CrtCheckMemory());
#endif
		Die("Out of memory, mymalloc(%u), curr %.3g bytes\n",
		  (unsigned) bytes, b);
		}
	return p;
	}

void myfree(void *p)
	{
	if (p == 0)
		return;
	free(p);
	//_free_dbg(p, _NORMAL_BLOCK);
	}

void myfree2(void *p, unsigned bytes)
	{
	++g_FreeCalls;
	g_TotalFreeBytes += bytes;
	g_NetBytes -= bytes;

	if (p == 0)
		return;
	free(p);
	}
#endif
