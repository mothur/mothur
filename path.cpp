//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#include "myutils.h"
#include "path.h"
#include "timing.h"

#define TRACE	0

const unsigned PathMagic = 0x9A783A16;

struct PathBuffer
	{
	unsigned Magic;
	char *Buffer;
	unsigned Size;
	bool InUse;
	};

static PathBuffer **g_PathBuffers;
static unsigned g_PathBufferSize;

static char *AllocBuffer(unsigned Size)
	{
	if (Size == 0)
		return 0;

// Is a free buffer that is big enough?
	for (unsigned i = 0; i < g_PathBufferSize; ++i)
		{
		PathBuffer *PB = g_PathBuffers[i];
		asserta(PB->Magic == PathMagic);
		if (!PB->InUse)
			{
			if (PB->Size >= Size)
				{
				PB->InUse = true;
				return PB->Buffer;
				}
			if (PB->Buffer == 0)
				{
				unsigned Size2 = Size + 1024;
				PB->Buffer = MYALLOC(char, Size2, Path);
				PB->Size = Size2;
				PB->InUse = true;
				return PB->Buffer;
				}
			}
		}

// No available buffer, must expand g_PathBuffers[]
	unsigned NewPathBufferSize = g_PathBufferSize + 1024;
	PathBuffer **NewPathBuffers = MYALLOC(PathBuffer *, NewPathBufferSize, Path);
	
	for (unsigned i = 0; i < g_PathBufferSize; ++i)
		NewPathBuffers[i] = g_PathBuffers[i];

	for (unsigned i = g_PathBufferSize; i < NewPathBufferSize; ++i)
		{
		PathBuffer *PB = MYALLOC(PathBuffer, 1, Path);
		PB->Magic = PathMagic;
		PB->Buffer = 0;
		PB->Size = 0;
		PB->InUse = false;
		NewPathBuffers[i] = PB;
		}

	PathBuffer *PB = NewPathBuffers[g_PathBufferSize];

	MYFREE(g_PathBuffers, g_PathBufferSize, Path);
	g_PathBuffers = NewPathBuffers;
	g_PathBufferSize = NewPathBufferSize;

	asserta(!PB->InUse && PB->Buffer == 0);

	unsigned Size2 = Size + 1024;
	PB->Buffer = MYALLOC(char, Size2, Path);
	PB->Size = Size2;
	PB->InUse = true;
	return PB->Buffer;
	}

static void FreeBuffer(char *Buffer)
	{
	if (Buffer == 0)
		return;

	for (unsigned i = 0; i < g_PathBufferSize; ++i)
		{
		PathBuffer *PB = g_PathBuffers[i];
		if (PB->Buffer == Buffer)
			{
			asserta(PB->InUse);
			PB->InUse = false;
			return;
			}
		}

	Die("FreeBuffer, not found");
	}

void PathData::Alloc(unsigned MaxLen)
	{
	if (MaxLen < Bytes)
		return;

	StartTimer(PathAlloc);
	if (Bytes > 0)
		{
		FreeBuffer(Front);
		}

	Bytes = MaxLen + 1;
	Front = AllocBuffer(Bytes);
	Back = Front + Bytes - 1;
	Start = 0;
	EndTimer(PathAlloc);
	}

void PathData::Free()
	{
	FreeBuffer(Front);
	Front = 0;
	Start = 0;
	Back = 0;
	}

void PathData::Copy(const PathData &rhs)
	{
	Alloc(rhs.Bytes);
	strcpy(Front, rhs.Front);
	Start = Front + (rhs.Start - rhs.Front);
	}

void PathData::FromStr(const char *PathStr)
	{
	asserta(PathStr != 0);
	unsigned NeededBytes = (unsigned) strlen(PathStr) + 1;
	Alloc(NeededBytes);
	strcpy(Front, PathStr);
	Start = Front;
	}

void LogPathStats()
	{
	Log("\n");
	unsigned Bytes = 0;
	for (unsigned i = 0; i < g_PathBufferSize; ++i)
		{
		const PathBuffer *PB = g_PathBuffers[i];
		Bytes += PB->Size;
		}
	Log("%u paths allocated, total memory %u bytes\n", g_PathBufferSize, Bytes);
	}
