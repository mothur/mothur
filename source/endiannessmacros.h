#ifndef EDIANNESSMACROS_H
#define EDIANNESSMACROS_H

/*
 *   endiannessmacros.h
 *  Mothur
 *
 *  Created by westcott on 7/9/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
/*********************************************************************/
/*********************************************************************/
// The following is copied from the staden io_lib-1.12.4 os.h - thanks!
/*********************************************************************/
/*********************************************************************/

/*
 * Author: 
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: operating system specific type definitions
 *
 */


/* Mac FAT binaries or unknown. Auto detect based on CPU type */
#if !defined(SP_BIG_ENDIAN) && !defined(SP_LITTLE_ENDIAN)
	
/*
 * x86 equivalents
 */
#if  defined(__i386) || defined(__i386__) || defined(__ia64__) ||  defined(WIN32) || defined(__arm__) || (defined(__mips__) && defined(__MIPSEL__)) || defined(__SYMBIAN32__) || \
     defined(__x86_64__) || defined(__x86_64) || defined(__i686__) || defined(__i686) || defined(__amd64__) || defined(__amd64) || defined(__LITTLE_ENDIAN__)
#define SP_LITTLE_ENDIAN
#else
#define SP_BIG_ENDIAN
#endif



/*
 * SUN Sparc
 */
#if defined(__sparc__) || defined(__sparc)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/* Some catch-alls */
#if defined(__LITTLE_ENDIAN__) || defined(__LITTLEENDIAN__)
#    define SP_LITTLE_ENDIAN
#endif

#if defined(__BIG_ENDIAN__) || defined(__BIGENDIAN__)
#    define SP_BIG_ENDIAN
#endif

#if defined(SP_BIG_ENDIAN) && defined(SP_LITTLE_ENDIAN)
#    error Both BIG and LITTLE endian defined. Fix os.h and/or Makefile
#endif

#if !defined(SP_BIG_ENDIAN) && !defined(SP_LITTLE_ENDIAN)
#    error Neither BIG nor LITTLE endian defined. Fix os.h and/or Makefile
#endif

#endif

/*-----------------------------------------------------------------------------
 * Byte swapping macros
 */

/*
 * Our new swap runs at the same speed on Ultrix, but substantially faster
 * (300% for swap_int4, ~50% for swap_int2) on an Alpha (due to the lack of
 * decent 'char' support).
 *
 * They also have the ability to swap in situ (src == dst). Newer code now
 * relies on this so don't change back!
 */
#define iswap_int8(x) \
    (((x & 0x00000000000000ffLL) << 56) + \
     ((x & 0x000000000000ff00LL) << 40) + \
     ((x & 0x0000000000ff0000LL) << 24) + \
     ((x & 0x00000000ff000000LL) <<  8) + \
     ((x & 0x000000ff00000000LL) >>  8) + \
     ((x & 0x0000ff0000000000LL) >> 24) + \
     ((x & 0x00ff000000000000LL) >> 40) + \
     ((x & 0xff00000000000000LL) >> 56))

#define iswap_int4(x) \
    (((x & 0x000000ff) << 24) + \
     ((x & 0x0000ff00) <<  8) + \
     ((x & 0x00ff0000) >>  8) + \
     ((x & 0xff000000) >> 24))

#define iswap_int2(x) \
    (((x & 0x00ff) << 8) + \
     ((x & 0xff00) >> 8))

#define swap_int8(src, dst) ((dst) = iswap_int8(src))
#define swap_int4(src, dst) ((dst) = iswap_int4(src))
#define swap_int2(src, dst) ((dst) = iswap_int2(src))


/*
 * Linux systems may use byteswap.h to get assembly versions of byte-swap
 * on intel systems. This can be as trivial as the bswap opcode, which works
 * out at over 2-times faster than iswap_int4 above.
 */
#if 0
#if defined(__linux__)
#    include <byteswap.h>
#    undef iswap_int8
#    undef iswap_int4
#    undef iswap_int2
#    define iswap_int8 bswap_64
#    define iswap_int4 bswap_32
#    define iswap_int2 bswap_16
#endif
#endif


/*
 * Macros to specify that data read in is of a particular endianness.
 * The macros here swap to the appropriate order for the particular machine
 * running the macro and return the new answer. These may also be used when
 * writing to a file to specify that we wish to write in (eg) big endian
 * format.
 *
 * This leads to efficient code as most of the time these macros are
 * trivial.
 */
#ifdef SP_BIG_ENDIAN
#define be_int8(x) (x)
#define be_int4(x) (x)
#define be_int2(x) (x)
#define be_int1(x) (x)

#define le_int8(x) iswap_int8((x))
#define le_int4(x) iswap_int4((x))
#define le_int2(x) iswap_int2((x))
#define le_int1(x) (x)
#endif

#ifdef SP_LITTLE_ENDIAN
#define be_int8(x) iswap_int8((x))
#define be_int4(x) iswap_int4((x))
#define be_int2(x) iswap_int2((x))
#define be_int1(x) (x)

#define le_int8(x) (x)
#define le_int4(x) (x)
#define le_int2(x) (x)
#define le_int1(x) (x)
#endif

#endif

