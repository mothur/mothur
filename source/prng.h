/*
 * prng.h
 *
 * $Id$
 *
 *****************************************************************************
 *
 * Copyright (c) 2004,  Luke Sheneman
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions 
 * are met:
 * 
 *  + Redistributions of source code must retain the above copyright 
 *    notice, this list of conditions and the following disclaimer. 
 *  + Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in 
 *    the documentation and/or other materials provided with the 
 *    distribution. 
 *  + The names of its contributors may not be used to endorse or promote 
 *    products derived  from this software without specific prior 
 *    written permission. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.  
 *
 *****************************************************************************
 *
 * Some function prototypes for the Mersenne Twister PRNG
 *
 *****************************************************************************
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_PRNG_H_
#define _INC_PRNG_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

#define NJ_RAND_MAX 0x7fffffffUL


/* some function prototypes */
void
init_genrand(unsigned long s);

void 
init_by_array(unsigned long init_key[],
	      int key_length);

unsigned long 
genrand_int32(void);

long int 
genrand_int31(void);

double 
genrand_real1(void);

double
genrand_real2(void);

double
genrand_real3(void);

double
genrand_res53(void);

long int
NJ_genrand_int31_top(long int top);

#ifdef __cplusplus
}
#endif

#endif /* _INC_PRNG_H_ */




