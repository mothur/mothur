/*
 * dayhoff.h
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
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_NJ_DAYHOFF_H_
#define _INC_NJ_DAYHOFF_H_ 1

/*
 * As sequence divergence increases, we need to correct for multiple hits
 * when using Kimura distance correction method for amino acid sequences.
 *
 * This matrix of values represents the estimated "Accepted Point Mutations"
 * or PAMs for a range of amino acid sequence divergence, starting at 75% 
 * up through 93% (in 0.1% increments).
 *
 * This model is derived from Dayhoff (1978).
 *
 * This Dayhoff matrix and the shortcut methods for dealing with Kimura
 * correction at high sequence divergence (> 75%) are derived from similar 
 * work in Clustal W:
 *
 *     Thompson, J.D., Higgins, D.G., Gibson, T.J., "CLUSTAL W:
 *     improving the sensitivity of progressive multiple sequence
 *     alignment through sequence weighting, position-specific gap
 *     penalties and weight matrix choice.", 
 *     Nucleic Acids Research, 22:4673-4680, 1994
 *
 */


int NJ_dayhoff[]={
  195,    196,    197,    198,    199,    200,    200,    201,    202,  203,
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,
  246,    248,    249,    250,    252,    253,    254,    255,    257,  258,
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
  988 
};



#endif /* _INC_NJ_DAYHOFF_H_ */



